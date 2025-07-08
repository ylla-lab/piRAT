import argparse
import sys
import os
import multiprocessing
import time
from pandas.errors import EmptyDataError
from typing import Tuple, List
import pirat.configuration_functions.configuration as conf
import pirat.configuration_functions.file_operations as f_op
import pirat.configuration_functions.user_input as u_in
import pirat.clustering_functions.data_preparation as data_prep
import pirat.clustering_functions.parameters as param
import pirat.clustering_functions.clustering as clust
import pirat.statistics_functions.statistics as statistics
import pirat.ping_pong.ping_pong_all as pin_pon
import pirat.configuration_functions.post_processing as pp
import pirat.configuration_functions.reports_generation as rg
from pirat._version import __version__

def clustering(raw_path_name: str,
               range_of_size: Tuple[int, int],
               threads: int,
               draw_plots: bool,
               s: int,
               k: int,
               eps: int,
               variation_threshold: int,
               which_type: str,
               sample_wise_analysis_files: List[str],
               plot_iter: bool,
               success: bool,
               saving_directory: str,
               automate: bool) -> Tuple[int, int]:
    """
        Perform the full pipeline of clustering

        Parameters:
            range_of_size: Tuple of ints of min and max length of piRNAs
            list_of_scaffolds: A list of scaffolds, where each scaffold is represented as a list of reads or data
                                            needed for clustering.
            k: Minimum number of reads within a specified distance (eps) to form a cluster.
            eps: Maximum allowed distance between reads for them to be considered part of the same cluster.
            range_of_size: A list specifying the minimum and maximum sizes of piRNAs ([min_size, max_size]).
            raw_path_name: The path name to the raw data directory or file used as input to the clustering function.
            threads: Number of CPU threads to use for parallel processing.
            variation_threshold: Threshold value used to filter or determine acceptable variability within clusters.
            saving_path_primary: Directory or file path for saving the primary clustering results.

        Returns:
            List: A list of results from the clustering process. Each element corresponds to the output of
                  `run_the_clustering` for a particular scaffold.

        Example:
            >>> multiprocess_clustering([["scaffold_1_data"], ["scaffold_2_data"]], k=10, eps=1000, range_of_size=[20,30],
            ...                                   raw_path_name="/path/to/raw", threads=4, variation_threshold=5,
            ...                                   saving_path_primary="/path/to/save")
            [["cluster_results_for_scaffold_1"], ["cluster_results_for_scaffold_2"]]

        Notes:
            - This function relies on `run_the_clustering` to process each scaffold.
            - Multiprocessing reduces computation time for large datasets by leveraging multiple CPU cores.
        """

    if not automate:
        if k is None or eps is None:
            autom_param = ''
            while autom_param not in ['y', 'n']:
                print("\nYou didn't specify the parameters, do you want to automatically find the parameters? [y/n]")
                autom_param = input()
                if autom_param == 'n':
                    print('\nPlease input your parameters as k,Eps. e.g: 8,700')
                    params = input()
                    k, eps = [int(x) for x in params.split(',')]
                elif autom_param == 'y':
                    k, eps = None, None
        else:
            k, eps = int(k), int(eps)

    list_of_files = data_prep.get_list_of_bam_files(raw_path_name)
    list_of_scaffolds: list[str] = data_prep.get_list_of_scaffolds(raw_path_name + list_of_files[0])
    if range_of_size is None:
        if which_type == 1 or success is None:
            range_of_size = param.find_scope_of_the_size(list_of_files, raw_path_name, saving_directory)
            if not automate:
                range_of_size = u_in.define_scope_range(saving_directory, range_of_size)
        elif which_type == 2:
            sample_wise_analysis_files = ['_'.join((file.split('.'))[:-1]) for file in sample_wise_analysis_files]
            try:
                range_of_size = param.find_scope_with_sample_analysis(sample_wise_analysis_files, saving_directory)
            except ValueError:
                print("\npiRAT coudln't automatically find range of size of piRNAs in the sample! Please input the range manually.")
                range_of_size = u_in.define_scope_range(saving_directory, [0,0])
            if not automate:
                range_of_size = u_in.define_scope_range(saving_directory, range_of_size)
        if range_of_size[0] < 26 and automate is True:
            print("\nDetected range of piRNAs in the sample is shorter than theoretical length of piRNAs.")
            range_of_size = u_in.define_scope_range(saving_directory, range_of_size)
        else:
            print(f"\nRange of size of piRNAs in the sample: {range_of_size}\n")
    if k is None and eps is None:
        print('Finding optimal clustering parameters...\n')
        data = [(scaffold, list_of_files, raw_path_name, variation_threshold, range_of_size)\
                for scaffold in list_of_scaffolds]
        with multiprocessing.Pool(threads) as pool:
            results = pool.starmap(param.find_parameters, data)
        appriopriate_dict_to_find_params = {4: [], 6: [], 8: [], 10: []}
        for what_scaffold in range(len(results)):
            for what_k in results[what_scaffold].keys():
                [appriopriate_dict_to_find_params[what_k].append(x) for x in results[what_scaffold][what_k]]
        dic_of_params, parameter_quality = param.finding_the_parameters(appriopriate_dict_to_find_params, saving_directory)
        k, eps = param.get_best_pair_of_parameters(dic_of_params, parameter_quality)
        if k is False:
            k, eps = u_in.define_params(k, eps, dic_of_params, saving_directory)
        if not automate:
            k, eps = u_in.define_params(k, eps, dic_of_params, saving_directory)
            print(f"\nParameters k and Eps are set as: {k}\t{eps}\n")
    list_of_scaffolds.remove('*') if '*' in list_of_scaffolds else None
    if s is None:
        results = clust.multiprocess_clustering(list_of_scaffolds, k, eps, range_of_size, raw_path_name, threads,
                                                variation_threshold, saving_directory)
    else:
        results = clust.multiprocess_clustering(list_of_scaffolds[:1], k, eps, range_of_size, raw_path_name, threads,
                                                variation_threshold, saving_directory)
    number_of_clusters = sum([len(x[1]) for x in results])
    print(f"Found {number_of_clusters} clusters!\n")
    try:
        f_op.save_clusters_to_gff(results, saving_directory)
        statistics.transform_gff_file(saving_directory)
        statistics.classify_clusters(saving_directory, eps)
        list_of_clusters, list_of_clusters_of_reads = statistics.get_cluster_reads_out_of_gff(saving_directory)

        list_of_clusters_for_statistics, list_of_scfs, cluster_length_distribution = statistics.read_output_gff(
        f'{saving_directory}gffs/clusters_out.gff')
        high_quality_clusters, low_quality_clusters = statistics.run_the_statistics(list_of_clusters_for_statistics,
                                                                                range_of_size[1],
                                                                                threads, list_of_files, draw_plots,
                                                                                list_of_clusters,
                                                                                list_of_clusters_of_reads,
                                                                                range_of_size, raw_path_name, plot_iter,
                                                                                    saving_directory)
        f_op.save_high_low_quality_to_gff([high_quality_clusters, low_quality_clusters], saving_directory)
        print(f"\n{len(high_quality_clusters)} clusters out of {number_of_clusters} found clusters are high quality!\n")
        statistics.statistics_of_found_clusters(raw_path_name, high_quality_clusters, low_quality_clusters, range_of_size, saving_directory)
        pin_pon.get_matrix_of_seqs([], range_of_size[1], f'{saving_directory}', f'{saving_directory}plots/', \
                                   f'{saving_directory}gffs/Clusters_reads.gff')
    except EmptyDataError:
        print("\nNo clusters have been found!")
    return k, eps, threads, variation_threshold, range_of_size

def sample_wise_analysis(raw_path_name: str,
                         variation_threshold: int,
                         list_of_files: List[str],
                         flag: int,
                         saving_path: str,
                         threads: int,
                         paired: bool,
                         current_time: str,
                         __version__: str) -> bool:
    path = raw_path_name
    pirna_sequences, overlapping_lengths, length_read_distribution, pirna_dict = \
        pin_pon.load_data(path, list_of_files, variation_threshold, [20, 40], flag, threads, paired)
    pirna_lengths = [len(x[2]) for x in pirna_sequences]
    file = list_of_files[0]
    bam_file_name = '_'.join((file.split('.'))[:-1]) if len(list_of_files) == 1 else ''
    saving_report_path = f"{saving_path}/reports/"
    saving_path = f"{saving_path}sample_wise_analysis/{bam_file_name}/"
    try:
        pin_pon.draw_plots(f"{saving_path}",
                           [["Initial_distribution", length_read_distribution, "Read length (nt)",
                             "Number of reads"],
                            ["5'-to-5' overlap", overlapping_lengths, "Length of 5'-to-5' overlap (nt)",
                             "Number of sequence pairs"],
                            ["Length distribution of the sequences\ndisplaying 10nts 5'-to-5' overlaps", pirna_lengths, "Length of the sequence (nt)",
                             "Number of sequences"]], saving_path)
        data_for_heatmap = pin_pon.prepare_data_for_heatmap(pirna_dict, (20, 40))
        pin_pon.generate_heatmap(saving_path, data_for_heatmap, [20, 40], True, 3)
        pin_pon.get_matrix_of_seqs([x[2] for x in pirna_sequences], 32, saving_path, saving_path)
        rg.save_report(saving_path, bam_file_name, current_time, raw_path_name, __version__, list_of_files, saving_report_path)
        return True
    except ValueError:
        print(f"There aren't any ping-pong signatures in the file {bam_file_name}!")
        rg.save_report_no_ping_pong(saving_path, bam_file_name, current_time, raw_path_name, __version__)
        return None

def ping_pong_annotate(raw_path_name: str,
                       range_of_size: Tuple[int, int],
                       variation_threshold: int,
                       list_of_files: List[str],
                       flag: int,
                       sample_wise_analysis_files: List[str],
                       saving_directory: str,
                       name: str,
                       threads: int,
                       automate: bool,
                       paired: bool,
                       current_time: str,
                       __version__: str) -> None:
    path = raw_path_name[:-1]
    saving_path = f"{saving_directory}"
    saving_path_plots = f"{saving_directory}/plots/"
    if range_of_size is None:
        try:
            sample_wise_analysis_files = ['_'.join((file.split('.'))[:-1]) for file in sample_wise_analysis_files]
            range_of_size = param.find_scope_with_sample_analysis(sample_wise_analysis_files, saving_directory)
        except ValueError:
            range_of_size = param.find_scope_of_the_size(list_of_files, raw_path_name, saving_directory)
        if automate is False:
            range_of_size = u_in.define_scope_range(raw_path_name, range_of_size)
            print(f"\nRange of size is set as: {range_of_size}\n")
        else:
            pass
    pirna_sequences, overlapping_lengths, length_read_distribution, pirna_dict = pin_pon.load_data(path, list_of_files,
                                                                                                   variation_threshold,
                                                                                                   range_of_size, flag, threads, paired)
    pin_pon.save_to_gff(saving_directory, pirna_sequences)
    try:
        data_for_2nd_heatmap = pin_pon.generete_venn(saving_directory + 'gffs/potential_pirna_reads.gff',
                                          saving_directory + 'gffs/ping_pong_reads.gff', saving_directory + 'gffs/clusters_out.gff', saving_path_plots, True)
        data_for_heatmap = pin_pon.prepare_data_for_heatmap(pirna_dict, (range_of_size[0], range_of_size[1]))
        pin_pon.generate_heatmap(saving_path_plots, data_for_heatmap, [range_of_size[0], range_of_size[1]], True, 1)
        pin_pon.generate_heatmap(saving_path_plots, data_for_2nd_heatmap, [range_of_size[0], range_of_size[1]], True, 2)
        pin_pon.get_matrix_of_seqs([x[2] for x in pirna_sequences], range_of_size[1], saving_path, saving_path_plots, None)
        pin_pon.draw_plots_annotate(f"{saving_directory}",
                               [["Initial_distribution", 'Initial_distribution.txt', "Read length (nt)", "Number of reads"],
                                ["5'-to-5' overlap", "5-to-5_overlap.txt", "Length of 5'-to-5' overlap (nt)", "Number of sequence pairs"],
                                ["Length distribution of the sequences\ndisplaying 10nts 5'-to-5' overlaps", 'Length_distribution_of_the_sequences_displaying_10nts_5-to-5_overlaps.txt', "Length of the sequence (nt)",
                                 "Number of sequences"]], ['_'.join((x.split('.'))[:-1]) for x in list_of_files], saving_path_plots)


        rg.save_report_annotate(saving_directory, name, current_time, raw_path_name, __version__, list_of_files)
    except:
        data_for_heatmap, data_for_2nd_heatmap = {}, {}
        pin_pon.generate_heatmap(saving_path_plots, data_for_heatmap, [range_of_size[0], range_of_size[1]], True, 1)
        pin_pon.generate_heatmap(saving_path_plots, data_for_2nd_heatmap, [range_of_size[0], range_of_size[1]], False, 2)
        #pin_pon.get_matrix_of_seqs([x[2] for x in pirna_sequences], range_of_size[1], saving_path, saving_path_plots, None)
        pin_pon.draw_plots_annotate(f"{saving_path}",
                                    [["Initial_distribution", 'Initial_distribution.txt', "Read length (nt))",
                                      "Number of reads"],
                                     ["5'-to-5' overlap", "5-to-5_overlap.txt",
                                      "Length of 5'-to-5' overlap (nt)", "Number of sequence pairs"],
                                     ["Length distribution of the sequences\ndisplaying 10nts 5'-to-5' overlaps",
                                      'Length_distribution_of_the_sequences_displaying_10nts_5-to-5_overlaps.txt', "Length of the sequence (nt)",
                                      "Number of sequences"]], ['_'.join((x.split('.'))[:-1]) for x in list_of_files], saving_path_plots)
        try:
            rg.save_report_annotate_wo_primary(saving_directory, name, current_time, raw_path_name, __version__)
        except FileNotFoundError:
            rg.save_report_annotate_wo_primary_and_secondary(saving_directory, name, current_time, raw_path_name, __version__)

    return None


def everything() -> None:
    current_time = time.strftime("%a, %d %b %Y %X", time.localtime())
    parser = argparse.ArgumentParser(description="Input",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", metavar="--path", type=str, help="Path to the directory with data, e.g /data/")
    parser.add_argument("-o", metavar="--output-path", type=str, help="Path to the directory where the output should be located")
    parser.add_argument("-k", metavar="--minreads", type=int, help="MinReads parameter")
    parser.add_argument("-e", metavar="--eps", type=int, help="Eps parameter")
    parser.add_argument("-r", metavar="--range_of_size", type=str, help="Range of size of piRNAs")
    parser.add_argument("-t", metavar="--threads", type=int, help="Number of threads to use for clustering")
    # Argument used for quick testing
    #parser.add_argument("-s", metavar="--testing", type=int, help="USED FOR TESTING - Define number of scaffolds to "
    #                                                               "perform clustering on them")
    parser.add_argument("-m", metavar="--module", type=str, help="Which module to run: 'primary' for clustering "
                                                                 "algorithm, 'secondary' for ping-pong piRNA detection,"
                                                                 " 'both' for both")
    parser.add_argument("-v", metavar="--variation_threshold", type=int, help="Variation threshold for cleansing the reads")
    parser.add_argument("-a", action='store_true', help="Flag to automate parameter choosing")
    parser.add_argument("-d", action='store_true', help="Generate plots for each found cluster")
    parser.add_argument("--plot_iter", type=int, help="Number of plots generated per iteration. More == faster generation"
                                                      " == more RAM use")
    #parser.add_argument("--z_score_threshold", type=float, help="Threshold for z-score cutoff for piRNA range of size "
    #                                                            "detection coming from sample-wise analysis")
    #parser.add_argument("--cluster_quality_t_freq_threshold", type=float, help="Threshold for cluster quality analysis"
    #                                                                           "of frequency of Thymine on the first position of reads"
    #                                                                           "In form of ratio  0.5 == 50%%")
    #parser.add_argument("--cluster_quality_dist_freq_threshold", type=float, help="Threshold for cluster quality analysis"
    #                                                                           "of frequency of reads lengths of piRNA size within cluster"
    #                                                                           "In form of ratio  0.5 == 50%%")
    parser.add_argument("--version", action="store_true")
    args, leftovers = parser.parse_known_args()
    range_of_size, success = None, None
    if args.version:
        print(f"piRAT {__version__}")
        quit()
    print("----------------piRAT----------------")
    print(f"piRAT {__version__}")
    if args.m is None:
        args.m = 'both'
    if args.p is None:
        print('\nYou need to specify directory with RNAseq data to analyze, e.g: /data/')
        raw_path_name = input()
        if raw_path_name == '':
            sys.exit("Specified directory path is wrong")
    else:
        if args.p[0] != '/':
            new_path = os.getcwd() + '/' + args.p
            path_exist = os.path.exists(new_path)
            if path_exist:
                if args.p == '.':
                    raw_path_name = os.getcwd() + '/'
                else:
                    raw_path_name = new_path
            else:
                sys.exit("Specified directory path is wrong")
        else:
            if os.path.exists(args.p):
                raw_path_name = args.p
            else:
                sys.exit("Specified directory path is wrong")
    if args.o is None:
        current_directory = os.getcwd()
    else:
        current_directory = args.o
    list_of_files = data_prep.get_list_of_bam_files(raw_path_name)
    if len(list_of_files) == 0:
        sys.exit("There are no BAM files in the specified path!")
    if raw_path_name.endswith('/'):
        pass
    else:
        raw_path_name += '/'
    path = raw_path_name[:-1]
    name_of_the_directory = path.split('/')[-1]
    try:
        if name_of_the_directory not in os.listdir(current_directory):
            os.mkdir(f'{current_directory}/{name_of_the_directory}')
        else:
            pass
    except FileNotFoundError:
        sys.exit("Output path is wrong")
    print(f"Input path: {raw_path_name}")
    print(f"Output path: {current_directory}/{name_of_the_directory}/")
    current_directory += f'/{name_of_the_directory}/'
    saving_directory = current_directory
    autom_range = ''
    if args.a is False:
        if args.r is None:
            autom_range = ''
            while autom_range not in ['y', 'n']:
                print("\nYou didn't specify the range of size of piRNAs, do you want to automatically find the range of size"
                      " of piRNAs? [y/n]")
                autom_range = input()
                if autom_range == 'n':
                    print('\nPlease input your range as lower_limit,upper_limit. e.g: 26,32')
                    ranges = input()
                    range_of_size = [int(x) for x in ranges.split(',')]
                elif autom_range == 'y':
                    range_of_size = None
        else:
            range_of_size = [int(x) for x in args.r.split(',')]
    else:
        range_of_size = None if args.r is None else [int(x) for x in args.r.split(',')]
        autom_range = 'y'

    variation_threshold = 3 if args.v is None else args.v
    if args.plot_iter is None:
        plot_iter = 16
    else:
        plot_iter = args.plot_iter
    threads = 1 if args.t is None else int(args.t)
    draw_plots = False if args.d is False else args.d
    start = time.time()
    paired = False
    for file in list_of_files:
        is_paired = data_prep.check_if_paired(f"{raw_path_name}/{file}")
        if is_paired:
            print(f"{file} is from paired-end sequencing. Only first reads are gonna be processed in the analysis.\n")
            paired = True
    if args.a is False:
        is_automatic = True if autom_range == 'y' else False
    else:
        is_automatic = True
    if args.m == 'primary':
        conf.create_dirs(current_directory, 1, is_automatic, draw_plots)
        k, eps, threads, variation_threhsold, range_of_size = clustering(raw_path_name, range_of_size, threads, draw_plots, None, args.k, args.e,
                             variation_threshold, 1, [], plot_iter, True, saving_directory, is_automatic)
        rg.save_final_report_clustering(f"{saving_directory}", f"{saving_directory}", name_of_the_directory, current_time, raw_path_name, __version__, list_of_files, (k, eps, threads, variation_threhsold, range_of_size))
    elif args.m == 'secondary':
        conf.create_dirs(current_directory, 2, is_automatic, draw_plots)
        for file in list_of_files:
            print(f"Performing analysis of file: {file}...")
            name_of_the_file = '_'.join((file.split('.'))[:-1])
            conf.create_dir_ping_pong(f"{current_directory}/", name_of_the_file)
            sample_wise_analysis(raw_path_name, variation_threshold, [file], 0, current_directory, threads, paired, current_time, __version__)
            print("Done!")
        print("Performing final annotation...")
        ping_pong_annotate(raw_path_name, range_of_size, variation_threshold, list_of_files, 1, list_of_files, saving_directory, name_of_the_directory, threads, is_automatic, paired, current_time, __version__)
    elif args.m == 'both':
        conf.create_dirs(current_directory, 3, args.a, draw_plots)
        for file in list_of_files:
            print(f"Performing analysis of file: {file}...")
            name_of_the_file = '_'.join((file.split('.'))[:-1])
            conf.create_dir_ping_pong(f"{current_directory}/", name_of_the_file)
            success = sample_wise_analysis(raw_path_name, variation_threshold, [file], 0, current_directory, threads, paired, current_time, __version__)
            print("Done!")
        k, eps, threads, variation_threhsold, range_of_size = clustering(raw_path_name, range_of_size, threads, draw_plots, None, args.k, args.e,
                             variation_threshold, 2, list_of_files, plot_iter, success, current_directory, is_automatic)
        rg.save_final_report_clustering(f"{saving_directory}", f"{saving_directory}", name_of_the_directory, current_time, raw_path_name, __version__, list_of_files, (k, eps, threads, variation_threhsold, range_of_size))
        print("Performing final annotation...")
        ping_pong_annotate(raw_path_name, range_of_size, variation_threshold, list_of_files, 1, list_of_files, saving_directory, name_of_the_directory, threads, is_automatic, paired, current_time, __version__)
        try:
            rg.save_final_report(saving_directory, f"{saving_directory}",
                                      f"{saving_directory}", name_of_the_directory, range_of_size, is_automatic, current_time, raw_path_name, __version__, list_of_files, (k, eps, threads, variation_threhsold, range_of_size))
        except ValueError:
            rg.generete_venn(saving_directory + 'gffs/potential_pirna_reads.gff',
                                  saving_directory + 'gffs/ping_pong_reads.gff', saving_directory + 'gffs/clusters_out.gff',
                                  saving_directory, False)
            saving_path_plots = f"{saving_directory}plots/"
            rg.draw_plots_annotate(f"{saving_directory}",
                                        [["Initial_distribution", 'Initial_distribution.txt',
                                          "Length of the sequence [BP]", "No. sequences"],
                                         ["Number of pairs", 'Number_of_overlapping_reads.txt', "Nt overlap",
                                          "No. sequences"]], ['_'.join((x.split('.'))[:-1]) for x in list_of_files], saving_path_plots)
            pin_pon.generate_heatmap(saving_directory, {}, [range_of_size[0], range_of_size[1]], False, 1)
            print("There aren't any ping-pong signatures in all provided files!")
            rg.save_final_report_no_ping_pong(saving_directory, f"{saving_directory}",
                                      f"{saving_directory}/", name_of_the_directory, range_of_size, is_automatic, current_time, raw_path_name, __version__)
        except FileNotFoundError as e:
            rg.save_report_annotate_wo_primary_and_secondary(saving_directory, name_of_the_directory, current_time, raw_path_name, __version__)

    print(f"piRNA clustering took {round(time.time() - start, 2)} seconds!")
    pp.display_random_quote()
    return None

import matplotlib.pyplot as plt
from kneed import KneeLocator
from scipy.stats import linregress
import pirat.clustering_functions.data_preparation as data_prep
import pirat.configuration_functions.file_operations as f_op
from typing import Tuple, List, Dict, Any
import pysam
import pandas as pd
import numpy as np

plt.rcParams["figure.figsize"] = (20, 10)


def find_scope_of_the_size(list_of_data_files: List[str], raw_path_name: str, saving_path_primary: str) -> Tuple[int, int]:

    """
    Detect piRNA size in dataset based on read length distribution

    This function is loading 100,000 random (mapped) reads from each file in the dataset, to check read length distribution.
    Then it look for peak density read length in range 26-31, and calculates weighted mean where each consicutive abudant
    read has weight lower by 0.1 - most abundant = weight 1.0, 2nd most abundant = weight 0.9, ...
    After calculating weighted mean, lengths neighboring the peak density length are checked whether they meet the abpve mean criteria
    All connected read lengths are defining piRNA size range.

    Parameters:
        list_of_data_files (List[str]) - list of names of files in working directory  - e.g. [file_1.bam, file_2.bam, ..]
        raw_path_name (str) - path to the file containing dataset files
        saving_path_primary (str) - path to the output directory where resulting plot should be saved

    Returns:
        Range of size of piRNAs in format of Tuple - e.g (28,29)
    """
    len_of_seq = []
    # Loading 100 000 random reads of each data file
    for file in list_of_data_files:
        number_of_reads = 100000 / int(pysam.view(raw_path_name + file, '-c'))
        if number_of_reads > 1:
            number_of_reads = 1
        # Load 100 000 reads with the usage of -s flag (by percentage of all the reads calculated in
        # number_of_reads), flag -b ensures, that we only get mapped reads
        raw_data = (pysam.view(raw_path_name + file, f'-s {float(number_of_reads)}', '-F 0x4')).split("\n")[:-2]
        formatted_data = [read.split("\t") for read in raw_data]
        for read in formatted_data:
            len_of_seq.append(len(read[9]))

    # Getting a list of number of lengths for finding out maximums and minimums
    list_of_y = [len_of_seq.count(x) for x in range(50)]

    heights_in_range_26_32 = []
    for k in range(25, 31):
        heights_in_range_26_32.append(list_of_y[k])

    height_of_peak = max(heights_in_range_26_32)
    peak = list_of_y.index(height_of_peak)
    upper_limit = 0
    lower_limit = 0
    list_for_weight = (sorted(heights_in_range_26_32))[::-1]
    list_of_weights = []
    weight = 1.0
    for i in range(len(list_for_weight)):
        list_of_weights.append(weight)
        weight -= 0.1
    avg = [list_for_weight[x] * list_of_weights[x] for x in range(len(list_for_weight))]
    avg = sum(avg) / sum(list_of_weights)
    ll = 'keep'
    ul = 'keep'
    for i in range(1, 5):
        if list_of_y[peak - i] > avg and ll == 'keep':
            lower_limit = i

        elif list_of_y[peak - i] < avg:
            ll = 'stop'

        if list_of_y[peak + i] > avg and ul == 'keep':
            upper_limit = i

        elif list_of_y[peak + i] < avg:
            ul = 'stop'

    lower_limit = peak - lower_limit
    upper_limit = peak + upper_limit

    list_of_x = [number for number in range(50)]

    labels_plot = ["Length of the read (nt)", "Number of reads"]

    f_op.draw_scope_plot("Size", labels_plot, list_of_x, list_of_y, lower_limit, upper_limit, raw_path_name, saving_path_primary)

    labels = ["lower_limit", "Peak", "upper_limit"]
    data_lst = [lower_limit, peak, upper_limit]
    file_name = f"{saving_path_primary}/other_data/list_of_scopes.txt"
    with open(file_name, 'w') as file:
        for label, data in zip(labels, data_lst):
            file.write(label + '\t' + str(data) + '\n')

    return lower_limit, upper_limit


def get_ranges_for_different_k(ready_list_of_reads: List[List]) -> Dict[int, List]:
    """
    Preparation of distances between k-reads for automatic clustering parameter detection.
    For each k (4, 6, 8, 10) - loop through reads and calculate distance between them.

    Parameters:
        ready_list_of_reads - it is a list of chosen reads within the predefined range (by default (0, 100))
                          first part of the list consists of sorted reads from negative strand,
                          the second one contains sorted reads from positive strand [[negative_reads], positive_reads]]

    Returns:
        Dictionary, where each key is k value, and values for keys are lists of calculated inter-read distances.
    """
    list_of_ranges = []
    dic_of_all_ranges = {}
    for k in range(4, 11, 2):
        for strand in range(2):
            for read in range(len(ready_list_of_reads[strand]) - k):
                try:
                    list_of_ranges.append(
                        ready_list_of_reads[strand][read + k][1] - ready_list_of_reads[strand][read][0])
                except:
                    print(f"{len(ready_list_of_reads[strand])}\n{read + k}\n"
                          f"{ready_list_of_reads[strand][read + k]}\n{ready_list_of_reads[strand][read]}\n")
                    break
        dic_of_all_ranges[k] = list_of_ranges
        list_of_ranges = []
    return dic_of_all_ranges


def find_parameters(one_scaffold: str, list_of_data_files: List[str], raw_path_name: str, variation_threshold: int, range_of_size: Tuple[int, int])\
        -> Dict[int, List]:

    """
    Loading reads from all files of dataset into lists to calculate inter-read distances.

    Parameters:
        one_scaffold (str) - scaffold to load reads from (function is used in multithread mode, where individual
        scaffold processing are distributed onto threads)
        list_of_data_files (List[str]) - list of names of files in working directory  - e.g. [file_1.bam, file_2.bam, ..]
        raw_path_name (str) - path to the file containing dataset files
        variation_threshold (int) - Threshold value used to filter or determine acceptable variability across reads with the same core.
        range_of_size: A list specifying the minimum and maximum sizes of piRNAs ([min_size, max_size]).

    Returns:
        Dictionary, where each key is k value, and values for keys are lists of calculated inter-read distances.
    """
    list_of_reads_for_parameters = [[], []]
    for file in list_of_data_files:
        temp_reads = (data_prep.list_of_reads_for_finding_parameters(raw_path_name + file,
                                                                     one_scaffold, variation_threshold, range_of_size))
        list_of_reads_for_parameters[0] += temp_reads[0]
        list_of_reads_for_parameters[1] += temp_reads[1]
    ranges_for_different_k = get_ranges_for_different_k(list_of_reads_for_parameters)
    return ranges_for_different_k


def get_best_pair_of_parameters(dic_of_pairs: Dict[int, List], parameter_quality: List[bool]) -> Tuple[Any, Any]:
    """
    Automatic detection of best pair of parameters.
    Checking in reverse order (from highest to lowest k), if the parameter pair is fulfilling the quality criteria.
    The quality criteria pass is passed from outher function if for of parameter_quality list.

    Parameters:
        dic_of_pairs (Dict[int, List]) - Dict of pairs of parameters, where keys and k-values, and eps are detected Eps
        values
        parameter_quality (List[bool]) - Lits containing quality check passes for corresponding parameter pairs

    Returns:
        List containing best pair of clustering parameters - [k, Eps] -e.g [8, 1000]
    """
    values = list(dic_of_pairs.values())
    param_to_choose = None
    for set_of_params in range(3, 0, -1):
        if parameter_quality[set_of_params] == True:
            param_to_choose = set_of_params
            break
        else:
            continue
    if param_to_choose == None:
        return False, False
    else:
        return list(dic_of_pairs.keys())[param_to_choose], values[param_to_choose]


def finding_the_parameters(dic_of_ranges: Dict[int, List], raw_path_name: str) -> (Dict[int, List], List[bool]):
    """
    Finding parameter pairs (Eps for each k).
    Lists of inter-read distances for each k are loaded, for which index of distance closest to 10,000 nt is detected
    (piRNA clusters shoudln't have inter-read distances that large, and it limits the resources needed for caluclation)
    After descedning sorting of distances,  the resulting inter-read distances (for distances < 10,000 nt) are split
    into 20 equal parts, for which slope is calculated. "Knee" - the sudden change in inter-read distances distribution
    is detected based on created list of slopes. From detected knee position, the eps value is tracked from original list
    of inter-read distances.
    After detecting eps value, the similar process is happening for evaluation of detcted eps - inter-read distances
    lower then eps are splitted into 20 equal parts, for which slope is calculated. From there, if the change in the
    slope is lower than 10% for at least half of parts, the eps values passes the quality check.

    Parameters:
        dic_of_ranges: Dictionary, where each key is k value, and values for keys are lists of calculated inter-read distances.
        raw_path_name (str) - path to the output directory to save resulting plots

    Returns:
        Dict of pairs of parameters, where keys and k-values, and eps are detected Eps values
        List of bool values stating the quality pass check for each k, Eps pair
    """
    dic_of_eps = {}
    which_k = 0
    color = ['red', 'green', 'blue', 'brown']
    list_of_k = [4, 6, 8, 10]
    parameter_good = []
    for ranges_for_each_k in dic_of_ranges.values():
        sorted_ranges = sorted(ranges_for_each_k)
        sorted_ranges = [x for x in sorted_ranges if x >= 0]
        list_to_find_scopes_indexes = [xs for xs in range(len(sorted_ranges))]
        list_to_find_scopes = sorted_ranges[::-1]
        for find_the_index in range(10000, 5000, -1):
            try:
                found_index_cutoff = list_to_find_scopes.index(find_the_index)
                break
            except:
                continue
        list_to_find_scopes = list_to_find_scopes[found_index_cutoff:]
        z = int(len(list_to_find_scopes) / 20)
        list_to_find_scopes_indexes = list_to_find_scopes_indexes[found_index_cutoff:]
        list_of_slopes = []
        for iterator in range(0, len(list_to_find_scopes_indexes), z):
            list_of_slopes.append(linregress(list_to_find_scopes_indexes[iterator:iterator + z], list_to_find_scopes[iterator:iterator + z])[0])
        kn = KneeLocator([x for x in range(len(list_of_slopes) - 2)], list_of_slopes[:-2], curve="concave",
                         direction='increasing', online=False, S=0)

        found_index_of_eps = list_to_find_scopes_indexes.index((int(kn.knee) + 1) * z + list_to_find_scopes_indexes[0])
        y_eps = list_to_find_scopes[found_index_of_eps:]
        found_eps_value = list_to_find_scopes[found_index_of_eps]
        dic_of_eps[list_of_k[which_k]] = found_eps_value
        fig_name = f"{raw_path_name}/plots/Distances_for_k_{list_of_k[which_k]}.svg"
        f_op.draw_parameters_plot(fig_name=fig_name, labels = [["No. scope", "Value of slope"], ["Read number", "Eps value"], 
        [f"Distances for k = {list_of_k[which_k]}", "Found cutoff Eps"]], x = [list_to_find_scopes_indexes, [0, len(list_to_find_scopes_indexes) + list_to_find_scopes_indexes[0]]], y=[list_of_slopes, list_to_find_scopes, [found_eps_value, found_eps_value]], low_lim=-250, up_lim=10000)
        list_to_find_scopes_indexes = list_to_find_scopes_indexes[found_index_of_eps:]
        z = int(len(y_eps) / 20)
        list_of_slopes2 = []
        for iterator in range(0, len(list_to_find_scopes_indexes), z):
            list_of_slopes2.append(
                linregress(list_to_find_scopes_indexes[iterator:iterator + z], y_eps[iterator:iterator + z])[0])
        is_the_eps_cutoff_good = sum(x > -0.1 for x in list_of_slopes2) >= len(list_of_slopes2) / 2
        parameter_good.append(is_the_eps_cutoff_good)
        plt.clf()
        which_k += 1
    return dic_of_eps, parameter_good

def find_scope_with_sample_analysis(list_of_data_files: List[str], raw_path_name: str) -> List[int]:
    """
    Finding range of size of piRNAs based on sample-wise analysis results.
    Load z-score distribution of read lengths of sequences displaying 10 nts 5'-5' overlaps from each datafile analyzed
    in sample wise analysis. Average the values, and pick range of size of piRNAs for values displaying values higher
    than 0.5 z_score.

    Parameters:
        list_of_data_files (List[str]) - list of names of files in working directory  - e.g. [file_1.bam, file_2.bam, ..]
        raw_path_name (str) - path to the file containing dataset files

    Returns:
        Range of size of piRNAs in format of Tuple - e.g (28,29)
    """
    dfs_list = []
    for file in list_of_data_files:
        file_path = f'{raw_path_name}sample_wise_analysis/{file}/Length_distribution_of_the_sequences_displaying_10nts_5-to-5_overlaps.txt'
        try:
            dfs_list.append(pd.read_csv(file_path, sep='\t'))
        except FileNotFoundError:
            continue
    all_data = pd.concat(dfs_list)
    average_z_scores = all_data.groupby('x')['z-score'].mean().reset_index()
    result_df = average_z_scores[average_z_scores['z-score'] > 0.5]
    result_list = result_df['x'].to_list()
    return [min(result_list), max(result_list)]
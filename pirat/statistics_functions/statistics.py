import pysam
from Bio.Seq import Seq
import numpy as np
import matplotlib.pyplot as plt
import PyPDF2
import seqlogo
import multiprocessing
import os
import pandas as pd
from tqdm import tqdm
import logomaker
from typing import Tuple, List, Dict, Any
import csv
from collections import Counter
import re
from matplotlib.patches import FancyArrowPatch


def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    """
        Generate a key for natural sorting of strings containing numbers.

        Splitting the input string into a list of text and numeric parts - converting numeric parts into integers

        Parameters:
            s (str): string to be split to generated sorting list

        Returns:
            tuple: A tuple of strings and integers usable as a sort key.

        """
    return tuple(int(text) if text.isdigit() else text.lower() for text in re.split(_nsre, s))


def check_continuity(temp_list, strand_of_next_cluster, start_of_next_cluster):
    """
        Check if two clusters are continuous

        If clusters are on different strands, the function checks whether the new cluster
        overlaps with any cluster from the previous list on the opposite strand. It uses
        the maximum end position (stop) of relevant clusters to determine continuity.

        Parameters:

            temp_list (List[Tuple[str, Tuple[int, int], str]]): A list of previous clusters,
            where each item is a tuple containing (ID, (start, end), strand).
            strand_of_next_cluster (str): Strand of the next cluster ('+' or '-').
            start_of_next_cluster (int): Start position of the next cluster.

        Returns:
            bool: True if the next cluster is continuous with the previous ones, otherwise False.

        """
    if len(temp_list) >= 3:
        stop_minus = max([x[1][1] for x in temp_list[1:] if x[2] == '-'])
        stop_plus = max([x[1][1] for x in temp_list[1:] if x[2] == '+'])
        stop_to_take = stop_plus if strand_of_next_cluster == '-' else stop_minus
        is_continous = True if start_of_next_cluster <= stop_to_take else False
        return is_continous
    else:
        return False

def return_uni_strand(cluster_info):
    """
        Format a cluster entry as a uni-strand feature.

        Parameters:
            cluster_info (List or Tuple): A list containing genomic cluster information

        Returns:
            List[str]: A formatted cluster entry with added 'ID' and 'type=uni-strand' attributes.

        """
    cluster = [cluster_info[0], cluster_info[1], cluster_info[2], cluster_info[3], cluster_info[4],
               cluster_info[5], cluster_info[6], cluster_info[7], f'ID={cluster_info[8]};', 'type=uni-strand']
    return cluster


def return_joined_clusters(temp_list, length_of_current_cluster, length_of_next_cluster, cluster_info, variation_threshold):
    """
        Create a new cluster entry by joining two nearby clusters on opposite strands.

        Parameters:
            temp_list (List): A list of cluster entries, where each entry is a tuple or list
            containing (scaffold, (start, end), strand, ...).
            length_of_current_cluster (int): Length of the current cluster.
            length_of_next_cluster (int): Length of the next cluster.
            cluster_info (List): A list containing metadata for formatting the joined cluster
                (e.g., feature type, score, ID).
            variation_threshold (int): A small tolerance value to allow near-continuity in base pairs.

        Returns:
            List[str]: A formatted GFF-like row representing the joined cluster.

        """
    which_smaller = length_of_current_cluster if length_of_current_cluster < length_of_next_cluster else length_of_next_cluster
    tag_1 = f"tag1=({temp_list[1][2]}){temp_list[1][1][0]}-{temp_list[1][1][1]}"
    tag_2 = f"tag2=({temp_list[2][2]}){temp_list[2][1][0]}-{temp_list[2][1][1]}"
    if -10 <= (temp_list[2][1][0] - temp_list[1][1][1]) / which_smaller <= -0.5:
        new_row = [temp_list[1][0], cluster_info[1], cluster_info[2], temp_list[1][1][0], max([temp_list[1][1][1], temp_list[2][1][1]]),
                   cluster_info[5], '.', cluster_info[7], f'ID={cluster_info[8]};',
                   f'type=dual-strand;{tag_1};{tag_2}']

    elif (temp_list[2][1][0] - temp_list[1][1][1]) <= variation_threshold:
        new_row = [temp_list[1][0], cluster_info[1], cluster_info[2], temp_list[1][1][0], max([temp_list[1][1][1], temp_list[2][1][1]]),
                   cluster_info[5], '.', cluster_info[7], f'ID={cluster_info[8]};',
                   f'type=bi-directional;{tag_1};{tag_2}']
    return new_row

def return_big_cluster(list, temp_list):
    """
        Merge multiple subclusters into a single dual-strand cluster entry.

        This function combines all subclusters in `temp_list` (excluding index 0) into one large
        cluster spanning the full start-to-stop range. It adds GFF-style tags for each subcluster
        and returns a formatted feature line.

        Parameters:
            list (List): A list containing metadata used to format the final output
            (e.g., feature type, score, ID, strand, etc.).
            temp_list (List): A list of subclusters, each represented as a tuple or list
            with scaffold, (start, end), and strand.

        Returns:
            List[str]: A formatted GFF-like entry representing the merged dual-strand cluster.

        """
    list_of_tags = []
    for subcluster in range(1, len(temp_list)):
        temp_tag = f"tag{subcluster}=({temp_list[subcluster][2]}){temp_list[subcluster][1][0]}-{temp_list[subcluster][1][1]}"
        list_of_tags.append(temp_tag)

    list_of_start_pos = [x[1][0] for x in temp_list[1:]]
    list_of_stop_pos = [x[1][1] for x in temp_list[1:]]
    start_of_new_cluster = min(list_of_start_pos)
    stop_of_new_cluster = max(list_of_stop_pos)
    new_row = [temp_list[1][0], list[1], list[2], start_of_new_cluster, stop_of_new_cluster,
               list[5], '.', list[7], f'ID={list[8]};',
               f'type=dual-strand;{";".join(list_of_tags)}']
    return new_row

def classify_clusters(raw_path_name: str, eps: int) -> None:
    """
        Classify and merge nearby piRNA clusters based on strand orientation and genomic proximity.

        This function reads a GFF-formatted cluster file, identifies clusters that are close to each
        other and either overlapping or nearly continuous, and merges them into larger clusters of
        type 'dual-strand' or 'bi-directional'. Clusters that don't meet merging criteria are labeled
        as 'uni-strand'. The updated cluster annotations are written back to the same file.

        Merging rules:
            - Clusters on opposite strands that overlap or are within `variation_threshold` distance may be merged.
            - Merged clusters receive GFF-style `tag` annotations to trace original subclusters.
            - If three subclusters are merged, they may be classified as 'dual-strand'.
            - If more than three subclusters are merged, they are considered as a 'big cluster'.
        Parameters:
            raw_path_name (str): Path to the working directory containing `gffs/clusters_out.gff`.
            eps (int): Maximum allowed gap (in base pairs) between clusters for merging.

        Returns:
            None: The merged and classified clusters are saved back into the original GFF file.

        """

    df = pd.read_csv(f'{raw_path_name}/gffs/clusters_out.gff', sep='\t', header=None)
    df = df.sort_values(by=[0, 3, 4])
    i = 0
    temp_list = []
    new_list_of_clusters = []
    found = False
    is_continous = False
    while i < len(df):
        if len(df) != 1:
            if found:
                if i == len(df)-1 or i == len(df):
                    if temp_list == []:
                        new_cluster = return_uni_strand(df.iloc[i])
                        found = False
                    elif len(temp_list) == 3:
                        new_cluster = return_joined_clusters(temp_list, length_of_current_cluster, length_of_next_cluster,
                                                             df.iloc[i], eps)
                        found = True
                    else:
                        new_cluster = return_big_cluster(df.iloc[i], temp_list)
                        temp_list = []
                        found = True

                    new_list_of_clusters.append(new_cluster)
                    found = False
                    is_continous = False
                    break
                else:
                    if len(temp_list) == 3:
                        new_cluster = return_joined_clusters(temp_list, length_of_current_cluster, length_of_next_cluster,
                                                             df.iloc[i], eps)
                        found = True
                    elif len(temp_list) > 3:
                        new_cluster = return_big_cluster(df.iloc[i], temp_list)
                    new_list_of_clusters.append(new_cluster)
                    temp_list = []
                    found = False
                    is_continous = False
                    i+=1
                    if i == len(df)-1:
                        new_cluster = return_uni_strand(df.iloc[i])
                        new_list_of_clusters.append(new_cluster)
                        break
                    if new_cluster[-1] == 'type=uni-strand':
                        new_list_of_clusters.append(new_cluster)



            scaffold_of_current_cluster, scaffold_of_next_cluster = df.iloc[i, 0], df.iloc[i + 1, 0]
            strand_of_current_cluster, strand_of_next_cluster = df.iloc[i, 6], df.iloc[i + 1, 6]
            start_of_current_cluster, stop_of_current_cluster = df.iloc[i, 3], df.iloc[i, 4]
            start_of_next_cluster, stop_of_next_cluster = df.iloc[i + 1, 3], df.iloc[i + 1, 4]
            while (strand_of_current_cluster != strand_of_next_cluster or len(temp_list) != 0) and \
                    df.iloc[i][0] == df.iloc[i+1][0] and len(df)-1 > i:
                start_of_current_cluster, stop_of_current_cluster = df.iloc[i, 3], df.iloc[i, 4]
                start_of_next_cluster, stop_of_next_cluster = df.iloc[i + 1, 3], df.iloc[i + 1, 4]
                scaffold_of_current_cluster, scaffold_of_next_cluster = df.iloc[i, 0], df.iloc[i + 1, 0]
                strand_of_current_cluster, strand_of_next_cluster = df.iloc[i, 6], df.iloc[i + 1, 6]

                is_continous = check_continuity(temp_list, strand_of_next_cluster, start_of_next_cluster)
                if start_of_next_cluster < stop_of_current_cluster or \
                    start_of_next_cluster - stop_of_current_cluster <= eps or \
                    (len(temp_list) >= 3 and is_continous):
                    stop_of_new_cluster_plus = stop_of_current_cluster if strand_of_current_cluster == '+' else stop_of_next_cluster
                    stop_of_new_cluster_reverse = stop_of_next_cluster if strand_of_next_cluster == '-' else stop_of_current_cluster
                    start_of_new_cluster_plus = start_of_current_cluster if strand_of_current_cluster == '+' else start_of_next_cluster
                    start_of_new_cluster_reverse = start_of_next_cluster if strand_of_next_cluster == '-' else start_of_current_cluster
                    length_of_current_cluster, length_of_next_cluster = stop_of_current_cluster - start_of_current_cluster, stop_of_next_cluster - start_of_next_cluster

                    if len(temp_list) == 0:
                        temp_list = [[start_of_new_cluster_plus, start_of_new_cluster_reverse, stop_of_new_cluster_plus,
                                      stop_of_new_cluster_reverse],
                                     [scaffold_of_current_cluster, [start_of_current_cluster, stop_of_current_cluster],
                                      strand_of_current_cluster],
                                     [scaffold_of_next_cluster, [start_of_next_cluster, stop_of_next_cluster],
                                      strand_of_next_cluster]]
                        i += 1
                        found = True

                    else:
                        if stop_of_new_cluster_plus > temp_list[0][2]:
                            temp_list[0][2] = stop_of_new_cluster_plus
                        elif stop_of_new_cluster_reverse > temp_list[0][3]:
                            temp_list[0][3] = stop_of_new_cluster_reverse
                        elif start_of_new_cluster_reverse > temp_list[0][0]:
                            temp_list[0][0] = start_of_new_cluster_plus
                        elif start_of_new_cluster_reverse > temp_list[0][1]:
                            temp_list[0][1] = start_of_new_cluster_reverse

                        temp_list.append([scaffold_of_next_cluster, [start_of_next_cluster, stop_of_next_cluster],
                                          strand_of_next_cluster])

                        i += 1
                        found = True
                    if i == len(df) - 1:
                        break
                elif len(temp_list) == 3 and is_continous is False:
                    new_cluster = return_joined_clusters(temp_list, length_of_current_cluster, length_of_next_cluster, df.iloc[i], eps)
                    found = True
                    break
                elif len(temp_list) > 3:
                    new_cluster = return_big_cluster(df.iloc[i], temp_list)
                    temp_list = []
                    found = True
                    break
                else:
                    new_cluster = return_uni_strand(df.iloc[i])
                    temp_list = []
                    found = True
                    break
        elif len(df) == 1:
            new_cluster = return_uni_strand(df.iloc[i])
            new_list_of_clusters.append(new_cluster)
            found = True
            i+=1
        if found is False and is_continous is False:
            new_cluster = return_uni_strand(df.iloc[i])
            found = True


    temp_df = pd.DataFrame(new_list_of_clusters)
    temp_df = temp_df.drop_duplicates(subset=[3, 4, 9])
    temp_df['sort_key'] = temp_df.iloc[:, 0].apply(natural_sort_key)
    temp_df = temp_df.sort_values(by=['sort_key', 3]).drop('sort_key', axis=1)
    temp_df[10] = [f'Cluster_{i}' for i in range(1, len(temp_df) + 1)]
    temp_df[10] = temp_df[10] + ';' + temp_df[9]
    all_list = temp_df[[0, 1, 2, 3, 4, 5, 6, 7, 10]].values.tolist()
    df_to_save = pd.DataFrame(all_list)
    df_to_save.to_csv(f'{raw_path_name}/gffs/clusters_out.gff', sep='\t', index=False, header=False)
    return new_list_of_clusters

def transform_gff_file(path: str):
    """
        Aggregate and rewrite a GFF file by counting duplicate piRNA read entries.

        Parameters:
            path (str): Path to the directory containing the `gffs/Clusters_reads.gff` file.

        Returns:
            None: The function overwrites the original file with aggregated entries.

        """
    counter = Counter()
    with open(f'{path}/gffs/Clusters_reads.gff', 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            # Count the occurrences of each unique combination of seqname, start, and end
            counter[(row[0], row[3], row[4], row[6], row[8])] += 1

    with open(f'{path}/gffs/Clusters_reads.gff', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        for (seqname, start, end, strand, sequence), count in counter.items():
            writer.writerow([seqname, 'piRNA_annotation', 'Read', start, end, '.', strand, '0', sequence, count+1])
    return None


def get_no_reads_for_given_range(file_list: List, scaffold: str, range_start: int, range_stop: int, path: str
                                 , size_of_pirna=None) -> Tuple[List, List, List]:
    """
        Count and collect piRNA reads mapped to a specific genomic region across multiple BAM files.

        This function processes a list of BAM files, extracts reads that align to a given scaffold
        and genomic range, and optionally filters them by length. It returns:
            - The list of piRNA sizes (15–39 nt),
            - The number of reads of each size found in the specified range,
            - The sequences of all matching reads (reverse-complemented if on the reverse strand).

        Parameters:
            file_list (List[str]): List of BAM file names (relative to `path`).
            scaffold (str): Target scaffold/chromosome name to filter reads.
            range_start (int): Start position of the region of interest.
            range_stop (int): End position of the region of interest.
            path (str): Path to the directory containing the BAM files.
            size_of_pirna (Tuple[int, int], optional): Minimum and maximum allowed read lengths.
                If None, no size filtering is applied.

        Returns:
            Tuple[List[int], List[int], List[str]]:
                - A list of piRNA sizes (15–39 nt),
                - A list of counts corresponding to each size,
                - A list of filtered read sequences (reverse-complemented if needed).

        """
    sizes = [x for x in range(15, 40)]
    list_of_reads = []
    list_of_seqs = []
    for file in file_list:
        samfile = pysam.AlignmentFile(path + file)
        for read in samfile.fetch(scaffold):
            seq = read.query_sequence
            start_pos = read.pos
            if size_of_pirna is None:
                if range_start <= start_pos < range_stop:
                    list_of_reads.append(len(seq))
                    reverse = read.is_reverse
                    if reverse:
                        list_of_seqs.append((Seq(seq)).reverse_complement())
                    else:
                        list_of_seqs.append(seq)
            elif size_of_pirna[0] <= len(seq) <= size_of_pirna[1]:
                if range_start <= start_pos < range_stop:
                    list_of_reads.append(len(seq))
                    reverse = read.is_reverse
                    if reverse:
                        list_of_seqs.append((Seq(seq)).reverse_complement())
                    else:
                        list_of_seqs.append(seq)

    return sizes, [list_of_reads.count(size) for size in sizes], list_of_seqs


def get_matrix_of_seqs(x: List, upper_param: int):
    """
        Obtain a list of scaffold names from the given BAM file.

        Parameters:
            path_to_file (str): The path to a BAM file.

        Returns:
            List[str]: A list of scaffold names extracted from the BAM file.

        """
    letters = ['A', 'C', 'G', 'T']
    max_length_of_read = max([len(z) for z in x])
    max_length_of_read = max_length_of_read if max_length_of_read <= upper_param else upper_param
    list_of_nucleotides_on_given_position = []
    for i in range(max_length_of_read):
        list_of_nucleotides_on_given_position.append(list(map(lambda y: y[i] if len(y) > i else None, x)))
    list_of_nucleotides_on_given_position = [[x.count(m) for m in letters] for x in
                                             list(list_of_nucleotides_on_given_position)]
    return np.array(list_of_nucleotides_on_given_position)


def read_output_gff(file: str) -> Tuple[List, List, List]:
    """
        Read a GFF-like file and extract cluster coordinates, scaffold names, and length distribution.

        This function parses a tab-delimited GFF file and returns:
            - A list of clusters with scaffold, start, end, strand, and attributes.
            - A list of scaffold names for all clusters.
            - A list of cluster lengths (end - start).


        Parameters:
           file (str): Path to the GFF file containing cluster annotations.

        Returns:
            Tuple[List, List, List]:
                - List of clusters as [scaffold, start, end, strand, attributes].
                - List of scaffold names (one per cluster).
                - List of cluster length values.

        """
    list_of_scaffolds = []
    list_of_found_clusters = []
    cluster_length_distribution = []
    with open(file) as f:
        raw_data = [cluster.split('\t') for cluster in f.read().split('\n')][:-1]
        [list_of_found_clusters.append([x[0], int(x[3]), int(x[4]), x[6], x[8]]) for x in raw_data]
        [cluster_length_distribution.append(x[2] - x[1]) for x in list_of_found_clusters]
        [list_of_scaffolds.append(x[0]) for x in list_of_found_clusters]
    return list_of_found_clusters, list_of_scaffolds, cluster_length_distribution



def get_cluster_reads_out_of_gff(path:str) -> Tuple[Dict, Dict]:
    """
        Parse a GFF file to extract piRNA read information and compute positional read densities.

        This function reads `Clusters_reads.gff` and builds two dictionaries:
            - One that stores piRNA read intervals by scaffold and strand.
            - One that tracks positional density (sum of read counts) across genome coordinates.

        Each read is stored with:
            [start, end, read_length, count]

        Parameters:
            path (str): Path to the directory containing `gffs/Clusters_reads.gff`.

        Returns:
            Tuple[Dict, Dict]:
                - dict_of_reads_of_found_clusters:
                    Dict[scaffold] = [reverse_reads_list, forward_reads_list]
                - dict_of_density_of_the_reads_of_found_clusters:
                    Dict[scaffold] = [reverse_density_dict, forward_density_dict]

        """
    dict_of_reads_of_found_clusters = {}
    dict_of_density_of_the_reads_of_found_clusters = {}
    gff_of_reads_of_found_clusters = (open(f"{path}/gffs/Clusters_reads.gff", 'r')).read()
    gff_of_reads_of_found_clusters = [x.split('\t') for x in gff_of_reads_of_found_clusters.split('\n')]
    for i in gff_of_reads_of_found_clusters[:-1]:
        scaffold = i[0]
        if scaffold in dict_of_reads_of_found_clusters.keys():
            if i[-4] == '+':
                dict_of_reads_of_found_clusters[scaffold][1].append([int(i[3]), int(i[4]), len(i[-2]) + 1, int(i[-1])])
            else:
                dict_of_reads_of_found_clusters[scaffold][0].append([int(i[3]), int(i[4]), len(i[-2]) + 1, int(i[-1])])
            for xrange in range(int(i[3]), int(i[3]) + len(i[-2]) + 1):
                if i[-4] == '+':
                    if xrange in dict_of_density_of_the_reads_of_found_clusters[scaffold][1].keys():
                        dict_of_density_of_the_reads_of_found_clusters[scaffold][1][xrange] += int(i[-1])
                    else:
                        dict_of_density_of_the_reads_of_found_clusters[scaffold][1][xrange] = int(i[-1])
                else:
                    if xrange in dict_of_density_of_the_reads_of_found_clusters[scaffold][0].keys():
                        dict_of_density_of_the_reads_of_found_clusters[scaffold][0][xrange] += int(i[-1])
                    else:
                        dict_of_density_of_the_reads_of_found_clusters[scaffold][0][xrange] = int(i[-1])
        else:
            dict_of_reads_of_found_clusters[scaffold] = [[], []]
            dict_of_density_of_the_reads_of_found_clusters[scaffold] = [{}, {}]
            if i[-4] == '+':
                dict_of_reads_of_found_clusters[scaffold][1].append([int(i[3]), int(i[4]), len(i[-2]) + 1, int(i[-1])])
            else:
                dict_of_reads_of_found_clusters[scaffold][0].append([int(i[3]), int(i[4]), len(i[-2]) + 1, int(i[-1])])
            for xrange in range(int(i[3]), int(i[3]) + len(i[-2]) + 1):
                if i[-4] == '+':
                    if xrange in dict_of_density_of_the_reads_of_found_clusters[scaffold][1].keys():
                        dict_of_density_of_the_reads_of_found_clusters[scaffold][1][xrange] += int(i[-1])
                    else:
                        dict_of_density_of_the_reads_of_found_clusters[scaffold][1][xrange] = int(i[-1])
                else:
                    if xrange in dict_of_density_of_the_reads_of_found_clusters[scaffold][0].keys():
                        dict_of_density_of_the_reads_of_found_clusters[scaffold][0][xrange] += int(i[-1])
                    else:
                        dict_of_density_of_the_reads_of_found_clusters[scaffold][0][xrange] = int(i[-1])
    return dict_of_reads_of_found_clusters, dict_of_density_of_the_reads_of_found_clusters


def draw_arrow(ax, x_location, sign, maximum, color, cluster_length):
    """
        Draw a directional arrow representing a piRNA cluster strand orientation.

        Depending on the strand sign, this function adds an arrow to a matplotlib Axes:
            - Arrows for '+' strand are drawn near the top, pointing right.
            - Arrows for '-' strand are drawn near the bottom, pointing left.

        Parameters:
            ax (matplotlib.axes.Axes): The axes to which the arrow will be added.
            x_location (int): The x-coordinate where the arrow starts.
            sign (str): Strand direction, either '+' or '-'.
            maximum (int): Maximum value used to determine vertical positioning.
            color (str): Color of the arrow.
            cluster_length (int): Total length of the cluster, used to scale the arrow.

        Returns:
            None

        """
    # Depending on the sign, set the direction of the arrow
    y_location = int(maximum * 0.95)
    x_step = int(cluster_length * 0.05)
    if sign == '+':
        # Draws an arrow on the upper side heading right
        arrow = FancyArrowPatch((x_location, y_location), (x_location + x_step, y_location),
                                arrowstyle='-|>', mutation_scale=30, color=color)

    else:
        # Draws an arrow on the bottom side heading left
        arrow = FancyArrowPatch((x_location, -y_location), (x_location - x_step, -y_location),
                                arrowstyle='-|>', mutation_scale=30, color=color)
    ax.add_patch(arrow)


def get_statistics_and_create_plots(analyzed_cluster: List, list_of_files: List,
                                    upper_param: int, draw_plots: bool, list_of_reads_for_both_strands: List,
                                    occurence_dict_for_both_strands: Dict, scope: Tuple[int, int], raw_path_name: str, saving_path_primary:str) \
                                    -> List:
    """
    Compute statistics and generate plots for a given piRNA cluster.

    This function:
        - Computes per-base coverage, read length distribution, and sequence logo statistics
        - Optionally draws and saves multiple plots:
            - Coverage and density bar plot with directional arrows
            - Read length distribution
            - Sequence logos for all and unique reads
        - Merges all generated plots into a single PDF
        - Classifies the cluster based on read composition (e.g., 1U bias)

    Parameters:
        analyzed_cluster (List): Cluster information (GFF-like row including scaffold, start, end, strand, and attributes).
        list_of_files (List[str]): List of BAM file names used to fetch reads.
        upper_param (int): Length of the longest read to consider for sequence logos.
        draw_plots (bool): Whether to draw and save plots.
        list_of_reads_for_both_strands (List[List]): List of reads for reverse and forward strands.
        occurence_dict_for_both_strands (Dict): Per-position read count dictionaries for each strand.
        scope (Tuple[int, int]): Read length window (e.g., (24, 32)) for filtering and density statistics.
        raw_path_name (str): Path to input data (e.g., BAM or metadata).
        saving_path_primary (str): Path where plots and PDFs should be saved.

    Returns:
        List: A pair of two elements:
            - If the cluster passes 1U + scope-based criteria: [cluster_info, '']
            - Otherwise: ['', cluster_info]
    """
    plt.rcdefaults()
    plt.rcParams.update({'font.size': 20})
    cluster_type = (analyzed_cluster[-1].split(';'))[1:]
    scaffold = analyzed_cluster[0]
    list_of_reads = [list_of_reads_for_both_strands[0],
                     list_of_reads_for_both_strands[1]]
    occurrence_dict = [dict(sorted(occurence_dict_for_both_strands[0].items())),
                       dict(sorted(occurence_dict_for_both_strands[1].items()))]
    list_of_reads_positions, occurence_of_scope_read, distribution_of_scope_reads_in_cluster, \
        length_of_reads_bar, no_of_reads_of_given_length, frequency_of_reads, \
        frequency_of_unique_reads, length_of_reads, average, subtype, clusters_info = prepare_data_for_plots(list_of_reads,
                                                                                     occurrence_dict,
                                                                                     analyzed_cluster,
                                                                                     list_of_files, scaffold,
                                                                                     analyzed_cluster[3],
                                                                                     upper_param, raw_path_name, cluster_type)

    if draw_plots:

        fig, (ax1, ax2) = plt.subplots(1, 2, clear=True, num=1)

        fig.set_size_inches(20, 10)

        ax1.spines['bottom'].set_position('center')

        ax1.spines['right'].set_color('none')

        ax1.spines['top'].set_color('none')

        check = int((analyzed_cluster[2] - analyzed_cluster[1]) / 30)

        ax1.set_xticks([analyzed_cluster[1] - check, analyzed_cluster[2] + check],
                       [analyzed_cluster[1], analyzed_cluster[2]])

        ax1.bar(list_of_reads_positions, occurence_of_scope_read, color='gray',
                label='Presence of putative piRNA read')
        ax1.bar(list_of_reads_positions, distribution_of_scope_reads_in_cluster,
                color='brown', label='Density of putative piRNA reads')
        maximum = max(distribution_of_scope_reads_in_cluster) if analyzed_cluster[3] == '+' \
            else -min(distribution_of_scope_reads_in_cluster)
        lim_value = max(distribution_of_scope_reads_in_cluster) if max(distribution_of_scope_reads_in_cluster) > np.abs(
            min(distribution_of_scope_reads_in_cluster)) else np.abs(min(distribution_of_scope_reads_in_cluster))

        if subtype == 1 or subtype == 2:
            counter = 0
            for subcluster in clusters_info:
                if subcluster[0] == '-':
                    pos, strand = subcluster[1][1], '-'
                else:
                    pos, strand = subcluster[1][0], '+'
                if counter == 0:
                    ax1.axvline(x=pos, ymin=-lim_value, ymax=lim_value, linestyle='--', color='green',
                            label=f'Start of pre-merged sub-cluster')
                    counter += 1
                else:
                    ax1.axvline(x=pos, ymin=-lim_value, ymax=lim_value, linestyle='--', color='green')
                draw_arrow(ax1, pos, strand, lim_value, 'green', analyzed_cluster[2] - analyzed_cluster[1])

            ax1.set_yticks([lim_value, -lim_value], [f'{lim_value}   -', f'{lim_value}   +'])

            ax1.set_ylim([-lim_value, lim_value])
        else:
            ax1.set_yticks([-max(distribution_of_scope_reads_in_cluster),
                            max(distribution_of_scope_reads_in_cluster)],
                           [f'{max(distribution_of_scope_reads_in_cluster)}   -',
                            f'{max(distribution_of_scope_reads_in_cluster)}   +']) \
                if analyzed_cluster[3] == '+' else ax1.set_yticks(
                [min(distribution_of_scope_reads_in_cluster),
                 -min(distribution_of_scope_reads_in_cluster)],
                [f'{-min(distribution_of_scope_reads_in_cluster)}   -',
                 f'{-min(distribution_of_scope_reads_in_cluster)}   +'])

            ax1.set_ylim([-max(distribution_of_scope_reads_in_cluster),
                          max(distribution_of_scope_reads_in_cluster)]) \
                if analyzed_cluster[3] == '+' else ax1.set_ylim(
                [min(distribution_of_scope_reads_in_cluster),
                 -min(distribution_of_scope_reads_in_cluster)])

        ax1.annotate(f'Per base read coverage of {(analyzed_cluster[4]).split(";")[0]}',
                     xy=(0.5, 1.10), xycoords='axes fraction', fontsize=30,
                     ha='center', transform=plt.gca().transAxes)

        ax1.annotate('Len. Cluster: %.3f kBP; Num. Reads [%s - %s nt]: %d; Read Density: %.3f' %
                     (round((analyzed_cluster[2] - analyzed_cluster[1]) / 1000, 3), str(scope[0]),
                      str(scope[1]), sum(no_of_reads_of_given_length[scope[0] - 15: scope[1] - 14]),
                      round(average, 3)), xy=(0.5, 1.05),
                     xycoords='axes fraction', fontsize=15, ha='center',
                     transform=plt.gca().transAxes)

        ax1.xaxis.set_tick_params(size=0)

        ax1.spines['bottom'].set_position('center')

        ax1.tick_params(axis='x', labelrotation=60)

        ax1.legend(loc='upper left') if analyzed_cluster[3] == '-' else ax1.legend(loc='lower left')

        ax1.set_ylabel('Density of putative piRNA reads')

        ax2.bar(length_of_reads_bar, no_of_reads_of_given_length,
                width=0.5, align='center')

        ax2.annotate(f'{(analyzed_cluster[4]).split(";")[0]} read len. distribution',
                     xy=(0.5, 1.10), xycoords='axes fraction', fontsize=30,
                     ha='center', transform=plt.gca().transAxes)

        ax2.annotate(f'Num. reads: {sum(no_of_reads_of_given_length)}',
                     xy=(0.5, 1.05), xycoords='axes fraction', fontsize=15,
                     ha='center', transform=plt.gca().transAxes)

        ax2.set_xlabel('Length of the read [BP]')

        ax2.set_ylabel('Number of the reads')

        plt.savefig(f'{saving_path_primary}/plots/Clusters/{(analyzed_cluster[4]).split(";")[0]}_distribution.pdf')

        fig, (ax3, ax4) = plt.subplots(2, 1, clear=True, num=1)

        fig.set_size_inches(20, 10)

        seqlogo_plot = logomaker.Logo(frequency_of_reads.ppm, ax=ax3)

        seqlogo_plot.ax.set_xlim(-0.5, upper_param - 0.5)

        seqlogo_plot.ax.set_xticks(range(0, upper_param + 1, 5))

        seqlogo_plot.ax.set_xticklabels([1, 5, 10, 15, 20, 25, 30, 35][:(int(upper_param / 5) + 1)])

        seqlogo_plot.ax.set_title("Seqlogo of all reads of given cluster range")

        seqlogo_plot.ax.set_ylabel('Frequency')

        seqlogo_plot = logomaker.Logo(frequency_of_unique_reads.ppm, ax=ax4)

        seqlogo_plot.ax.set_xlim(-0.5, upper_param - 0.5)

        seqlogo_plot.ax.set_xticks(range(0, upper_param + 1, 5))

        seqlogo_plot.ax.set_xticklabels([1, 5, 10, 15, 20, 25, 30, 35][:(int(upper_param / 5) + 1)])

        seqlogo_plot.ax.set_title("Seqlogo of unique reads of given cluster range")

        seqlogo_plot.ax.set_ylabel("Frequency")
        plt.savefig(f'{saving_path_primary}/plots/Clusters/{(analyzed_cluster[4]).split(";")[0]}_seqlogo.pdf')

        merger = PyPDF2.PdfMerger()

        for pdf in [f'{saving_path_primary}/plots/Clusters/{(analyzed_cluster[4]).split(";")[0]}_distribution.pdf',
                    f'{saving_path_primary}/plots/Clusters/{(analyzed_cluster[4]).split(";")[0]}_seqlogo.pdf']:
            merger.append(pdf)
            os.remove(pdf)
        merger.write(f'{saving_path_primary}/plots/Clusters/{(analyzed_cluster[4]).split(";")[0]}.pdf')
        merger.close()

    analyzed_cluster += [sum(no_of_reads_of_given_length), len(length_of_reads)]

    freq_of_t = frequency_of_reads.pm['T'][0]

    all_of_occ = sum(no_of_reads_of_given_length)

    scope_occ = sum(no_of_reads_of_given_length[scope[0] - 15: scope[1] - 14])
    if freq_of_t > 0.5 and scope_occ / all_of_occ > 0.5:
        quality = [analyzed_cluster, '']
    else:
        quality = ['', analyzed_cluster]
    return quality


def prepare_data_for_plots(list_of_reads: List, occurrence_dict: Dict, analyzed_cluster: List, list_of_files: List,
                           scaffold: str, strand: str, upper_param: int, raw_path_name: str, cluster_type: str) -> Tuple[List, List, \
                           List, List, List, List, Any, Any, float]:
    """
    Process read and cluster data to prepare coverage and sequence logo plot inputs.

    This function extracts various statistics for a given piRNA cluster:
        - Read coverage and density across the cluster region
        - Read length distribution
        - Sequence logo data for all and unique reads
        - Subtype classification (uni-strand, dual-strand, bi-directional)

    Parameters:
        list_of_reads (List): Reads from both strands, formatted as:
            [ [rev_reads], [fwd_reads] ], where each read is [start, end, length, count].
        occurrence_dict (Dict): Two dictionaries (per strand) with positional read counts.
        analyzed_cluster (List): Cluster information [scaffold, start, end, strand, attributes].
        list_of_files (List[str]): List of BAM file names to extract read sequences.
        scaffold (str): Scaffold/chromosome name for the cluster.
        strand (str): Strand of the cluster ('+', '-', or '.').
        upper_param (int): Max read length to consider for sequence logo construction.
        raw_path_name (str): Path to directory containing GFF and BAM files.
        cluster_type (str): Type of cluster ('uni-strand', 'bi-directional', 'dual-strand', with tags).

    Returns:
        Tuple containing:
            - list_of_reads_positions (List[int]): Genomic positions of reads across cluster.
            - occurence_of_scope_read (List[int]): Placeholder markers for presence of reads (strand-aware).
            - distribution_of_scope_reads_in_cluster (List[int]): Per-base read density values.
            - length_of_reads_bar (List[int]): Read lengths used in histogram.
            - no_of_reads_of_given_length (List[int]): Frequency of each read length.
            - frequency_of_reads (seqlogo.CompletePm): Position matrix for all reads (for logo).
            - frequency_of_unique_reads (seqlogo.CompletePm): Matrix for unique reads (for logo).
            - length_of_reads (List[int]): Individual read lengths used.
            - average (float): Average read coverage across the cluster.
            - subtype (int): Cluster subtype (0: uni-strand, 1: bi-directional, 2: dual-strand).
            - clusters_info (List): List of subcluster metadata (strand, [start, end]) if subtype > 0.
    """
    length_of_reads = []
    if 'uni-strand' in cluster_type[0]:
        subtype = 0
        clusters_info = []
    else:
        pattern = r'tag\d=\(([-+])\)(\d+)-(\d+)'
        clusters_info = []
        for item in cluster_type:
            match = re.match(pattern, item)
            if match:
                sub_strand = match[1]  # strand is the first capture group
                sub_start = int(match[2])  # start coordinate is the second capture group
                sub_end = int(match[3])  # end coordinate is the third capture group
                clusters_info.append([sub_strand, [sub_start, sub_end]])

        if 'bi' in cluster_type[0]:
            subtype = 1
        elif 'dual' in cluster_type[0]:
            subtype = 2
    if strand == '.':
        is_strand_forward = 2
        for read in list_of_reads[0] + list_of_reads[1]:
            length_of_reads.append(int(read[2]) - 1) if analyzed_cluster[1] <= read[1] <= analyzed_cluster[2] else None
        joined_dict = dict(Counter(occurrence_dict[0]) + Counter(occurrence_dict[1]))
    else:
        is_strand_forward = 1 if strand == '+' else 0
        other_strand = 1 if is_strand_forward == 0 else 0

        for read in list_of_reads[is_strand_forward]:
            length_of_reads.append(int(read[2]) - 1) if analyzed_cluster[1] <= read[1] <= analyzed_cluster[2] else None

        for read in list_of_reads[other_strand]:
            length_of_reads.append(int(read[2]) - 1) if analyzed_cluster[1] <= read[1] <= analyzed_cluster[2] else None

    if is_strand_forward == 2:
        list_of_reads_positions = sorted(position for position in list(occurrence_dict[0].keys()) \
                                                         if analyzed_cluster[1] <= position <= analyzed_cluster[2])
        list_of_reads_positions_other_strand = sorted(position for position in list(occurrence_dict[1].keys()) \
                                                         if analyzed_cluster[1] <= position <= analyzed_cluster[2])
        list_of_reads_positions_all = list_of_reads_positions + list_of_reads_positions_other_strand
    else:
        list_of_reads_positions = sorted(position for position in list(occurrence_dict[is_strand_forward].keys()) \
                                                         if analyzed_cluster[1] <= position <= analyzed_cluster[2])
        list_of_reads_positions_other_strand = sorted(position for position in list(occurrence_dict[other_strand].keys()) \
                                                         if analyzed_cluster[1] <= position <= analyzed_cluster[2])



    if is_strand_forward == 2:
        distribution_of_scope_reads_in_cluster =[occurrence_dict[0][position]
                                                               for position in list_of_reads_positions]


        distribution_of_scope_reads_in_cluster_other_strand = [occurrence_dict[1][position]
                                                               for position in list_of_reads_positions_other_strand]


        average = sum(distribution_of_scope_reads_in_cluster + distribution_of_scope_reads_in_cluster_other_strand) / (
                analyzed_cluster[2] - analyzed_cluster[1])

        distribution_of_scope_reads_in_cluster = [-y for y in \
                                                               distribution_of_scope_reads_in_cluster]

        maxim = max(list_of_reads_positions_all)
        #lim_value = max(distribution_of_scope_reads_in_cluster) if max(distribution_of_scope_reads_in_cluster) > np.abs(
        #    min(distribution_of_scope_reads_in_cluster)) else np.abs(min(distribution_of_scope_reads_in_cluster))
        occurence_of_scope_read = [-maxim if y != 0 else 0 for y  in
                                   distribution_of_scope_reads_in_cluster]
        occurence_of_scope_read += [maxim if y != 0 else 0 for y in distribution_of_scope_reads_in_cluster_other_strand]


    else:
        distribution_of_scope_reads_in_cluster = [occurrence_dict[is_strand_forward][position]
                                                               for position in list_of_reads_positions]

        distribution_of_scope_reads_in_cluster_other_strand = [occurrence_dict[other_strand][position]
                                                               for position in list_of_reads_positions_other_strand]

        distribution_of_scope_reads_in_cluster = [-y for y in distribution_of_scope_reads_in_cluster] \
                                                    if is_strand_forward == 0 else \
                                                    [y for y in distribution_of_scope_reads_in_cluster]
        distribution_of_scope_reads_in_cluster_other_strand = [-y for y in \
                                                               distribution_of_scope_reads_in_cluster_other_strand] \
            if other_strand == 0 else [y for y in distribution_of_scope_reads_in_cluster_other_strand]
        average = sum(distribution_of_scope_reads_in_cluster) / (analyzed_cluster[2] - analyzed_cluster[1]) \
            if strand == '+' else -sum(distribution_of_scope_reads_in_cluster) / (
                analyzed_cluster[2] - analyzed_cluster[1])


        occurence_of_scope_read = [max(distribution_of_scope_reads_in_cluster) for y in
                                   distribution_of_scope_reads_in_cluster if y != 0] \
            if is_strand_forward == 1 else [min(distribution_of_scope_reads_in_cluster) for y in
                                            distribution_of_scope_reads_in_cluster if y != 0]
        occurence_of_scope_read += [-max(distribution_of_scope_reads_in_cluster) for y in
                                   distribution_of_scope_reads_in_cluster_other_strand if y != 0] \
            if is_strand_forward == 1 else [-min(distribution_of_scope_reads_in_cluster) for y in
                                            distribution_of_scope_reads_in_cluster_other_strand if y != 0]

    list_of_reads_positions += list_of_reads_positions_other_strand
    distribution_of_scope_reads_in_cluster += distribution_of_scope_reads_in_cluster_other_strand

    length_of_reads_bar, no_of_reads_of_given_length, read_seqs = get_no_reads_for_given_range(list_of_files,
                                                                                               scaffold,
                                                                                               analyzed_cluster[1],
                                                                                               analyzed_cluster[2],
                                                                                               raw_path_name)

    matrix_of_basis = get_matrix_of_seqs(read_seqs, upper_param)
    matrix_of_unique_basis = get_matrix_of_seqs(list(set(read_seqs)), upper_param)
    frequency_of_reads = seqlogo.CompletePm(matrix_of_basis)
    frequency_of_unique_reads = seqlogo.CompletePm(matrix_of_unique_basis)
    return list_of_reads_positions, occurence_of_scope_read, distribution_of_scope_reads_in_cluster, \
        length_of_reads_bar, no_of_reads_of_given_length, frequency_of_reads, frequency_of_unique_reads, \
        length_of_reads, average, subtype, clusters_info


def statistics_of_found_clusters(raw_path_name: str, high_quality: List, low_quality: List, scope:Tuple, saving_path_primary: str) -> None:
    """
    Generate summary statistics and visualizations for high- and low-quality piRNA clusters.

    This function:
        - Computes basic statistics for cluster length, total reads, and in-scope reads (e.g., 24–32 nt).
        - Plots histograms of cluster lengths for all, high-, and low-quality clusters.
        - Saves statistics to both CSV and HTML formats for reporting and downstream use.

    Parameters:
        raw_path_name (str): Path to the directory with raw input files (not directly used here).
        high_quality (List): List of high-quality clusters (each entry: [scaffold, start, end, strand, ...]).
        low_quality (List): List of low-quality clusters (same format).
        scope (Tuple[int, int]): Range of read lengths considered as "in-scope" (e.g., (24, 32)).
        saving_path_primary (str): Path to save the output plots and statistics files.

    Returns:
        None: Output is written to disk as plots and summary tables.
    """
    
    plt.rcParams.update({'font.size': 26})
    length_distribution_high_quality_clusters = [x[2] - x[1] for x in high_quality]
    length_distribution_low_quality_clusters = [x[2] - x[1] for x in low_quality]
    if len(length_distribution_high_quality_clusters) != 0 and len(length_distribution_low_quality_clusters) != 0:
        max_length = max(max(length_distribution_high_quality_clusters), max(length_distribution_low_quality_clusters))
    elif len(length_distribution_high_quality_clusters) == 0 and len(length_distribution_low_quality_clusters) != 0:
        max_length = max(length_distribution_low_quality_clusters)
    else:
        max_length = max(length_distribution_high_quality_clusters)
    bins = [x for x in range(0, max_length, 10)]
    plt.close('all')
    fig, ax = plt.subplots(1,1,clear=True,num=1)
    plt.rcParams.update({'font.size': 26})
    ax.hist(length_distribution_high_quality_clusters + length_distribution_low_quality_clusters, bins=bins, lw=1.2)
    fig.set_size_inches(10, 10)
    ax.set_title(f'All clusters')
    ax.set_ylabel('Number of clusters')
    ax.set_xlabel('Cluster length [BP]')
    plt.savefig(f'{saving_path_primary}/plots/Cluster_length_distribution.svg')
    fig, (ax1, ax2) = plt.subplots(1, 2, clear=True)
    fig.set_size_inches(20, 10)
    ax1.hist(length_distribution_high_quality_clusters, bins=bins, lw=1.2)
    ax2.hist(length_distribution_low_quality_clusters, bins=bins, lw=1.2)
    ax1.set_ylabel('Number of clusters')
    ax2.set_ylabel('Number of clusters')
    ax2.set_xlabel('Cluster length [BP]')
    ax1.set_xlabel('Cluster length [BP]')
    ax2.set_title('Low quality clusters')
    ax1.set_title('High quality clusters')
    plt.savefig(f'{saving_path_primary}/plots/Cluster_length_distribution_quality.svg')
    plt.clf()
    number_of_clusters = [len(length_distribution_high_quality_clusters), len(length_distribution_low_quality_clusters)]
    if number_of_clusters[0] == 0:
        avg_length_of_clusters = [np.mean(length_distribution_low_quality_clusters)]
        number_of_reads_low_quality_all, number_of_scope_reads_low_quality = get_number_of_reads_for_found_clusters(
            low_quality)
        statistics_of_number_of_reads = [(sum(x), np.mean(x), np.std(x), sum(x)) for x in
                                         [
                                          number_of_reads_low_quality_all, number_of_scope_reads_low_quality]]
        df_of_statistics = {"High quality":
                                [f"0", f"0",
                                 f"0",
                                 f"0",
                                 f"0",
                                 f"0"],
                            "Low quality":
                                [f"{number_of_clusters[1]}", f"{round(avg_length_of_clusters[0], 3)}",
                                 f"{statistics_of_number_of_reads[0][0]}",
                                 f"{round(statistics_of_number_of_reads[0][1], 3)}",
                                 f"{statistics_of_number_of_reads[1][0]}",
                                 f"{round(statistics_of_number_of_reads[1][1], 3)}"],
                            "Total":
                                [f"{number_of_clusters[1]}", f"{round(avg_length_of_clusters[0], 3)}",
                                 f"{statistics_of_number_of_reads[0][0]}",
                                 f"{round(statistics_of_number_of_reads[0][1], 3)}",
                                 f"{statistics_of_number_of_reads[1][0]}",
                                 f"{round(statistics_of_number_of_reads[1][1], 3)}"]
                            }
    elif number_of_clusters[1] == 0:
        avg_length_of_clusters = [np.mean(length_distribution_high_quality_clusters)]
        number_of_reads_high_quality_all, number_of_scope_reads_high_quality = get_number_of_reads_for_found_clusters(
            high_quality)
        statistics_of_number_of_reads = [(sum(x), np.mean(x), np.std(x), sum(x)) for x in
                                         [
                                             number_of_reads_high_quality_all, number_of_scope_reads_high_quality]]
        df_of_statistics = {"High quality":
                                [f"{number_of_clusters[0]}", f"{round(avg_length_of_clusters[0], 3)}",
                                 f"{statistics_of_number_of_reads[0][0]}",
                                 f"{round(statistics_of_number_of_reads[0][1], 3)}",
                                 f"{statistics_of_number_of_reads[1][0]}",
                                 f"{round(statistics_of_number_of_reads[1][1], 3)}"],
                            "Low quality":
                                [f"0", f"0",
                                 f"0",
                                 f"0",
                                 f"0",
                                 f"0"],
                            "Total":
                                [f"{number_of_clusters[0]}", f"{round(avg_length_of_clusters[0], 3)}",
                                 f"{statistics_of_number_of_reads[0][0]}",
                                 f"{round(statistics_of_number_of_reads[0][1], 3)}",
                                 f"{statistics_of_number_of_reads[1][0]}",
                                 f"{round(statistics_of_number_of_reads[1][1], 3)}"]
                            }
    else:
        avg_length_of_clusters = [np.mean(length_distribution_high_quality_clusters),
                                  np.mean(length_distribution_low_quality_clusters)]
        number_of_reads_high_quality_all, number_of_scope_reads_high_quality = get_number_of_reads_for_found_clusters(
            high_quality)
        number_of_reads_low_quality_all, number_of_scope_reads_low_quality = get_number_of_reads_for_found_clusters(
            low_quality)
        statistics_of_number_of_reads = [(sum(x), np.mean(x), np.std(x), sum(x)) for x in
                                         [number_of_reads_high_quality_all, number_of_scope_reads_high_quality,
                                          number_of_reads_low_quality_all, number_of_scope_reads_low_quality]]
        df_of_statistics = {"High quality":
                                [f"{number_of_clusters[0]}", f"{round(avg_length_of_clusters[0], 3)}",
                                 f"{statistics_of_number_of_reads[0][0]}",
                                 f"{round(statistics_of_number_of_reads[0][1], 3)}",
                                 f"{statistics_of_number_of_reads[1][0]}",
                                 f"{round(statistics_of_number_of_reads[1][1], 3)}"],
                            "Low quality":
                                [f"{number_of_clusters[1]}", f"{round(avg_length_of_clusters[1], 3)}",
                                 f"{statistics_of_number_of_reads[2][0]}",
                                 f"{round(statistics_of_number_of_reads[2][1], 3)}",
                                 f"{statistics_of_number_of_reads[3][0]}",
                                 f"{round(statistics_of_number_of_reads[3][1], 3)}"],
                            "Total":
                                [f"{sum(number_of_clusters)}", f"{round(np.mean(length_distribution_high_quality_clusters + length_distribution_low_quality_clusters), 3)}",
                                 f"{statistics_of_number_of_reads[0][0] + statistics_of_number_of_reads[2][0]}",
                                 f"{round(np.mean(number_of_reads_low_quality_all + number_of_reads_high_quality_all), 3)}",
                                 f"{len(number_of_scope_reads_low_quality + number_of_scope_reads_high_quality)}",
                                 f"{round(np.mean(number_of_scope_reads_low_quality + number_of_scope_reads_high_quality), 3)}"]
                            }
    df_of_statistics = pd.DataFrame(df_of_statistics, index=["Number of clusters", "Average cluster length",
                                                             "Number of all reads in clusters",
                                                             "Average number of reads per cluster",
                                                             f"Number of {scope[0]}-{scope[1]}nts reads in clusters",
                                                             f"Average number {scope[0]}-{scope[1]}nts reads per cluster"
                                                             f" per cluster"])
    df_of_statistics.to_csv(f'{saving_path_primary}/other_data/Statistics_of_found_clusters.csv')
    df_of_statistics.to_html(f'{saving_path_primary}/other_data/Statistics_of_found_clusters.html')

    return None


def get_number_of_reads_for_found_clusters(list_of_clusters: List) -> List[List]:
    """
    Extract total and in-scope read counts from annotated piRNA clusters.

    Each cluster is expected to contain precomputed metadata in positions:
        - x[5]: Total number of reads in the cluster
        - x[6]: Number of reads within a specific size range (e.g., 24–32 nt)

    Parameters:
        list_of_clusters (List): A list of cluster entries, where each cluster is a list
                                 containing read count information at specific indices.

    Returns:
        List[List]: A two-element list:
            - List of total read counts for each cluster.
            - List of in-scope (size-filtered) read counts for each cluster.
    """
    reads_all = [x[5] for x in list_of_clusters]
    reads_of_found_size = [x[6] for x in list_of_clusters]
    return [reads_all, reads_of_found_size]


def run_the_statistics(list_of_clusters_for_statistics: List, upper_param: int, threads: int, list_of_files: List,
                       draw_plots: bool, list_of_clusters: List, list_of_clusters_of_reads: List, scope: Tuple[int, int],
                       raw_path_name: str, plot_iter: int, saving_path_primary:str) -> Tuple[List, List]:
    """
    Perform statistical analysis and quality classification of detected piRNA clusters.

    This function evaluates each cluster using `get_statistics_and_create_plots`, either
    with or without graphical output. Clusters are processed in batches and in parallel
    using multiple threads. Based on read composition (e.g., 1U bias and size distribution),
    clusters are classified into high or low quality.

    Parameters:
        list_of_clusters_for_statistics (List): List of clusters to analyze. Each entry is assumed
            to be a cluster tuple containing an identifier and other metadata.
        upper_param (int): Maximum read length to consider in sequence logo plots.
        threads (int): Number of parallel threads for processing.
        list_of_files (List[str]): List of BAM files used to extract read data.
        draw_plots (bool): Whether to generate and save plots during analysis.
        list_of_clusters (List): Mapping from cluster ID to full cluster data.
        list_of_clusters_of_reads (List): Mapping from cluster ID to corresponding read data.
        scope (Tuple[int, int]): Read length range to consider as "in-scope" (e.g., (24, 32)).
        raw_path_name (str): Path to the raw data directory.
        plot_iter (int): Number of clusters to process per plotting batch.
        saving_path_primary (str): Directory path to save output plots and results.

    Returns:
        Tuple[List, List]:
            - high_all (List): All clusters classified as high-quality.
            - low_all (List): All clusters classified as low-quality.
    """
    low_all = []
    high_all = []
    stat = 'statistical analysis with graphical output' if draw_plots is True else 'statistical analysis'
    print(f"Performing {stat} of found clusters...")
    number_of_plots_per_iteration = plot_iter
    for i in tqdm(range(0, len(list_of_clusters_for_statistics), number_of_plots_per_iteration)):
        if len(list_of_clusters_for_statistics[i:]) < number_of_plots_per_iteration:
            k = [(cluster, list_of_files, upper_param, draw_plots, list_of_clusters[cluster[0]],
                  list_of_clusters_of_reads[cluster[0]], scope, raw_path_name, saving_path_primary) for cluster in
                 list_of_clusters_for_statistics[i:]]
            pool = multiprocessing.Pool(threads)
            results = pool.starmap(get_statistics_and_create_plots, k)
            pool.close()
        else:
            k = [(cluster, list_of_files, upper_param, draw_plots, list_of_clusters[cluster[0]],
                  list_of_clusters_of_reads[cluster[0]], scope, raw_path_name, saving_path_primary) for cluster in
                 list_of_clusters_for_statistics[i:i + number_of_plots_per_iteration]]
            pool = multiprocessing.Pool(threads)
            results = pool.starmap(get_statistics_and_create_plots, k)
            pool.close()
        high_quality = [x[0] for x in results if x[0] != '']
        low_quality = [x[1] for x in results if x[1] != '']
        high_all += high_quality
        low_all += low_quality
    return high_all, low_all

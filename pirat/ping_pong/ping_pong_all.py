import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import seaborn as sns
from collections import Counter
from typing import Tuple, List, Dict
import seqlogo
from Bio.Seq import Seq
import multiprocessing
from intervaltree import IntervalTree
from itertools import groupby
from operator import itemgetter
from matplotlib_venn import venn2
from pandas.errors import EmptyDataError
import pyranges as pr
import matplotlib.patches as mpatches
import pirat.clustering_functions.data_preparation as data_prep

def load_data(path, list_of_bam_files, epsilon, range_of_size, flag, threads, paired):
    """
    Loads BAM files and processes piRNA sequence data.

    Parameters:
        path (str): Path to the directory containing BAM files.
        list_of_bam_files (list): List of BAM file names.
        epsilon (int): Allowed positional error margin for read overlaps.
        range_of_size (tuple): Tuple indicating min and max allowed read lengths.
        flag (int): Flag to filter reads based on BAM flags.
        threads (int): Number of threads for multiprocessing.
        paired (bool): Indicates if the BAM files contain paired-end data.

    Returns:
        tuple: All piRNA sequences, overlapping sequences, read lengths, and a dictionary of piRNA mappings.
    """
    list_of_bam_files = [f"{path}/{bam_file}" for bam_file in list_of_bam_files]
    list_of_scaffolds = pysam.AlignmentFile(list_of_bam_files[0], 'rb').references
    all_pirna_sequences, all_overlapping_sequences, length_of_all_reads, all_pirna_dict = manage_data(list_of_bam_files, list_of_scaffolds, epsilon, range_of_size, flag, threads, paired)
    return all_pirna_sequences, all_overlapping_sequences, length_of_all_reads, all_pirna_dict


def cleanse_reads(list_of_reads, dict_of_abundace, epsilon, range_of_size):
    """
    Cleans and filters sequencing reads, removing duplicates and resolving overlaps.

    Parameters:
        list_of_reads (list): List of read tuples with genomic coordinates and sequences.
        dict_of_abundace (dict): Dictionary containing read abundance counts.
        epsilon (int): Allowed positional tolerance when resolving overlaps.
        range_of_size (tuple): Allowed read length range.

    Returns:
        tuple: Filtered list of reads and list of overlaps removed.
    """
    overlaps = []
    final_list_pos = []
    list_of_reads = list(set(list_of_reads))
    sorted_data = sorted(list_of_reads, key=itemgetter(0))
    grouped_data = groupby(sorted_data, key=itemgetter(0))
    result = [list(group) for _, group in grouped_data]
    for group_of_reads in result:
        if len(group_of_reads) > 1:
            list_to_add, list_to_del, chosen_read = chose_the_final_list(group_of_reads, dict_of_abundace, 1, epsilon,
                                                                         range_of_size)
            final_list_pos += list_to_add
        else:
            final_list_pos += group_of_reads
    final_list = []
    sorted_data = sorted(final_list_pos, key=itemgetter(1))
    grouped_data = groupby(sorted_data, key=itemgetter(1))
    result = [list(group) for _, group in grouped_data]
    for group_of_reads in result:
        if len(group_of_reads) > 1:
            list_to_add, list_to_del, chosen_read = chose_the_final_list(group_of_reads, dict_of_abundace, 0, epsilon,
                                                                         range_of_size)
            final_list += list_to_add
        else:
            final_list += group_of_reads
    return final_list, overlaps


def chose_the_final_list(temp_dict, dict_of_abundace, which_strand, epsilon, range_of_size):
    """
    Selects representative reads from groups of overlapping reads based on abundance and length.

    Parameters:
        temp_dict (list): Grouped overlapping reads.
        dict_of_abundace (dict): Abundance information for reads.
        which_strand (int): Indicates strand orientation (0 or 1).
        epsilon (int): Allowed positional tolerance.
        range_of_size (tuple): Allowed read length range.

    Returns:
        tuple: Final selected reads, indicator of chosen reads, and the minimal representative read.
    """
    chosen_reads = 0
    not_good = []
    temp = {read: dict_of_abundace[read] for read in temp_dict}
    temp_list_of_reads = []
    temp = {key: value for key, value in sorted(temp.items())}
    if int(sum(list(temp.values()))) > len(list(temp.values())):
        values = list(temp.values())
        max_count = values.count(max(values))
        if max_count > 1:
            keys = []
            minimal = sorted(list(temp.keys()))[0]
            for key, value in temp.items():
                if value == max_count:
                    keys.append(key)
            for key in keys:
                if range_of_size[0] <= len(key[2]) <= range_of_size[1]:
                    temp_list_of_reads.append(key)
                    for read in temp_list_of_reads:
                        if len(read[2]) > len(minimal[2]):
                            minimal = read
                            break
        else:
            max_key = max(temp, key=temp.get)
            minimal = max_key
            chosen_reads = 1
    else:
        minimal = sorted(list(temp.keys()))[0]
        sorted_keys = sorted(list(temp.keys()))
        for mm in sorted_keys:
            if range_of_size[0] <= len(mm[2]) <= range_of_size[1]:
                temp_list_of_reads.append(mm)
        for read in temp_list_of_reads:
            if len(read[2]) > len(minimal[2]):
                minimal = read
                break


    for x in temp_dict:
        if x != minimal and -epsilon > int(x[which_strand] - minimal[which_strand]) <= epsilon:
            not_good.append(x)
    final_list = [minimal]
    for item in not_good:
        if item not in final_list:
            final_list.append(item)
    final_list = list(set(final_list))
    return final_list, chosen_reads, minimal

def check_if_paired(file):
    """
    Checks if a BAM file predominantly contains paired-end sequencing reads.

    Parameters:
        file (str): Path to the BAM file.

    Returns:
        bool: True if >50% of reads are paired, otherwise False.
    """
    bam = pysam.AlignmentFile(file, "rb")
    paired_reads = 0
    total_reads = 0
    for read in bam:
        total_reads += 1
        if read.is_paired:
            paired_reads += 1
        if total_reads > 10000:
            break
    return paired_reads / total_reads > 0.5


def manage_data(list_of_bam_files, list_of_scaffolds, epsilon, range_of_size, flag, threads, paired):
    """
    Manages data processing across scaffolds using parallel multiprocessing.

    Parameters:
        list_of_bam_files (list): List of BAM file paths.
        list_of_scaffolds (list): Scaffold identifiers from BAM files.
        epsilon (int): Positional tolerance for overlaps.
        range_of_size (tuple): Allowed read lengths.
        flag (int): BAM file filtering flag.
        threads (int): Number of threads for parallel processing.
        paired (bool): Indicates paired-end data.

    Returns:
        tuple: All piRNA sequences, overlapping sequences, read lengths, and piRNA dictionary.
    """
    k = [(scaffold, epsilon, range_of_size, list_of_bam_files, flag, paired) for scaffold in list_of_scaffolds]
    with multiprocessing.Pool(threads) as p:
        results = p.starmap(multiproces, k)
    all_pirna_sequences = [item for sublist in [x[0] for x in results] for item in sublist]
    all_overlapping_sequences = [item for sublist in [x[1] for x in results] for item in sublist]
    pirna_dict = [x[2] for x in results]
    updated_dict = {}
    [updated_dict.update(x) for x in pirna_dict]
    length_of_all_reads = [item for sublist in [x[3] for x in results] for item in sublist]
    return all_pirna_sequences, all_overlapping_sequences, length_of_all_reads, updated_dict

def multiproces(scaffold, epsilon, range_of_size, list_of_bam_files, flag, paired):
    """
    Processes reads from a specific scaffold in parallel.

    Parameters:
        scaffold (str): Scaffold identifier.
        epsilon (int): Positional overlap tolerance.
        range_of_size (tuple): Read length constraints.
        list_of_bam_files (list): BAM file paths.
        flag (int): Filtering flag.
        paired (bool): Indicates paired-end data.

    Returns:
        tuple: piRNA sequences, overlapping sequences, piRNA dictionary, and read lengths.
    """
    length_of_reads_list = []
    wyniki_, wyniki_2 = mp_m_d(scaffold, list_of_bam_files, range_of_size, flag, paired)
    wyniki = wyniki_
    length_of_reads_list += wyniki_2

    reads_rev, reads_pos, ab_rev, ab_pos = wyniki[0], wyniki[1], wyniki[2], wyniki[3]
    data_of_scaffold_rev, data_of_scaffold_pos = cleanse_reads(reads_rev, ab_rev, epsilon, range_of_size), \
        cleanse_reads(reads_pos, ab_pos, epsilon, range_of_size)
    data_of_scaffold = data_of_scaffold_rev[0] + data_of_scaffold_pos[0]
    all_pirna_sequences, all_overlapping_sequences, pirna_dict = find_overlapping_sequences(data_of_scaffold)
    return all_pirna_sequences, all_overlapping_sequences, pirna_dict, length_of_reads_list

def get_abundance(list_of_reads):
    """
    Calculates abundance for each unique sequencing read.

    Parameters:
        list_of_reads (list): Reads for abundance calculation.

    Returns:
        dict: Dictionary of read abundances.
    """
    dict_of_abundance = {}
    for read in list_of_reads:
        if read in dict_of_abundance.keys():
            dict_of_abundance[read] += 1
        else:
            dict_of_abundance[read] = 1
    return dict_of_abundance

def mp_m_d(scaffold, list_of_bam_files, range_of_size, flag, paired):
    """
    Prepares read data for a single scaffold across multiple BAM files.

    Parameters:
        scaffold (str): Scaffold identifier.
        list_of_bam_files (list): BAM file paths.
        range_of_size (tuple): Allowed read lengths.
        flag (int): Filtering flag.
        paired (bool): Indicates paired-end data.

    Returns:
        tuple: Scaffold reads organized by strand, abundances, and read lengths.
    """
    number = 0
    length_of_all_reads = []
    data_of_scaffold_rev, data_of_scaffold_pos = [], []
    for bam_file in list_of_bam_files:
        temp_data_of_scaffold, temp_length_of_all_reads = data_prep.get_reads(bam_file, scaffold, flag, range_of_size, paired)
        data_of_scaffold_rev += temp_data_of_scaffold[0]
        data_of_scaffold_pos += temp_data_of_scaffold[1]
        length_of_all_reads += temp_length_of_all_reads
    number += len(data_of_scaffold_rev) + len(data_of_scaffold_pos)
    abun_rev, abun_pos = get_abundance(data_of_scaffold_rev), get_abundance(data_of_scaffold_pos)
    return [data_of_scaffold_rev, data_of_scaffold_pos, abun_rev, abun_pos], length_of_all_reads

def find_overlapping_sequences(data_of_scaffold):
    """
    Identifies overlapping piRNA sequences from given scaffold data.

    Parameters:
        data_of_scaffold (list): Scaffold reads with genomic coordinates.

    Returns:
        tuple: piRNA sequences, overlapping potential sequences, and piRNA dictionary.
    """

    col = ['Start', 'End', 'Sequence', 'Strand', 'Chromosome']
    df = pd.DataFrame(data_of_scaffold, columns = col)
    df['Strand'] = df['Strand'].apply(lambda x: '+' if x is True else '-')
    if df.shape[0] < 2:
        return [], [], {}
    else:
            gr = pr.PyRanges(df)
            gr_plus = gr[gr.Strand == '+']
            gr_minus = gr[gr.Strand == '-']
            if gr_plus.length == 0 or gr_minus.length == 0:
                return [], [], {}
            else:
                overlaps = gr_plus.join(gr_minus, strandedness=False, slack=10, apply_strand_suffix=True)
                if overlaps.length == 0 :
                    return [], [], {}
                else:
                    df = overlaps.df
                    # Calculate the length of each overlap
                    df['OverlapLength_1'] = df['End_b'] - df['Start'] + 1
                    df['OverlapLength_2'] = df['End'] - df['Start_b'] + 1
                    df['Strand'] = df['Strand'].apply(lambda x: str_to_bool(x))
                    df['Strand_b'] = df['Strand_b'].apply(lambda x: str_to_bool(x))
                    # Filter overlaps that have an overlap length of exactly 10
                    exact_10_df = df[(df['OverlapLength_2'] == 10)]

                    # Extract the sequences of the reads that overlap by exactly 10
                    pirna_dict = pd.Series(exact_10_df[['Start_b', 'End_b', 'Sequence_b', 'Strand_b',
                                                        'Chromosome']].apply(tuple, axis=1).values, index=exact_10_df[
                        ['Start', 'End', 'Sequence', 'Strand', 'Chromosome']].apply(tuple, axis=1)).to_dict()
                    df['MinOverlapLength'] = df[['OverlapLength_2']].min(axis=1)
                    df = df[df['MinOverlapLength'] > 0]
                    overlapping_potential_list = df['MinOverlapLength'].tolist()
                    pirna_list = [item for sublist in list(pirna_dict.items()) for item in sublist]
                    return pirna_list, overlapping_potential_list, pirna_dict



def str_to_bool(s):
    """
    Converts strand indicator from string to boolean.

    Parameters:
        s (str): Strand representation as '+' or '-'.

    Returns:
        bool: True if '+' else False.
    """
    return s.lower() == "+"


def build_interval_tree(data):
    """
    Builds an interval tree from genomic coordinate data for efficient querying.

    Parameters:
        data (list): List of genomic intervals with associated data.

    Returns:
        IntervalTree: Interval tree of provided genomic intervals.
    """
    tree = IntervalTree()
    for item in data:
        start, end, _, val, _ = item
        tree[start:end] = item
    return tree


def find_overlapping_pairs(data, tree, overlap_threshold):
    """
    Finds overlapping read pairs based on an interval tree and overlap threshold.

    Parameters:
        data (list): List of reads with genomic coordinates.
        tree (IntervalTree): Interval tree built from genomic data.
        overlap_threshold (int): Minimum overlap required between reads.

    Returns:
        list: List of overlapping read pairs.
    """
    overlapping_pairs = []
    for item in data:
        start, end, _, val1, _ = item
        overlaps = tree.envelop(start + overlap_threshold - 1, end - overlap_threshold)
        for interval in overlaps:
            _, _, _, val2, _ = interval.data
            if val1 != val2:
                overlapping_pairs.append((item, interval.data))
    return overlapping_pairs

def get_matrix_of_seqs(x: List, upper_param: int, path: str, saving_path, path_to_gff = None):
    """
    Generates and saves sequence logos for provided piRNA sequences.

    Parameters:
        x (list): piRNA sequences.
        upper_param (int): Maximum length for sequences.
        path (str): Base path for reading data.
        saving_path (str): Path for saving sequence logo plots.
        path_to_gff (str, optional): Path to GFF file for alternative input.

    Returns:
        None
    """
    if path_to_gff is not None:
        gff_df = pd.read_csv(path_to_gff, sep='\t', header=None)
        x = gff_df[8].tolist()
        title = 'SeqLogo of primary piRNAs'
        saving_name = "seqlogo_clusters.svg"
    else:
        title = "SeqLogo of ping-pong piRNAs"
        saving_name = "seqlogo_pingpong.svg"
    letters = ['A', 'C', 'G', 'T']
    max_length_of_read = max([len(z) for z in x])
    max_length_of_read = max_length_of_read if max_length_of_read <= upper_param else upper_param
    list_of_nucleotides_on_given_position = []
    for i in range(max_length_of_read):
        list_of_nucleotides_on_given_position.append(list(map(lambda y: y[i] if len(y) > i else None, x)))
    list_of_nucleotides_on_given_position = [[x.count(m) for m in letters] for x in
                                             list(list_of_nucleotides_on_given_position)]
                                             
    plt.rcParams.update({'font.size': 20})
    fig, ax1 = plt.subplots(1, 1, clear=True, num=1)
    fig.set_size_inches(15,5)
    seqlogo_plot = logomaker.Logo(seqlogo.CompletePm(np.array(list_of_nucleotides_on_given_position)).ppm, ax=ax1)
    seqlogo_plot.ax.set_xlim(-0.5, upper_param - 0.5)
    seqlogo_plot.ax.set_xticks([0, 4, 9, 14, 19, 24, 29, 34][:(int(upper_param / 5) + 1)])
    seqlogo_plot.ax.set_xticklabels([1, 5, 10, 15, 20, 25, 30, 35][:(int(upper_param / 5) + 1)])
    seqlogo_plot.ax.set_title(title)
    seqlogo_plot.ax.set_ylabel('Frequency')
    plt.tight_layout()
    fig.savefig(f'{saving_path}/{saving_name}')
    del fig
    plt.clf()
    return None


def draw_plots(path, list_of_data, saving_path):
    """
    Generates and saves distribution plots for piRNA data.

    Parameters:
        path (str): Base data path.
        list_of_data (list): Data details for plotting.
        saving_path (str): Path for saving plot images.

    Returns:
        None
    """
    for data in list_of_data:
        if len(data) == 0:
            data_to_plot = {x:0 for x in range(20, 41)}
        else:
            data_to_plot = dict(Counter(data[1]))
        sorted_data_to_plot = sorted(data_to_plot.items())
        x, y = zip(*sorted_data_to_plot)
        y_z_score = [round((elem - np.mean(y)) / np.std(y), 2) for elem in y]
        create_plots(path, data[0], x, y, data[2], data[3], saving_path)
        create_plots(path, data[0], x, y_z_score, data[2], "z-score", saving_path)

        data_for_txt = pd.DataFrame({'x': x, 'y': y, 'z-score': y_z_score})
        save_raw_data(data_for_txt, data[0], f"{path}")
    return None


def draw_plots_annotate(path, list_of_data, files, saving_path):
    """
    Creates annotated plots comparing piRNA data across multiple samples.

    Parameters:
        path (str): Path to the base data.
        list_of_data (list): Details for plotting.
        files (list): List of sample file identifiers.
        saving_path (str): Path to save the annotated plots.

    Returns:
        None
    """
    cmap = plt.cm.get_cmap('viridis', len(files))
    for data in list_of_data:
        plt.rcParams.update({'font.size': 16})
        fig, ax1 = plt.subplots(1, 1, clear=True)
        fig2, ax2 = plt.subplots(1, 1, clear=True)
        fig2.set_size_inches(7, 5)
        fig.set_size_inches(7, 5)
        failed = 0
        if 'Length_distribution' in data[1]:
            name = "Length_distribution_of_the_sequences_displaying_10nts_5-to-5_overlaps"
        elif "5-to-5_overlap" in data[1]:
            name = "5-to-5_overlap"
        else:
            name = "Initial_distribution"
        for idx, i in enumerate(files):
            try:
                all_data = pd.read_csv(f"{path}sample_wise_analysis/{i}/{name}.txt", sep='\t')
                x, y, z_score = all_data['x'].to_numpy().tolist(), all_data['y'].to_numpy().tolist(), all_data[
                    'z-score'].to_numpy().tolist()
                color = cmap(idx)
                ax1.plot(x, y, label=i, marker='o', color=color)
                ax1.axis(xmin=min(x), xmax=max(x))
                max_y = max(y)
                x_to_y = x[y.index(max_y)]
                ax1.grid()
                ax1.text(x_to_y, max_y, str(x_to_y), fontsize=16)
                if data[0] == "Initial_distribution":
                    ax1.set_title("Initial distribution")
                else:
                    ax1.set_title(data[0])
                ax1.set_xlabel(data[2])
                ax1.set_ylabel(data[3])
                ax2.plot(x, z_score, label=i, marker='o', color=color)
                ax2.axis(xmin=min(x), xmax=max(x))
                max_y = max(z_score)
                x_to_y = x[z_score.index(max_y)]
                ax2.grid()
                ax2.text(x_to_y, max_y, str(x_to_y), fontsize=16)
                if data[0] == "Initial_distribution":
                    ax2.set_title("Initial distribution")
                else:
                    ax2.set_title(data[0])
                ax2.set_xlabel(data[2])
                ax2.set_ylabel('z-score')
            except:
                failed += 1
                continue
        if failed == len(files):
            continue
        else:
            ax1.legend(loc='upper left', bbox_to_anchor=(1, 1))
            fig.savefig(f'{saving_path}{name}.svg', bbox_inches='tight')
            ax2.legend(loc='upper left', bbox_to_anchor=(1, 1))
            fig2.savefig(f'{saving_path}{name}_z_score.svg', bbox_inches='tight')
    return None


def save_raw_data(data, name, path)->None:
    """
    Saves processed data to a text file.

    Parameters:
        data (DataFrame): Data to save.
        name (str): Name of the file.
        path (str): Directory path to save the data.

    Returns:
        None
    """
    name = name.replace(' ','_')
    name = name.replace("'", "")
    if 'Length_distribution_of_the_sequences' in name:
        name = "Length_distribution_of_the_sequences_displaying_10nts_5-to-5_overlaps"
    data.to_csv(f"{path}/{name}.txt", sep='\t', index=False)
    return None

def save_to_gff(path: str, pirna_list: List[List]) -> None:
    """
    Writes piRNA sequence data into a GFF file format.

    Parameters:
        path (str): Output path for the GFF file.
        pirna_list (list): List of piRNA sequences and coordinates.

    Returns:
        None
    """
    path = f"{path}gffs/ping_pong_reads.gff"
    counter, pair_id = 1, 1
    with open(path, "w") as file:
        for read in pirna_list:
            strand = '+' if read[3] == False else '-'
            row = f"{read[-1]}\t.\tread\t{read[0]}\t{read[1]}\t.\t{strand}\t.\tID=pingpong_pirna_{counter};PingPongPair=pingpong_pair_{pair_id};Seq={str(read[2])}\n"
            pair_id += 1 if counter % 2 == 0 else 0
            file.write(row)
            counter += 1

def prepare_data_for_heatmap(pirna_dict: Dict, range_: Tuple) -> Dict:
    """
    Prepares piRNA data for generating heatmaps showing pairwise sequence length relationships.

    Parameters:
        pirna_dict (dict): Dictionary mapping paired piRNAs.
        range_ (tuple): Range of piRNA lengths to consider.

    Returns:
        dict: Counts for each piRNA length pair.
    """
    dict_pairs_counter_for_heatmap = {}
    unique_lengths_lst = list(range(range_[0], range_[1]+1))
    for i in unique_lengths_lst:
        for j in unique_lengths_lst:
            my_key = str(i) + str(j)
            dict_pairs_counter_for_heatmap[my_key] = 0
    for key, values in pirna_dict.items():
            if len(key[2]) == len(values[2]):
                created_key = str(len(key[2])) + str(len(values[2]))
                if created_key in dict_pairs_counter_for_heatmap:
                    dict_pairs_counter_for_heatmap[created_key] += 1
                else:
                    dict_pairs_counter_for_heatmap[created_key] = 1
            else:
                created_key1 = str(len(key[2])) + str(len(values[2]))
                created_key2 = str(len(values[2])) + str(len(key[2]))
                if created_key1 in dict_pairs_counter_for_heatmap and created_key2 in dict_pairs_counter_for_heatmap:
                    dict_pairs_counter_for_heatmap[created_key1] += 1
                    dict_pairs_counter_for_heatmap[created_key2] += 1
                else:
                    dict_pairs_counter_for_heatmap[created_key1] = 1
                    dict_pairs_counter_for_heatmap[created_key2] = 1
    return dict_pairs_counter_for_heatmap


def generate_heatmap(path: str, dict_pairs: Dict, range_, status, type) -> None:
    """
    Creates and saves a heatmap representing piRNA pair length distributions.

    Parameters:
        path (str): Path to save the heatmap image.
        dict_pairs (dict): Pair counts for heatmap.
        range_ (tuple): Range of lengths to display.
        status (bool): Flag indicating presence of data.
        type (int): Heatmap style or mode indicator.

    Returns:
        None
    """
    plt.close('all')
    if type == 2:
        for i in range(range_[0], range_[1]+ 1):
            for j in range(range_[0], range_[1]+ 1):
                text = f'{i}{j}'
                if text in dict_pairs.keys():
                    continue
                else:
                    dict_pairs[text] = 0
    dict_pairs = dict(sorted(dict_pairs.items()))
    unique_lengths_lst = list(range(range_[0], range_[1] + 1))
    number = range_[1] - range_[0] + 1
    if status:
        matrix = np.zeros((number, number)).astype(int)
        dic_val = list(dict_pairs.values())
        z = 0
        for i in range(number):
            for j in range(number):
                matrix[j][i] = int(dic_val[z])
                z += 1
    else:
        matrix = np.array({k : 0 for k in range(range_[0], range_[1])})
    mask = np.ones_like(matrix, dtype=bool)
    try:
        for i in range(mask.shape[0]):
            for j in range(mask.shape[1]):
                if j < i + 1:
                    mask[i, j] = False

        matrix = np.where(mask, 0, matrix) if type == 1 else matrix
    except IndexError or ValueError:
        matrix = [[0, 1], [2, 3]]
        mask = [[False, True], [True, False]]  # True to hide, False to show
    matrix = np.where(mask, 0, matrix) if type == 1 else matrix
    plt.rcParams.update({'font.size': 14})
    if type == 1:
        hm = sns.heatmap(matrix, cmap="rocket_r", xticklabels=unique_lengths_lst, yticklabels=unique_lengths_lst,
                     mask=mask, vmin=0, vmax=np.max(matrix), annot=True, fmt="d")
    elif type == 2:
        hm = sns.heatmap(matrix, cmap="rocket_r", xticklabels=unique_lengths_lst, yticklabels=unique_lengths_lst,
                         vmin=0, vmax=np.max(matrix), annot=True, fmt="d")
    else:
        hm = sns.heatmap(matrix, cmap="rocket_r", xticklabels=unique_lengths_lst, yticklabels=unique_lengths_lst,
                         mask=mask, vmin=0, vmax=np.max(matrix))

    ax = hm.axes
    all_x_labels = ax.get_xticks()
    step=2
    if len(all_x_labels) > 7:
        ax.set_xticks(all_x_labels[::step])
        ax.set_yticks(all_x_labels[::step])
    fig = hm.get_figure()
    fig.set_size_inches(6, 5)
    if type == 1 or type == 3:
        plt.xlabel('Length of the potential piRNA #2')
        plt.ylabel('Length of the potential piRNA #1')
        plt.title('Length of the 2 potential piRNAs\nforming the ping-pong signature')
        plt.tight_layout()
        plt.savefig(f"{path}/heatmap.svg")
    elif type == 2:
        plt.xlabel('Potential primary piRNA\n(T in the 1st position)')
        plt.ylabel('Potential secondary piRNA\n(A in the 10th position)')
        plt.title('Length of the 2 potential piRNAs separating \neach pair in primary and secondary based on 1T/10U')
        plt.tight_layout()
        plt.savefig(f"{path}/heatmap_pairs.svg")
    return None



def create_plots(path, name, x, y, x_name, y_name, saving_path):
    """
    Generates and saves individual line plots for provided data.

    Parameters:
        path (str): Data directory.
        name (str): Plot title and file name.
        x (list): X-axis data points.
        y (list): Y-axis data points.
        x_name (str): X-axis label.
        y_name (str): Y-axis label.
        saving_path (str): Path to save plot image.

    Returns:
        None
    """
    plt.rcParams.update({'font.size': 16})
    fig, ax1 = plt.subplots(1, 1, clear=True, num=1)
    fig.set_size_inches(7, 5)
    ax1.plot(x, y, marker='o')
    ax1.axis(xmin=min(x), xmax=max(x))
    if name == "Initial_distribution":
        name = "Initial distribution"
    ax1.set_title(name)
    ax1.set_xlabel(x_name)
    ax1.set_ylabel(y_name)
    max_y = max(y)
    x_to_y = x[y.index(max_y)]
    ax1.grid()
    ax1.text(x_to_y, max_y, str(x_to_y), fontsize=16)
    plt.tight_layout()
    if 'Length distribution of the sequences' in name:
        name2 = "Length_distribution_of_sequences_displaying_10nts_5-to-5_overlaps"
    elif "5'-to-5' overlap" in name:
        name2 = "5-to-5_overlap"
    else:
        name2 = "Initial_distribution"
    if y_name == 'z-score':
        plt.savefig(f"{saving_path}/{name2}_z_score.svg")
    else:
        plt.savefig(f"{saving_path}/{name2}.svg")
        plt.savefig(f"{saving_path}/{name2}.svg")
    return None


def get_mask(strand, df2, df1):
    """
    Generates a mask identifying reads within piRNA clusters based on strand orientation.

    Parameters:
        strand (str): '+' or '-' strand indicator.
        df2 (DataFrame): Read data.
        df1 (DataFrame): Cluster data.

    Returns:
        ndarray: Boolean mask indicating reads within clusters.
    """
    grouped_clusters = df1.groupby(6)
    clusters_same_strand = grouped_clusters.get_group(strand)
    start = np.asarray(clusters_same_strand[3])
    stop = np.asarray(clusters_same_strand[4])
    reads_same_strand = df2[df2[6] == strand]
    read_start = np.asarray(reads_same_strand[3]).reshape(-1, 1)
    read_stop = np.asarray(reads_same_strand[4]).reshape(-1, 1)
    within_start = (read_start >= start) & (read_stop <= stop)
    within_cluster = np.any(within_start, axis=1)
    return within_cluster

def generete_venn(prim, sec, clusters, path, status):
    """
    Creates a Venn diagram illustrating overlap between piRNA clusters and ping-pong piRNAs.

    Parameters:
        prim (str): Primary piRNA file path.
        sec (str): Secondary piRNA file path.
        clusters (str): Cluster data file path.
        path (str): Output path for Venn diagram.
        status (bool): Indicator for data processing method.

    Returns:
        dict or None: Dictionary of piRNA pair counts or None if no pairs exist.
    """

    sets = Counter()
    try:
        a = pd.read_csv(prim, sep='\t', header=None)
    except EmptyDataError:
        a = np.array([0])
    if status:
        b = pd.read_csv(sec, sep='\t', header=None)
        c = pd.read_csv(clusters, sep='\t', header=None)
        d = pd.read_csv(sec, sep='\t', header=None)
        mask_plus = get_mask('+', a, c)
        mask_minus = get_mask('-', a, c)
        mask = np.concatenate((mask_plus, mask_minus))
        #a = a[mask]
        b[9] = b[8].str.extract('PingPongPair=(.*);')
        b[8] = b[8].str.extract('Seq=(.*)')
        d[9] = d[8].str.extract('ID=(.*);Pi')
        d[8] = d[8].str.extract('Seq=(.*)')
        merged = pd.merge(a, b, on=[8])
        merged = merged.drop_duplicates(subset=['3_y', '4_y', 8, '0_y'])
        sets['10'] = f"piRNAs from clusters (primary piRNAs) (n = {a.shape[0]})"
        sets['01'] = f"piRNAs in ping-pong (n = {b.shape[0]})"
        sets['11'] = f"Intersection (n = {merged.shape[0]})"
        dic_of_pairs = {}
        for i in range(0, len(d), 2):
            first = d.iloc[i, 8]
            second = d.iloc[i + 1, 8]
            which_first = ''
            if first[0] == 'T' and second[9] == 'A':
                text = f"{len(first)}{len(second)}"
                if text in dic_of_pairs.keys():
                    dic_of_pairs[text] += 1
                else:
                    dic_of_pairs[text] = 1
            elif first[9] == 'A' and second[0] == 'T':
                text = f"{len(second)}{len(first)}"
                if text in dic_of_pairs.keys():
                    dic_of_pairs[text] += 1
                else:
                    dic_of_pairs[text] = 1
    else:
        b = np.array([0])
        merged = np.array([0])
        sets['10'] = f"piRNAs from clusters (primary piRNAs) (n = {a.shape[0] - 1})"
        sets['01'] = f"piRNAs in ping-pong (n = {b.shape[0] - 1})"
        sets['11'] = f"Intersection (n = {merged.shape[0] - 1})"
    setlabels = ['piRNAs from clusters (primary piRNAs)', 'piRNAs in ping-pong']
    plt.figure(figsize=(8,5))
    ax = plt.gca()
    v = venn2(subsets = (a.shape[0] - merged.shape[0], b.shape[0] - merged.shape[0], merged.shape[0]),
              set_labels = setlabels, ax=ax)
    h, l = [], []

    for i in sets:
        label = v.get_label_by_id(i)
        if label is not None:
            label.set_text("")

        handle = v.get_patch_by_id(i)
        if handle is None:
            # Create a proxy patch if you need a legend entry anyway
            handle = mpatches.Patch(color='gray', alpha=0.5)
        h.append(handle)
        l.append(sets[i])
    ax.legend(handles=h, labels=l, title="counts", bbox_to_anchor=(1.0, 0.9))
    plt.title("Venn diagram")
    plt.savefig(f"{path}/venn.svg", bbox_inches='tight')
    return dic_of_pairs if status else None


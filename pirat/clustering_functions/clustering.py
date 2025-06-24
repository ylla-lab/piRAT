import pirat.clustering_functions.data_preparation as data_prep
import multiprocessing
from typing import Tuple, List
import pandas as pd
import numpy as np
import os

def multiprocess_clustering(list_of_scaffolds: List[List],
                            k: int,
                            eps: int,
                            range_of_size: List[int],
                            raw_path_name: str,
                            threads: int,
                            variation_threshold: int,
                            saving_path_primary: str) -> List:
    """
    Distribute clustering tasks across multiple CPU cores using a multiprocessing pool.

    This function prepares a list of parameters for each scaffold and then uses Python's `multiprocessing.Pool`
    to run the clustering process in parallel. Each scaffold is processed by the `run_the_clustering` function.

    Parameters:
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
\
    """

    # Print the parameters being used for clustering, for user reference.
    print(f"\nPerforming cluster finding with parameters:"
          f"\n\trange of size of piRNAs: {range_of_size},"
          f"\n\tk: {k},"
          f"\n\teps: {eps}"
          f"\n\tvariation threshold: {variation_threshold}\n")

    # Prepare a list of tuples, where each tuple contains the parameters for processing a single scaffold.
    variables_for_clustering = [
        (scaffold, k, eps, range_of_size, raw_path_name, variation_threshold, saving_path_primary)
        for scaffold in list_of_scaffolds
    ]

    # Create a pool of worker processes.
    with multiprocessing.Pool(threads) as pool:
        # Distribute the work among the worker processes using starmap, which unpacks the parameters for each call.
        results = pool.starmap(run_the_clustering, variables_for_clustering)

    return results


def run_the_clustering(which_scaffold: str,
                       k: int,
                       eps: int,
                       range_of_size: Tuple[int, int],
                       path: str,
                       variation_threshold: int,
                       saving_path_primary: str) -> List:
    """
    Perform clustering of piRNA reads for a given scaffold and write results to GFF files.

    Parameters:
        which_scaffold (str): Identifier of the scaffold to be processed.
        k (int): Min number of reads required to form a cluster.
        eps (int): Max allowed distance between reads for them to be considered part of the same cluster.
        range_of_size (Tuple[int,int]): (min_size, max_size) for piRNAs.
        path (str): Directory containing BAM files.
        variation_threshold (int): Threshold to filter reads by variation.
        saving_path_primary (str): Directory to save output GFF files.

    Returns:
        List: [scaffold_identifier, list_of_clusters]
    """
    # Get all reads for this scaffold
    raw_reads = get_all_reads_for_scaffold(path, which_scaffold, range_of_size, variation_threshold)

    # Remove duplicates and prepare reads for clustering
    processed_reads = data_prep.check_for_duplicates_in_reads_combined(
        raw_reads,
        variation_threshold,  # using variation_threshold as clustering_epsilone
        range_of_size,
        return_length_only=True
    )
    # Cluster the reads
    clusters = finding_cluster(k, eps, processed_reads)

    # Write the clustered reads to GFF
    # First, prepare DataFrame for all reads
    df_reads = prepare_dataframe_for_gff(raw_reads, which_scaffold)

    # Assign clusters and write to GFF file
    clustered_reads = assign_clusters_to_reads(df_reads, clusters)
    write_gff(f'{saving_path_primary}/gffs/Clusters_reads.gff', clustered_reads)

    # Check duplicates again for potential piRNA reads and prepare DataFrame
    potential_pirna_reads = data_prep.check_for_duplicates_in_reads_combined(
        raw_reads,
        variation_threshold,  # using variation_threshold as clustering_epsilone
        range_of_size,
        return_length_only=False
    )
    df_pirna = prepare_dataframe_for_gff(potential_pirna_reads, which_scaffold)

    # Assign clusters and write to another GFF file
    clustered_pirna_reads = assign_clusters_to_reads(df_pirna, clusters)
    write_gff(f'{saving_path_primary}/gffs/potential_pirna_reads.gff', clustered_pirna_reads)

    return [which_scaffold, clusters]


def get_all_reads_for_scaffold(path: str, which_scaffold: str, range_of_size: Tuple[int, int],
                               variation_threshold: int) -> List[List]:
    """
    Retrieve and consolidate reads from all BAM files for a given scaffold.
    """
    bam_files = data_prep.get_list_of_bam_files(path)
    all_reads = []

    for bam_file in bam_files:
        full_path = os.path.join(path, bam_file)
        reads, _ = data_prep.get_reads(full_path, which_scaffold, range_of_size, variation_threshold, 1)
        all_reads.append(reads)

    # Flatten lists
    pos_reads = [item for sublist in [x[0] for x in all_reads] for item in sublist]
    neg_reads = [item for sublist in [x[1] for x in all_reads] for item in sublist]

    # Sort by start position
    pos_reads = sorted(pos_reads, key=lambda x: x[0])
    neg_reads = sorted(neg_reads, key=lambda x: x[0])

    return [pos_reads, neg_reads]


def prepare_dataframe_for_gff(reads: List[List], scaffold: str) -> pd.DataFrame:
    """
    Convert read data into a DataFrame suitable for GFF output.

    The input 'reads' can be:
    - A nested list structure such as [[(start, stop, seq), (start, stop, seq), ...], [(start, stop, seq), ...]]
      or it could already be a single list of tuples.

    Expected formats:
    1. Five-column data:
       (piRNA_annotation, strand, sequence, start_position, stop_position)

       Example of one read:
       ("piRNA_annotation", "-", "ACGT...", 100, 128)

    2. Three-column data:
       (start_position, stop_position, sequence)

       Example of one read:
       (100, 128, "ACGT...")

    If only three columns are provided, the function will add missing columns:
    - piRNA_annotation: "piRNA_annotation"
    - strand: "-" (or "+"; adjust as needed)
    - Read: "Read"

    After creating the DataFrame, this function adds additional columns required for GFF output
    and arranges them in a standard GFF-compatible order.

    Parameters:
        reads (List[List]): The read data, possibly nested.
        scaffold (str): The scaffold name.

    Returns:
        pd.DataFrame: A DataFrame containing columns suitable for GFF output.
    """

    # Flatten the input if it's nested
    # If reads looks like [[(start, stop, seq), ...], [(start, stop, seq), ...]]
    # flatten it to [(start, stop, seq), ...]
    if len(reads) > 0 and isinstance(reads[0], list):
        combined_reads = [r for sublist in reads for r in sublist]
    else:
        combined_reads = reads

    if not combined_reads:
        # If no reads are provided, return an empty DataFrame with the expected columns
        columns = ['Scaffold', 'piRNA_annotation', 'Read', 'start_position', 'stop_position', '.', 'strand', '0',
                   'sequence']
        return pd.DataFrame(columns=columns)

    # Determine how many columns the combined reads have by checking the first element
    num_columns = len(combined_reads[0])

    if num_columns == 5:
        # Data should be in the format: (piRNA_annotation, strand, sequence, start_position, stop_position)
        #df = pd.DataFrame(combined_reads,
        #                  columns=['piRNA_annotation', 'strand', 'sequence', 'start_position', 'stop_position'])
        df = pd.DataFrame(combined_reads,
                          columns=['start_position', 'stop_position', 'sequence', 'strand', 'piRNA_annotation'])
        df['strand'] = np.where(df['strand'] == True, '-', '+')

    elif num_columns == 3:
        # Data is in the format: (start_position, stop_position, sequence)
        # Add missing columns with default values.
        # Adjust these defaults based on your needs.
        processed_data = []
        for (start, stop, seq) in combined_reads:
            # You can pick a default strand or infer it from context. Let's say '-'
            # If you need strand differentiation, pass it as a parameter or handle it elsewhere.
            processed_data.append(('piRNA_annotation', '-', seq, start, stop))

        df = pd.DataFrame(processed_data,
                          columns=['piRNA_annotation', 'strand', 'sequence', 'start_position', 'stop_position'])


    # Add universal columns
    df['Scaffold'] = scaffold
    df['Read'] = 'Read'
    df['.'] = '.'
    df['0'] = '.'

    # Reorder columns to standard GFF fields:
    # GFF typically has 9 columns: seqid, source, type, start, end, score, strand, phase, attributes
    # Here we map them as:
    # Scaffold -> seqid
    # piRNA_annotation -> source
    # Read -> type
    # start_position -> start
    # stop_position -> end
    # '.' -> score
    # strand -> strand
    # '0' -> phase
    # sequence -> attributes (or an attribute field)
    columns_order = ['Scaffold', 'piRNA_annotation', 'Read', 'start_position', 'stop_position', '.', 'strand', '0',
                     'sequence']
    # Make sure all these columns exist in df
    df = df[[col for col in columns_order if col in df.columns]]

    return df


def assign_clusters_to_reads(df: pd.DataFrame, clusters: List[List[int]]) -> List[List[str]]:
    """
    Assign cluster IDs to reads in the DataFrame based on cluster ranges.
    Returns a list of lists representing rows suitable for GFF output.
    """
    cluster_id = 1
    clustered_rows = []

    for cluster in clusters:
        cluster_reads = df[(df['start_position'] >= cluster[0]) & (df['stop_position'] <= cluster[1])].copy()
        cluster_reads.loc[:, 'which_cluster'] = cluster_id
        clustered_rows.extend(cluster_reads.values.tolist())
        cluster_id += 1

    return clustered_rows


def write_gff(filepath: str, rows: List[List]) -> None:
    """
    Write the given rows to a GFF file at the specified filepath.
    Each row should be an array of fields corresponding to GFF columns.
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'a+') as gff_file:
        for row in rows:
            # Assuming the last column might contain which_cluster or other fields,
            # write only the first 9 columns of GFF and ignore extra columns if they exist.
            gff_file.write("\t".join(map(str, row[0:9])) + "\n")


def finding_cluster(k: int,
                    eps: int,
                    list_of_reads: List[List[Tuple[int, int, int]]]) -> List[List[int]]:
    """
    Execute the clustering function for each strand and combine results.

    Parameters:
        k (int): Minimum number of reads within a specified distance (eps) to form a cluster.
        eps (int): Maximum allowed distance between reads.
        list_of_reads (List[List[Tuple[int,int,int]]]): A list of two sublists of reads:
            - list_of_reads[0]: Reads from the positive strand
            - list_of_reads[1]: Reads from the reverse strand
          Each read is a tuple: (read_length, start_position, stop_position).

    Returns:
        List[List[int]]: A combined list of clusters for both strands.
                         Each cluster is represented by [start_position, end_position, strand].
    """
    positive_clusters = localize_the_cluster(k, eps, list_of_reads[0], '-')
    reverse_clusters = localize_the_cluster(k, eps, list_of_reads[1], '+')
    return positive_clusters + reverse_clusters


def localize_the_cluster(k: int,
                         eps: int,
                         sublist_of_reads: List[Tuple[int, int, int]],
                         strand: str) -> List[List[int]]:
    """
    Identify clusters of reads using a modified DBSCAN-like approach. Clusters are formed
    by consecutive reads that are close to each other in genomic coordinates.

    Parameters:
        k (int): Minimum number of reads required to form a cluster.
        eps (int): Maximum allowed distance to group reads into the same cluster.
        sublist_of_reads (List[Tuple[int,int,int]]): List of reads (length, start, stop).
                                                     This list will be sorted by stop position.
        strand (str): '+' or '-' indicating the strand orientation of these reads.

    Returns:
        List[List[int]]: A list of clusters, each represented by [start_position, end_position, strand].
    """
    # Sort reads by their stop position (x[2])
    sublist_of_reads = sorted(sublist_of_reads, key=lambda x: x[2])

    clusters = []
    current_cluster_start = 0
    current_cluster_end = -eps
    current_index = 0
    total_reads = len(sublist_of_reads)

    while current_index <= total_reads - k - 1:
        if can_start_cluster(sublist_of_reads, current_index, k, eps, current_cluster_start):
            # Initialize a new cluster
            current_cluster_start = sublist_of_reads[current_index][1]
            current_cluster_end = sublist_of_reads[current_index + k][2]
            current_index += k  # Jump ahead by k to skip over the cluster core
        elif can_extend_cluster(sublist_of_reads, current_index, k, eps, current_cluster_start, current_cluster_end):
            # Extend the current cluster
            # Check if a bigger cluster can be formed by including read at current_index + k
            if current_index + k < total_reads and \
               (sublist_of_reads[current_index + k][2] - sublist_of_reads[current_index][1] <= eps):
                # Extend by k reads if they fit within eps
                current_cluster_end = sublist_of_reads[current_index + k][2]
                current_index += k
            else:
                # Extend cluster end by the current read stop position
                current_cluster_end = sublist_of_reads[current_index][2]
        else:
            # Current read does not continue or start a cluster
            # If we are currently in a cluster, check if it's large enough to finalize
            if current_cluster_start != 0 and (current_cluster_end - current_cluster_start > eps / 2):
                clusters.append([current_cluster_start, current_cluster_end, strand])
            # Reset cluster start to indicate no active cluster
            current_cluster_start = 0

        current_index += 1

    # After the loop, check if there's an ongoing cluster that meets the size condition
    if current_cluster_start != 0 and (current_cluster_end - current_cluster_start > eps / 2):
        clusters.append([current_cluster_start, current_cluster_end, strand])

    return clusters


def can_start_cluster(reads: List[Tuple[int, int, int]],
                      idx: int,
                      k: int,
                      eps: int,
                      current_start: int) -> bool:
    """
    Determine if a new cluster can start at the given index.

    Conditions:
    - No current cluster is active (current_start == 0).
    - The set of k reads starting at idx fits within eps distance:
      reads[idx + k].stop - reads[idx].start <= eps
    """
    if current_start != 0:
        return False
    if idx + k >= len(reads):
        return False
    return (reads[idx + k][2] - reads[idx][1]) <= eps


def can_extend_cluster(reads: List[Tuple[int, int, int]],
                       idx: int,
                       k: int,
                       eps: int,
                       current_start: int,
                       current_end: int) -> bool:
    """
    Check if the current read at idx can be part of the existing cluster.
    There are several conditions from the original code:
    - The distance between reads[idx].stop and reads[idx - k].start <= eps
      (ensuring reads behind still form a cluster)
    - Or reads[idx + k].stop - reads[idx].start <= eps (if we can skip forward by k and still be within eps)
    - Or reads[idx].stop - current_end <= eps (the read is close enough to the current cluster end)

    Before checking idx - k or idx + k, ensure these indices are valid.
    """
    if current_start == 0:
        return False  # No active cluster to extend

    total = len(reads)
    read_stop = reads[idx][2]
    # Attempt to check reads behind
    can_check_behind = (idx - k) >= 0
    # Attempt to check reads ahead
    can_check_ahead = (idx + k) < total

    behind_condition = can_check_behind and ((read_stop - reads[idx - k][1]) <= eps)
    ahead_condition = can_check_ahead and ((reads[idx + k][2] - reads[idx][1]) <= eps)
    end_extension_condition = (read_stop - current_end) <= eps

    return behind_condition or ahead_condition or end_extension_condition
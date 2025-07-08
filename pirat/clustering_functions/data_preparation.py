import os
import pysam
from pathlib import Path
import logging
from typing import List, Tuple
from Bio.Seq import Seq
import pirat.ping_pong.ping_pong_all as pin_pon

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def get_list_of_bam_files(raw_path_name: str) -> List[str]:
    """
    Retrieve a list of BAM files from a given directory. If the corresponding index files
    (.bai or .csi) are missing, they are created automatically.

    Parameters:
        raw_path_name (str): The path to the directory containing BAM files.

    Returns:
        List[str]: A list of BAM filenames found in the given directory.
    """
    raw_path = Path(raw_path_name)
    if not raw_path.is_dir():
        logging.error(f"{raw_path_name} is not a valid directory.")
        return []

    bam_files = []
    all_files = os.listdir(raw_path)

    for data_file in all_files:
        if data_file.lower().endswith(".bam"):
            bam_files.append(data_file)
            # Check for index files
            bai_exists = f"{data_file}.bai" in all_files
            csi_exists = f"{data_file}.csi" in all_files
            if not (bai_exists or csi_exists):
                # Create index if missing
                try:
                    pysam.index(str(raw_path / data_file))
                except Exception as e:
                    sys.exit(f"Failed to index BAM file {data_file} in {raw_path_name}: {e}")
    return bam_files

def get_list_of_scaffolds(path_to_file: str) -> List[str]:
    """
    Obtain a list of scaffold names from the given BAM file.

    Parameters:
        path_to_file (str): The path to a BAM file.

    Returns:
        List[str]: A list of scaffold names extracted from the BAM file.

    """
    try:
        idxstats_result = pysam.idxstats(path_to_file)
    except Exception as e:
        logging.error(f"Error reading idxstats from {path_to_file}: {e}")
        return []

    lines = idxstats_result.strip().split("\n")
    formatted = [line.split("\t") for line in lines if line]
    # The last line might be empty or contain a reference with zero length
    # Typically idxstats output ends with an asterisk line indicating unmapped reads
    # We'll exclude it if it's there.
    scaffolds = [sc[0] for sc in formatted if len(sc) > 1 and sc[0] != '*']
    return scaffolds


def check_for_duplicates_in_reads_combined(
        list_of_reads_for_scaffold: List[List],
        variation_threshold: int,
        range_of_size: List[int],
        return_length_only: bool = True
) -> List[List[Tuple]]:
    """
    Cleans and filters reads for a given scaffold, checking for duplicates and removing
    those that don't fit the desired size range. The output structure and sorting criteria
    depend on the `return_length_only` flag.

    Parameters:
        list_of_reads_for_scaffold (List[List]): Two lists of reads,
            one for reverse and one for forward:
            - list_of_reads_for_scaffold[0] are the reverse reads
            - list_of_reads_for_scaffold[1] are the forward reads
        variation_threshold (int): Distance threshold for clustering.
        range_of_size (List[int]): [min_size, max_size] range of allowed read sizes.
        return_length_only (bool):
            If True, returns data as (read_length, start, end) tuples and sorts by start position.
            If False, returns data as (start, end, sequence) tuples and sorts by start position.

    Returns:
        Filtered list containing sublist of reads from reverse and positive strand
    """

    # Extract reads for reverse and forward
    rev_reads = [(int(r[0]), int(r[1]), r[2]) for r in list_of_reads_for_scaffold[0]]
    pos_reads = [(int(r[0]), int(r[1]), r[2]) for r in list_of_reads_for_scaffold[1]]

    # Get read abundances
    rev_abundance = pin_pon.get_abundance(rev_reads)
    pos_abundance = pin_pon.get_abundance(pos_reads)

    # Remove duplicates by using set
    rev_reads = list(set(rev_reads))
    pos_reads = list(set(pos_reads))

    # Cleanse reads using pin_pon function
    rev_reads, _ = pin_pon.cleanse_reads(rev_reads, rev_abundance, variation_threshold, range_of_size)
    pos_reads, _ = pin_pon.cleanse_reads(pos_reads, pos_abundance, variation_threshold, range_of_size)

    # Filter by size
    rev_filtered = [r for r in rev_reads if range_of_size[0] <= len(r[2]) <= range_of_size[1]]
    pos_filtered = [r for r in pos_reads if range_of_size[0] <= len(r[2]) <= range_of_size[1]]

    # If return_length_only is True, transform tuples to (length, start, end)
    # and sort by start (which is at index 1 in this transformed tuple).
    if return_length_only:
        rev_filtered = [(len(r[2]), r[0], r[1]) for r in rev_filtered]
        pos_filtered = [(len(r[2]), r[0], r[1]) for r in pos_filtered]

        rev_filtered.sort(key=lambda x: x[1])  # sort by start position
        pos_filtered.sort(key=lambda x: x[1])  # sort by start position

    else:
        # If return_length_only is False, keep (start, end, sequence) and sort by start (index 0).
        # We already have (start, end, seq) in rev_filtered and pos_filtered.
        # Just confirm structure and then sort:
        # They are currently (start, end, sequence), no transformation needed.
        rev_filtered.sort(key=lambda x: x[0])  # sort by start position
        pos_filtered.sort(key=lambda x: x[0])  # sort by start position

    return [rev_filtered, pos_filtered]

def get_reads(file: str, one_scaffold: str, range_of_size: Tuple[int, int], variation_threshold: int, flag: int) -> List[List]:
    """
    Load reads from datafiles in specific way based on part of analysis they are gonna be used for.
    Making sure that reference length is not Null
    If reads are used for clustering/ping-pong annotation: flag 1 - load the reads from defined size of range of piRNAs extended
    by variation_threshold parameter, additionally, if they are coming from reverse strand - rev compl them
    If reads are used for sample wise analysis - flag: 0 - wide range for initial analysis of individual file
    At the end save information about every read for read length distribution

    Parameters:
       file (str) - path to the file to be checked if contains paired reads
       one_scaffold (str) - scaffold name of reads to be loaded from
       variation_threshold (int) - Threshold value used to filter or determine acceptable variability across reads with the same core.
       range_of_size: A list specifying the minimum and maximum sizes of piRNAs ([min_size, max_size]).
       flag (int) - in what way the reads should be loaded
    Returns:
        List containing sublist of reads from reverse and positive strand

    """
    length_of_all_reads = []
    list_of_proper_reads = [[], []]
    paired = check_if_paired(file)
    for row in pysam.AlignmentFile(file, 'rb').fetch(one_scaffold):
        if paired:
            if row.is_paired and row.is_read1:
                if row.reference_length is not None and row.reference_length != 0:
                    if flag == 1:
                        if len(row.seq) == row.reference_length and range_of_size[0] - variation_threshold <= len(row.seq) <= range_of_size[1] + variation_threshold:
                            reverse = row.is_reverse
                            sequence = row.seq if row.is_reverse is False else (Seq(row.seq)).reverse_complement()
                            read_info = (row.pos + 1, row.reference_end, sequence, row.is_reverse, row.reference_name)
                            if reverse:
                                list_of_proper_reads[0].append(read_info)
                            else:
                                list_of_proper_reads[1].append(read_info)
                    elif flag == 0:
                        if len(row.seq) == row.reference_length and 20 <= len(row.seq) <= 40:
                            reverse = row.is_reverse
                            sequence = row.seq if row.is_reverse is False else (Seq(row.seq)).reverse_complement()
                            read_info = (row.pos + 1, row.reference_end, sequence, row.is_reverse, row.reference_name)
                            if reverse:
                                list_of_proper_reads[0].append(read_info)
                            else:
                                list_of_proper_reads[1].append(read_info)
            elif row.is_paired is False:
                if row.reference_length is not None and row.reference_length != 0:
                    if flag == 1:
                        if len(row.seq) == row.reference_length and range_of_size[0] - variation_threshold <= len(row.seq) <= range_of_size[1] + variation_threshold:
                            reverse = row.is_reverse
                            sequence = row.seq if row.is_reverse is False else (Seq(row.seq)).reverse_complement()
                            read_info = (row.pos + 1, row.reference_end, sequence, row.is_reverse, row.reference_name)
                            if reverse:
                                list_of_proper_reads[0].append(read_info)
                            else:
                                list_of_proper_reads[1].append(read_info)
                    elif flag == 0:
                        if len(row.seq) == row.reference_length and 20 <= len(row.seq) <= 40:
                            reverse = row.is_reverse
                            sequence = row.seq if row.is_reverse is False else (Seq(row.seq)).reverse_complement()
                            read_info = (row.pos + 1, row.reference_end, sequence, row.is_reverse, row.reference_name)
                            if reverse:
                                list_of_proper_reads[0].append(read_info)
                            else:
                                list_of_proper_reads[1].append(read_info)
            else:
                continue
        else:
                if row.reference_length is not None and row.reference_length != 0:
                    if flag == 1:
                        if len(row.seq) == row.reference_length and range_of_size[0] - variation_threshold <= len(row.seq) <= range_of_size[1] + variation_threshold:
                            reverse = row.is_reverse
                            sequence = row.seq if row.is_reverse is False else (Seq(row.seq)).reverse_complement()
                            read_info = (row.pos + 1, row.reference_end, sequence, row.is_reverse, row.reference_name)
                            if reverse:
                                list_of_proper_reads[0].append(read_info)
                            else:
                                list_of_proper_reads[1].append(read_info)
                    elif flag == 0:
                        if len(row.seq) == row.reference_length and 20 <= len(row.seq) <= 40:
                            reverse = row.is_reverse
                            sequence = row.seq if row.is_reverse is False else (Seq(row.seq)).reverse_complement()
                            read_info = (row.pos + 1, row.reference_end, sequence, row.is_reverse, row.reference_name)
                            if reverse:
                                list_of_proper_reads[0].append(read_info)
                            else:
                                list_of_proper_reads[1].append(read_info)

        length_of_all_reads += [len(row.seq)]

    return [list_of_proper_reads, length_of_all_reads]


def check_if_paired(file: str) -> bool:
    """
    Check if file is containing paired reads.

    Parameters:
        file (str) - path to the file to be checked if contains paired reads

    Returns:
        bool - True if paired reads are present, False otherwise
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
    bam.close()
    return paired_reads / total_reads > 0.5

def list_of_reads_for_finding_parameters(file: str, one_scaffold: str, variation_threshold: int, range_of_size: List[int]) -> List[List]:
    """
    Loading and processing the reads for cluster parameter finding.
    First loading the reads for given scaffold using get_reads function, then finding out abundance of the reads, and then cleansing based
    on variation threshold parameter

    Parameters:
        file (str) - path to the file to be checked if contains paired reads
        one_scaffold (str) - scaffold name of reads to be loaded from
        variation_threshold (int) - Threshold value used to filter or determine acceptable variability across reads with the same core.
        range_of_size: A list specifying the minimum and maximum sizes of piRNAs ([min_size, max_size]).

    Returns:
        List containing sublist of reads from reverse and positive strand
    """

    list_of_read_of_a_scaffold, _ = get_reads(file, one_scaffold, (range_of_size[0] - variation_threshold,range_of_size[1] + variation_threshold),
                                              variation_threshold=variation_threshold, flag=1)
    list_for_checking_rev = [(int(read[0]), int(read[1]), read[2]) for read in
                             list_of_read_of_a_scaffold[0]]
    list_for_checking_pos = [(int(read[0]), int(read[1]), read[2]) for read in
                             list_of_read_of_a_scaffold[1]]
    rev_abundance = pin_pon.get_abundance(list_for_checking_rev)
    pos_abundance = pin_pon.get_abundance(list_for_checking_pos)
    list_for_checking_rev = list(set(list_for_checking_rev))
    list_for_checking_pos = list(set(list_for_checking_pos))
    list_for_checking_rev.sort(key=lambda read: read[1])
    list_for_checking_pos.sort(key=lambda read: read[1])
    list_for_checking_rev = pin_pon.cleanse_reads(list_for_checking_rev, rev_abundance, epsilon=variation_threshold,
                                                  range_of_size =range_of_size)
    list_for_checking_pos = pin_pon.cleanse_reads(list_for_checking_pos, pos_abundance, epsilon=variation_threshold,
                                                  range_of_size=range_of_size)
    list_for_checking_rev = [read for read in list_for_checking_rev[0] \
                             if range_of_size[0] <= len(read[2]) <= range_of_size[1]]
    list_for_checking_pos = [read for read in list_for_checking_pos[0] \
                             if range_of_size[0] <= len(read[2]) <= range_of_size[1]]
    return [list_for_checking_rev, list_for_checking_pos]

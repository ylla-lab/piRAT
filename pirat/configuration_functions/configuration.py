import sys
import os
import pirat.configuration_functions.user_input as u_ip

def get_threads() -> int:
    """
    Determines and returns the number of available CPU threads.

    Returns:
        int: Number of CPU threads available.
    """
    if sys.platform == 'win32':
        return int(os.environ['NUMBER_OF_PROCESSORS'])
    else:
        return int(os.popen('grep -c cores /proc/cpuinfo').read())

def create_dirs(raw_path_name: str, opt: int, automate: bool, draw_plots: bool) -> None:
    """
        Creates necessary directory structures for analysis results and plots, and removes existing result files if needed.

        Parameters:
            raw_path_name (str): Path to the directory where subdirectories will be created.
            opt (int): Operation mode option determining directory/file cleanup behavior.
            automate (bool): Flag indicating whether file removal confirmation should be automated.
            draw_plots (bool): Flag indicating whether directories for plots should be created.

        Returns:
            None
    """
    list_of_files = os.listdir(raw_path_name)
    if 'reports' not in list_of_files:
        os.mkdir(f"{raw_path_name}/reports")
    if 'plots' not in list_of_files:
        os.mkdir(f"{raw_path_name}/plots")
    if 'gffs' not in list_of_files:
        os.mkdir(f"{raw_path_name}/gffs")
    if 'sample_wise_analysis' not in list_of_files:
        os.mkdir(f"{raw_path_name}/sample_wise_analysis")
    if 'other_data' not in list_of_files:
        os.mkdir(f"{raw_path_name}/other_data")
    list_of_files = os.listdir(raw_path_name)
    if draw_plots and "Clusters" not in os.listdir(f"{raw_path_name}/plots/"):
        os.mkdir(f"{raw_path_name}/plots/Clusters")
    if 'gffs' in list_of_files:
        list_of_results = os.listdir(f"{raw_path_name}/gffs")
    if ("Clusters_reads.gff" in list_of_results or "potential_pirna_reads.gff"  in list_of_results) \
            and (opt == 1 or opt == 3):
        if automate is False:
            proceed = u_ip.delete_existing_gff(raw_path_name)
        else:
            proceed = True
        if proceed is True and "Clusters_reads.gff" in os.listdir(f"{raw_path_name}/gffs"):
            os.remove(f"{raw_path_name}/gffs/Clusters_reads.gff")
        if proceed is True and "potential_pirna_reads.gff" in os.listdir(f"{raw_path_name}/gffs"):
            os.remove(f"{raw_path_name}/gffs/potential_pirna_reads.gff")

def create_dir_ping_pong(raw_path_name, name_of_the_file):
    """
        Creates a dedicated directory for storing ping-pong analysis results for a specific sample.

        Parameters:
            raw_path_name (str): Base directory for storing results.
            name_of_the_file (str): Name of the sample-specific directory to create.

        Returns:
            None
    """
    list_of_files = os.listdir(f"{raw_path_name}/sample_wise_analysis")
    if name_of_the_file not in list_of_files:
        os.mkdir(f"{raw_path_name}/sample_wise_analysis/{name_of_the_file}")
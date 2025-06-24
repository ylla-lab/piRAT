import matplotlib.pyplot as plt
from typing import List


def save_raw_data(data_lst: list, labels: List[str], file_name: str) -> None:
    """
        Saves provided datasets to individual text files with standardized filenames.

        Parameters:
            data_lst (list): List of pandas DataFrames containing data to be saved.
            labels (List[str]): List of labels to be used as filenames.
            file_name (str): Common suffix for output filenames.

        Returns:
            None
    """
    i = 0
    for label, data in zip(labels, data_lst):
        labels[i] = labels[i].replace(' ', '_')
        labels[i] = labels[i].replace("'", "")
        data.to_csv(labels[i] + file_name + '.txt', sep='\t', index=False)
        i += 1


def draw_scope_plot(fig_name: str, labels: List, x: List, y: List, low_lim: int, up_lim: int, raw_path_name: str, saving_path_primary: str) \
        -> None:
    """
        Generates and saves a bar plot indicating sequence length distribution with marked range boundaries.

        Parameters:
            fig_name (str): Filename for saving the plot.
            labels (List[str]): Labels for y-axis and x-axis respectively.
            x (List): X-axis values (sequence lengths).
            y (List): Corresponding Y-axis values (frequencies).
            low_lim (int): Lower limit of the highlighted range.
            up_lim (int): Upper limit of the highlighted range.
            raw_path_name (str): Base path for raw data (currently unused).
            saving_path_primary (str): Directory to save the generated plot.

        Returns:
            None
    """
    plt.rcParams["figure.figsize"] = (12, 8)
    plt.rcParams.update({'font.size': 26})
    plt.rcParams["axes.linewidth"] = 3
    plt.bar(x, y, align='center')
    plt.plot([low_lim - 0.5, low_lim - 0.5], [0, max(y)], color='red')
    plt.plot([up_lim + 0.5, up_lim + 0.5], [0, max(y)], color='red')
    plt.xlim([9.5, 50.5])
    plt.ylabel(labels[0])
    plt.xlabel(labels[1])
    plt.savefig(f"{saving_path_primary}/plots/{fig_name}.svg")
    return None


def draw_parameters_plot(fig_name: str, labels: List, x: List, y: List, low_lim: int, up_lim: int) -> None:
    """
        Creates and saves comparative line plots illustrating parameter dependencies and thresholds.

        Parameters:
            fig_name (str): Filename to save the plot.
            labels (List): Nested list containing axis and plot labels.
                           Format: [[x-axis1, y-axis1], [x-axis2, y-axis2], [legend_labels]]
            x (List): List containing data series for the x-axis.
            y (List): List containing multiple data series for the y-axis.
            low_lim (int): Lower limit for the y-axis.
            up_lim (int): Upper limit for the y-axis.

        Returns:
            None
    """
    
    plt.rcParams["figure.figsize"] = (20, 10)
    plt.rcParams["axes.linewidth"] = 3
    plt.rcParams.update({'font.size': 30})
    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax2.plot(y[0], label=labels[2][0], lw=2)
    ax2.set_xlabel(labels[0][0])
    ax2.set_ylabel(labels[0][1])
    ax1.plot(x[0], y[1], label="Eps", lw=2)
    ax1.plot(x[1], y[2], label=labels[2][1], lw=2)
    ax1.set_xlabel(labels[1][0])
    ax1.set_ylabel(labels[1][1])
    ax1.set_ylim([low_lim, up_lim])
    plt.legend()
    plt.tight_layout(pad=1.2)
    plt.savefig(f"{fig_name}")
    plt.clf()
    del fig
    return None
    
def save_clusters_to_gff(results: List, raw_path_name: str) -> None:
    """
        Writes identified piRNA clusters to GFF file ('clusters_out.gff').

        Parameters:
            results (List): List containing tuples of cluster information (chromosome, positions, and strand).
            raw_path_name (str): Base directory to save the output GFF files.

        Returns:
            None
    """
    which_cluster = 1
    file = open(f"{raw_path_name}/gffs/clusters_out.gff", "w+")
    for x in results:
        a = sorted(x[1], key=lambda k: k[0])
        for i in a:
            file.write(f"{x[0]}\tpiRNA_annotation\tCluster\t{i[0]}\t{i[1]}\t.\t{i[2]}\t.\tCluster_{which_cluster}\n")
            which_cluster += 1
    file.close()


def save_high_low_quality_to_gff(results: List, raw_path_name: str) -> None:
    """
        Saves clusters classified by quality ('high_quality' and 'low_quality') into separate GFF files.

        Parameters:
            results (List): List containing two sublists of cluster data, each for high and low-quality clusters.
            raw_path_name (str): Base directory for saving the quality-classified GFF files.

        Returns:
            None
    """
    it_cl = 0
    for quality in ['high_quality', 'low_quality']:
        file = open(f"{raw_path_name}/gffs/{quality}_clusters_out.gff", "w+")
        for i in results[it_cl]:
            file.write(f"{i[0]}\tpiRNA_annotation\tCluster\t{i[1]}\t{i[2]}\t.\t{i[3]}\t.\tID={i[4]}\n")
        it_cl += 1
        file.close()
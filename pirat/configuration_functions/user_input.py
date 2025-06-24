from typing import Tuple, List, Dict

def define_scope_range(raw_path_name: str, range_of_size: List[int]) -> List[int]:
    """

    Input:

    range_of_size - the found range of the scope

    Output:

    Based on user choice it returns either the updated range_of_size variable or the previous one.

    """
    print(
        f"\nRange of the size of piRNAs in bam file is found to be {range_of_size[0]} to {range_of_size[1]}.\n\nDo you "
        f"want to input your own range based on plot found in {raw_path_name}? [y/n]:")
    adjust_the_range = input()
    if adjust_the_range.lower() == "y":
        print("\nPlease input your range in this format:\tlower_limit,upper_limit\n")
        try:
            range_of_size = [int(el) for el in (input()).split(",")]
            print(f"\nAdjusted range is from {range_of_size[0]} to {range_of_size[1]}")
        except IndexError:
            print("New range is in the wrong format. The default range remains.")
    else:
        pass
    return range_of_size


def define_params(k: int, eps: int, dic_of_eps: Dict, raw_path_name: str) -> Tuple[int, int]:
    list_of_k = [4, 6, 8, 10]
    user_input = ''
    good = False
    if k is False:
        print(f"Automatic detection of clustering parameters failed. \n"
              f"Please input parameters based on plots found in: {raw_path_name}plots/\n"
              f"as k,Eps. e.g: 8,700")
        while good != True:
            params = input()
            try:
                k, eps = [int(x) for x in params.split(',')]
                good = True
            except ValueError:
                print('\nPlease input your parameters as k,Eps. e.g: 8,700')
    else:
        while user_input not in ['y', 'n']:

            print(f"\nFound parameters:\n")
            [print(f"{k_} : {dic_of_eps[k_]}") for k_ in list_of_k]
            print(
                f"\nPlots for above parameters can be found in: {raw_path_name}.\nAs default, chosen parameters "
                f"are k: {k}, Eps: {eps}\nDo you want to change the parameters? [y/n]]\n")
            user_input = input()

            if user_input == 'y':
                print('\nPlease input your parameters as k,Eps. e.g: 8,700')
                params = input()
                k, eps = [int(x) for x in params.split(',')]

    return k, eps

def delete_existing_gff(raw_path_name: str) -> bool:
    user_input = ""
    while user_input != "y":
        print(f"\nThere is existing results file: {raw_path_name}gffs/Cluster_reads.gff. "
              f"It will be deleted before the continuation. Enter 'y' to proceed:")
        user_input = input()
    return True

from extracting_phylogenies.utilities import utilities as ut 
import argparse # argument parsing 
from argparse import RawTextHelpFormatter
import datetime # for error messages and logging
import re # parsing newick string
import nltk # Levenshtein
from ete3 import Tree # turn newick into tree structure
import os # for checking outfile paths of the tsv
import sys
import logging
from extracting_phylogenies.utilities import utilities as ut
from itertools import combinations
from statistics import mean, median

# create logger
console_logger = logging.getLogger(__name__)

def get_newick_from_file(path):
    """
    Given a path to a newick file returns the newick.

    Args:
        path (str): path to newick file

    Raises:
        FileNotFoundError: if file doesnt exist
        IsADirectoryError: if path is a directory path
        ValueError: if Newick is not valid   
    """
    if not os.path.exists(path):
        raise FileNotFoundError("File doesn't exist.")
    elif os.path.isdir(path):
        raise IsADirectoryError("Expected filepath but got directory path.")
    with open(path, "r") as nwk_file:
        newick = nwk_file.read()
    # check if newick is valid 
    try:
        ut.get_newick_format(newick)
        return newick
    except Exception as e:
        raise e 

def get_filename(filepath):
    """
    Given a filepath returns the files name.
    Throws exceptions if path doesnt exist or isnt a filepath. 

    Args:
        filepath (str): path to a file

    Raises:
        IsADirectoryError: If the path is a directory path
        FileNotFoundError: if the path doesnt exist

    Returns:
        str: name of the file
    """
    if os.path.exists(filepath): 
        if os.path.isfile(filepath):
            head, tail = os.path.split(filepath)
            return tail
        else:
            raise IsADirectoryError("Expected a filepath but got a directory path.")
    else:
        raise FileNotFoundError("Path is not valid.")
    
# for comparison of taxa especially those with single character substitutions
# problem: insertions and deletions of characters
def hamming_distance(str1, str2):
    """
    Given two strings of equal length calculates the hamming distance.
    Returns None if strings arent of equal length.

    Args:
        str1 (str): first string
        str2 (str): second string

    Returns:
        int: hamming distance  
    """
    if not len(str1) == len(str2):
        raise ValueError(f"Strings have to be of equal length.")
    hamming_distance = 0
    for index in range(len(str1)):
        if not str1[index] == str2[index]:
            hamming_distance += 1
    return hamming_distance

def edit_distance_ratio(original_taxon, generated_taxon):
    return 1 - (nltk.edit_distance(original_taxon, generated_taxon, substitution_cost=1, transpositions=False)
                /max(len(original_taxon),len(generated_taxon)))
    
def get_taxon_pairs_greedy(original_taxa, generated_taxa):
    """
    Given a dictionary with two lists: mismatched original taxa and mismatched generated taxa,
    returns a dict with each key being an original taxon and the value being the corresponding generated taxon.
    
    Iterates through the original taxa and greedily assigns the best generated taxon to the current original taxon the pair 
    has a edit distance ratio of at least 75%.
    
    If there are more original than generated taxa unpaired original taxa are assigned an empty string.
    If there are more generated than original taxa then generated taxa with lowest similarity are ignored.
    
    Not reliable if amount of original taxa and generated_taxa doesnt match 
    
    Args:
        original_taxa (List(str)): List of original taxa 
        generated_taxa (List(str)): List of generated taxa 
    Returns:
        dict: dict with original_taxon/generated_taxa pairs
    """
    if not (orig:=len(original_taxa)) == (gen:=len(generated_taxa)):
        console_logger.info(
            f"Amount of taxa in original and generated newick do not match. Original: {orig}. Generated: {gen}."
        )
    pair_dict = dict()
    for original_taxon in original_taxa:
        score_dict = dict()
        for generated_taxon in generated_taxa:
            score_dict[generated_taxon] = edit_distance_ratio(original_taxon, generated_taxon)
        if score_dict and max(score_dict.values()) > 0.75:    
            pair_dict[original_taxon] = max(score_dict, key=score_dict.get)
            # remove the generated taxon that was just assigned from the list of generated taxa
            generated_taxa.remove(pair_dict[original_taxon])
        else: 
            pair_dict[original_taxon] = ""
    return pair_dict 

class Comparison_Job():
    def __init__(
        self,
        original_newick_path = None,
        generated_newick_path = None,
        original_newick = None,
        generated_newick= None,
        format_original = None, # remember formats i.e. 0,5,9 or 100
        format_generated = None, 
    ):
        self.original_newick_path = original_newick_path
        self.generated_newick_path = generated_newick_path
        self.original_newick = original_newick
        self.generated_newick = generated_newick
        self.format_original = format_original
        self.format_generated = format_generated
        
    def get_tsv_header(self, info_header):
        """
        Function that returns the tsv header

        Returns:
            str: .tsv header
        """
        tsv_header = "newick1\tnewick2\tformat1\tformat2\t#_taxa1\t#_taxa2\tcorrect_taxa_ratio\t#_equal_length\t"
        tsv_header += "#_unequal_length\tmean_ham_dist\tmean_ham_ratio\tmean_edit_dist\tmean_edit_dist_total\t"
        tsv_header += "mean_edit_ratio\tmean_edit_ratio_total\trf_dist\tmax_rf_dist\trf_ratio\t"
        tsv_header += "#_orig_edges\t#_gen_edges\t#_common_edges\t#_missing_edges\tcorrect_edge_ratio\t"
        if info_header:
            tsv_header += info_header
        # make sure there is only one newline at the end of the header
        tsv_entry = tsv_entry.rstrip("\n") + "\n"
        return tsv_header

    def get_tsv_entry(self, param_entry, taxa_comp, dist_comp, topo_comp):
        tsv_entry = f"{get_filename(self.original_newick_path)}\t{get_filename(self.generated_newick_path)}\t"
        tsv_entry+=f"{self.format_original}\t{self.format_generated}\t"
        if taxa_comp:
            tsv_entry+=f"{taxa_comp["taxa_count_original"]}\t{taxa_comp["taxa_count_generated"]}\t"
            tsv_entry+=f"{taxa_comp["correct_taxa_ratio"]}\t{taxa_comp["count_equal_length"]}\t{taxa_comp["count_unequal_length"]}\t"
            tsv_entry+=f"{taxa_comp["mean_hamming_distance"]}\t{taxa_comp["mean_hamming_distance_ratio"]}\t"
            tsv_entry+=f"{taxa_comp["mean_edit_distance"]}\t{taxa_comp["mean_edit_distance_total"]}\t{taxa_comp["mean_edit_ratio"]}\t"
            tsv_entry+=f"{taxa_comp["mean_edit_ratio_total"]}\t"
        else:
            tsv_entry+="None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t"
        if dist_comp:
            tsv_entry+=f"{dist_comp["mean_abs_diff"]}\t{dist_comp["median_abs_diff"]}\t"
            tsv_entry+=f"{dist_comp["mean_neg_diff"]}\t{dist_comp["median_neg_diff"]}\t"
            tsv_entry+=f"{dist_comp["mean_pos_diff"]}\t{dist_comp["median_pos_diff"]}\t"
            tsv_entry+=f"{dist_comp["mean_pairwise_diff"]}\t{dist_comp["median_pairwise_diff"]}\t"
        else:
            tsv_entry+="None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t"
        if topo_comp:
            tsv_entry+=f"{topo_comp["rf"]}\t{topo_comp["max_rf"]}\t{topo_comp["rf_ratio"]}\t{topo_comp["count_original_edges"]}\t"
            tsv_entry+=f"{topo_comp["count_missing_edges"]}\t{topo_comp["count_common_edges"]}\t"
            tsv_entry+=f"{topo_comp["correct_edges_ratio"]}\t{topo_comp["treeko_dist"]}\t"
            tsv_entry+=f"{topo_comp["count_multifurcations_original"]}\t{topo_comp["count_multifurcations_original"]}"
        else:
            tsv_entry+="None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t"
        if param_entry:
            tsv_entry+=param_entry
        # make sure there is only one newline at the end of the entry
        tsv_entry = tsv_entry.rstrip("\n") + "\n"
        return tsv_entry
        
    def compare_taxa(self):
        """
        Given the strings of two newicks returns a dictionary of statistics of the taxa of both

        Args:
            original (str): original newick that the generated one is compared to
            generated (str): (AI-generated) newick whose quality is calculated

        Returns:
            dict: statistics of how similar the taxa of both newicks are
        """
        # dictionary for the comparison of taxa
        taxa_dict = dict()
        original_taxa = ut.get_taxa(self.original_newick)
        generated_taxa = ut.get_taxa(self.generated_newick)
        # count taxa in both newicks
        original_taxa_count = len(original_taxa)
        generated_taxa_count = len(generated_taxa)
        taxa_dict["taxa_count_original"] = original_taxa_count
        taxa_dict["taxa_count_generated"] = generated_taxa_count
        # percentage of correct taxa, average hamming distance, average levenshtein distance
        match_counter = 0
        taxon_pairs = get_taxon_pairs_greedy(original_taxa, generated_taxa)
        for original, generated in taxon_pairs.items():
            if original == generated:
                match_counter += 1
        taxa_dict["correct_taxa_ratio"] = round(match_counter/len(original_taxa), 4) 
        # count number of taxa with equal length (implies substitutions), unequal length (implies insertions/deletions)
        equal_length_counter = 0
        unequal_length_counter = 0
        # save which original taxa dont have a matching generated taxon
        mismatches = []
        # calculate hamming distance and ratio over all taxa pairs, taxa pairs that cant be compared are maximally 
        # "punished"
        accu_hamming_distances = 0
        accu_hamming_ratios = 0
        for original, generated in taxon_pairs.items():
            if generated: 
                if len(original) == len(generated):
                    hdist = hamming_distance(original, generated)
                    accu_hamming_distances += hdist
                    accu_hamming_ratios += 1 - (hdist/len(original))
                    equal_length_counter += 1
                else:
                    unequal_length_counter += 1
                    # pairs that dont have matching lengths get max hamming distance (length of original newick)
                    accu_hamming_ratios += 0.0
                    accu_hamming_distances += len(original)
            else:
                # original taxon without matching generated one also gets max hamming distance
                accu_hamming_ratios += 0.0
                accu_hamming_distances += len(original)
                mismatches.append(original)
        # save counts of matching length, mismatchingl length and taxon mismatches/missing taxa
        taxa_dict["count_equal_length"] = equal_length_counter
        taxa_dict["count_unequal_length"] = unequal_length_counter
        taxa_dict["mismatches"] = mismatches
        # save mean hamming distance and ratio
        taxa_dict["mean_hamming_distance"] = round(accu_hamming_distances/original_taxa_count, 4)
        taxa_dict["mean_hamming_distance_ratio"] = round(accu_hamming_ratios/original_taxa_count, 4)
        # calculate average edit distance and ratio for all taxa pairs
        accu_edit_distances = 0
        accu_edit_ratios = 0
        for original, generated in taxon_pairs.items():
            edit_distance = nltk.edit_distance(original, generated, substitution_cost=1, transpositions=False)
            accu_edit_distances += edit_distance
            # edit distance ratio = 1 - edit distance/max length => totally different strings get 0.00 
            accu_edit_ratios += 1 - (edit_distance/max(len(original), len(generated)))
        # mean distance calculated only over pairs where the generated taxon isnt "" 
        # => should generally be higher than mean_distance_total
        taxa_dict["mean_edit_distance"] = round(accu_edit_distances/len(taxon_pairs), 4)
        # mean distance calculated over all pairs (even those where the generated taxon is "" => max distance)
        # => should generally be lower than mean_edit_distance
        taxa_dict["mean_edit_distance_total"] = round(accu_edit_distances/(len(taxon_pairs)- len(mismatches)), 4)
        # also calculate a ratio over pairs with generated newick and over all pairs
        taxa_dict["mean_edit_ratio"] = round(accu_edit_ratios/(len(taxon_pairs)- len(mismatches)), 4)
        taxa_dict["mean_edit_ratio_total"] = round(accu_edit_ratios/len(taxon_pairs), 4)
        return taxa_dict  

    def compare_topology(self):
        """
        Given 2 Newick strings compares the topologies of their resulting trees.
        To compare the topologies, taxa of the generated newick are assigned to their counterpart in the original newick.
        However the AI might have spelling mistakes and thus taxa from the generated newick are greedily assigned to the
        taxa of the original newick

        Args:
            original (str): original newick
            generated (str): generated newick
        """
        original = self.original_newick
        generated = self.generated_newick
        # TODO add optimal pairing for edge cases like taxon1, taxon2
        taxa_pairs = get_taxon_pairs_greedy(ut.get_taxa(original),ut.get_taxa(generated))
        # correct spelling mistakes by AI
        for original_taxon, generated_taxon in taxa_pairs.items():
            generated.replace(generated_taxon, original_taxon)
        original_tree = Tree(original)
        generated_tree = Tree(generated)
        original = original_tree.write()
        generated = generated_tree.write()
        # dictionary for the comparison of taxa
        topo_dict = dict()
        # before calculating RF distance check if both trees have the same taxa
        comp_dict = original_tree.compare(generated_tree)
        topo_dict["rf"] = comp_dict["rf"]
        topo_dict["max_rf"] = comp_dict["max_rf"]
        # modified normalized rf distance = 1 - rf / max_rf 
        topo_dict["rf_ratio"] = 1 - comp_dict["rf"] / comp_dict["max_rf"]
        # compare edges
        topo_dict["ref_edges_in_source"] = len(list(comp_dict["ref_edges_in_source"]))
        count_original_edges = len(list(comp_dict["ref_edges"]))
        count_common_edges = len(list(comp_dict["common_edges"]))
        count_missing_edges = count_original_edges - count_common_edges
        topo_dict["count_original_edges"] = count_original_edges
        topo_dict["count_missing_edges"] = count_missing_edges
        topo_dict["count_common_edges"] = count_common_edges
        topo_dict["correct_edges_ratio"] = round(count_common_edges/count_original_edges, 4)
        topo_dict["treeko_dist"] = comp_dict["treeko_dist"]
        # compare multifurcations
        topo_dict["count_multifurcations_original"] = ut.count_multifurcations(original)
        topo_dict["count_multifurcations_generated"] = ut.count_multifurcations(generated)
        return topo_dict
    
    def compare_distances(self):    
        """
        Compares distances of the original newick to the distances of the generated newick and returns dict with 
        following keys: \n
        mean_abs_leaf_dist_diff: mean of the absolute differences in leaf branch lengths of corresponding leaves,
        leaf branch lengths refer to the length of the branch that connects a leaf to its parent and the difference is 
        calculated for corresponding leaves i.e. branch length of taxon1 in original newick - branch length of taxon1 in 
        generated newick
        median_abs_leaf_dist_diff: median 
        mean_neg_leaf_dist_diff: mean of the negative differences of the above leaf branch lengths i.e. those branches 
        the model made longer than in the original newick 
        median_neg_leaf_dist_diff: median 
        mean_pos_leaf_dist_diff: mean of the positive differences i.e. those branches the model made shorter
        median_pos_leaf_dist_diff: median
        mean_abs_pairwise_leaf_dist_diff: mean of the differences of all pairwise leaf distances. For each pair of 
        leaves in the original tree the length of the path connecting them is computed. The same is done for the 
        corresponding pair in the generated tree and then substracted from the length of the original path.  
        median_abs_pairwise_leaf_dist_diff: median
        

        Returns:
            dict: dict with comparisons
        """
        dist_dict = dict()
        # get newicks
        original = self.original_newick
        generated = self.generated_newick
        # get taxa pairs
        taxa_pairs = get_taxon_pairs_greedy(ut.get_taxa(original),ut.get_taxa(self.generated_newick))
        for original_taxon, generated_taxon in taxa_pairs.items():
            generated.replace(generated_taxon, original_taxon)
        # create trees with updated taxa
        original_tree = Tree(original)
        generated_tree = Tree(generated)
        # save every leaf-parent distance difference
        leaf_parent_dist_diffs = []
        # get leaf nodes
        orig_leaf_nodes = [node for node in original_tree.traverse() if node.is_leaf()]
        gen_leaf_nodes = [node for node in generated_tree.traverse() if node.is_leaf()]
        # iterate over leaf nodes and save the differences
        for orig_leaf in orig_leaf_nodes:
            for gen_leaf in gen_leaf_nodes:
                if orig_leaf.name == gen_leaf.name:
                    leaf_parent_dist_diffs.append(orig_leaf.dist - gen_leaf.dist) 
                else:
                    # if the leaf doesnt exist in the generated tree then append the maximum difference
                    leaf_parent_dist_diffs.append(orig_leaf.dist) 
        # calculate mean and median over absolute differences
        abs_leaf_parent_dist_diffs = map(abs,leaf_parent_dist_diffs)
        dist_dict["mean_abs_diff"] = mean(abs_leaf_parent_dist_diffs)
        dist_dict["median_abs_diff"] = median(abs_leaf_parent_dist_diffs)
        # calculate mean diff over negative values i.e. generated edge is longer than original one
        neg_leaf_parent_dist_diffs = [diff for diff in leaf_parent_dist_diffs if diff < 0]
        # "on average how much shorter are shorter edges i.e. edges that the model made shorter than they actually are?"
        dist_dict["mean_neg_diff"] = mean(neg_leaf_parent_dist_diffs)
        dist_dict["median_neg_diff"] = median(neg_leaf_parent_dist_diffs)
        # calculate mean diff over positive values i.e. generated edge is shorter than original one
        pos_leaf_parent_dist_diffs = [diff for diff in leaf_parent_dist_diffs if diff > 0]
        dist_dict["mean_pos_diff"] = mean(pos_leaf_parent_dist_diffs)
        dist_dict["median_pos_diff"] = median(pos_leaf_parent_dist_diffs)
        # list for pairwise distances
        abs_pairwise_distances = []
        generated_taxa = ut.get_taxa(generated_tree)
        # get all leaf to leaf distances and calculate difference to corresponding distance in generated tree
        for orig_leaf1, orig_leaf2 in combinations(orig_leaf_nodes):
            ############# DEBUGGING
            console_logger.info(f"leaf 1: {orig_leaf1}, leaf 2: {orig_leaf2}")
            # get distance between both leaves in the original tree
            orig_dist = original_tree.get_distance(orig_leaf1, orig_leaf2)
            # get distance between both leaves in the generated tree
            if orig_leaf1 in generated_taxa and orig_leaf2 in generated_taxa:
                generated_dist = generated_tree.get_distance(orig_leaf1, orig_leaf2)
            else:
                # if the pair doesnt exist in the generated tree then let the difference be the maximum distance possible
                generated_dist = 0
            abs_pairwise_distance = abs(orig_dist-generated_dist)
            ####### DEBUGGING
            console_logger.info(f"distance: {abs_pairwise_distance}")
            abs_pairwise_distances.append(abs_pairwise_distance)
        dist_dict["mean_pairwise_diff"] = mean(abs_pairwise_distances)
        dist_dict["median_pairwise_diff"] = median(abs_pairwise_distances)
        return dist_dict

def main():
    argument_parser = argparse.ArgumentParser(
        description="Module for comparing two newicks e.g. an original newick of a dataset and an AI-generated newick.",
    )
    argument_parser.add_argument("-n", "--original_newick", required=True, type=str,
                                 help="Path to the .nwk of the original Newick.")
    argument_parser.add_argument("-g", "--generated_newick", required=True, type=str,
                                 help="Path to the .nwk of the AI-generated Newick or any newick that should be "
                                 "compared to the original one.")
    argument_parser.add_argument(
        "-o", 
        "--outfile", 
        required=False, 
        type=str,
        help="Filepath where the .tsv containing the comparisons will be saved at. If there already exists a file at "
        "the filepath then the .tsv entry will be added to the existing file. If no path is provided the result will "
        "be printed to console.")
    argument_parser.add_argument(
        "--quiet", 
        required=False, 
        action="store_true", 
        default=False, 
        help="On/Off flag. If --quiet is specified all console logs will be disabled and only the output printed out. "
        "This is useful inside a pipeline where the newick is piped into another application.")
    argument_parser.add_argument(
        "-p",
        "--params",
        required=False,
        help="Option for merging a given parameter .tsv with the comparison .tsv in order to have all values necessary "
        "for the complete performance analysis in one file. This appends the given .tsv's header to the comparison "
        "header and adds the first entry of the parameter .tsv to the comparison entry."
    )
    ############### ARGUMENTS ###############
    args = argument_parser.parse_args()
    path_original_newick = args.original_newick
    path_generated_newick = args.generated_newick
    outfile_path = args.outfile
    quiet = args.quiet
    params_path = args.params
    ########## CONFIGURE LOGGER ##########
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(logFormatter)
    console_logger.addHandler(stream_handler)
    if quiet:
        console_logger.setLevel(logging.WARNING)
    else:
        console_logger.setLevel(logging.INFO)
    console_logger.info("Starting newick_comparison")
    ########## CREATE COMPARISON_JOB OBJECT ##########
    # get newicks from their files
    original_newick = get_newick_from_file(path_original_newick)
    generated_newick = get_newick_from_file(path_generated_newick)
    # check if newicks are valid and set their format
    format_original = ut.get_newick_format(original_newick)
    format_generated = ut.get_newick_format(generated_newick)
    console_logger.info(f"Newick 1: {original_newick}")
    console_logger.info(f"Newick 2: {generated_newick}")
    # set flags, if one of the two newicks e.g. doesnt have distances, their distances wont be compared
    if format_original == 100 or format_generated == 100:
        topo_only = True
    else:
        topo_only = False
    if format_original == 9 or format_generated == 9:
        taxa_only = True
    else: 
        taxa_only = False
    # create job object
    comparison_job = Comparison_Job(
        original_newick_path=path_original_newick,
        generated_newick_path=path_generated_newick,
        original_newick=original_newick,
        generated_newick=generated_newick,
        format_original=format_original,
        format_generated=format_generated,
    )
    ########## CHECKS ##########
    # if both newicks have different formats warn user and only compare attributes both newicks have
    if not format_original == format_generated:
        console_logger.warning(f"Newicks have different formats. Format original newick: {format_original}. Format "
                         f"generated newick: {format_generated}.")
    console_logger.info(
        f"Comparing {"just topology" if topo_only else "just taxa" if taxa_only else "taxa and branch lengths"}."
    )
    if outfile_path:
        if os.path.isdir(outfile_path):
            raise IsADirectoryError("--outfile: Expected filepath but got directory path.")
    if params_path:
        if not os.path.exists(params_path):
            raise FileNotFoundError("--info: File does not exist")
        elif os.path.isdir(params_path):
            raise IsADirectoryError("--info: Expected filepath but got directory path.")
    ########## CREATE OUTPUT ##########
    # get the info of the passed tsv
    info_header = None
    info_entry = None
    if params_path:
        with open(params_path, "r") as tsv:
            info_header = tsv.readline()
            info_entry = tsv.readline()
    # calculate comparisons for taxa-only pairs
    if taxa_only:
        taxa_comp_dict = comparison_job.compare_taxa()
        dist_comp_dict = comparison_job.compare_distances()
        topo_comp_dict = comparison_job.compare_topology()
    # calculate only topology comparison
    elif topo_only:
        # TODO add comparison for topology-only trees
        taxa_comp_dict = None
        dist_comp_dict = None
        topo_comp_dict = None
    else:
        taxa_comp_dict = comparison_job.compare_taxa()
        dist_comp_dict = comparison_job.compare_distances()
        topo_comp_dict = comparison_job.compare_topology()
    # create tsv header and entry
    tsv_header = comparison_job.get_tsv_header(info_header)
    tsv_entry = comparison_job.get_tsv_entry(info_entry, taxa_comp=taxa_comp_dict,dist_comp=dist_comp_dict,topo_comp=topo_comp_dict)
    # either print the results or write them to file 
    if outfile_path:
        if not os.path.exists(outfile_path):
            console_logger.info("File doesn't exist. Creating .tsv with header.")
            with open(outfile_path, "w") as tsv_file:
                tsv_file.write(tsv_header)
        else: 
            console_logger.info("File exists. Writing .tsv entry directly.")
        with open(outfile_path, "a") as tsv_file:
            tsv_file.write(tsv_entry)
    else: 
        print(tsv_header)
        print(tsv_entry)
    console_logger.info("Finished newick_comparison")
# execute main method
if __name__ == "__main__":
    main()
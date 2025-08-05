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

# create logger
console_logger = logging.getLogger(__name__)

def get_newick_from_file(path):
    """
    Method for retrieving the Newick string from its file
    Args:
        path (str): path to newick

    Raises:
        ValueError: Newick is not valid
        
    Returns:
        str: newick 
    """
    if not os.path.exists(path) or not os.path.isfile(path):
        raise FileNotFoundError("Path to newick is not valid.")
    with open(path, "r") as nwk_file:
        return nwk_file.read()
    
    
def get_taxa(newick):
    """
    Returns all taxa from the given string in left to right order

    Args:
        newick (str): Newick string

    Returns:
        List(str): List of all taxa found 
    """
    return re.findall(r"(?<=[,\(])[A-Za-z0-9_.]+(?=:)", newick)

def get_filename(filepath):
    """
    Given a filepath returns the files name.
    Throws exceptions if path doesnt exist or isnt a filepath. 

    Args:
        filepath (_type_): path to a file

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
        taxa_only = False,
        topo_only = False, # flag for when newick is topology-only
    ):
        self.original_newick_path = original_newick_path
        self.generated_newick_path = generated_newick_path
        self.original_newick = original_newick
        self.generated_newick = generated_newick
        self.taxa_only = taxa_only
        self.topo_only = topo_only
        
    def get_tsv_header(self, info_header):
        """
        Function that returns the tsv header

        Returns:
            str: .tsv header
        """
        tsv_header = "newick1\tnewick2\ttaxa_only\ttopo_only\t#_taxa1\t#_taxa2\tcorrect_taxa_ratio\t#_equal_length\t"
        tsv_header += "#_unequal_length\tmean_ham_dist\tmean_ham_ratio\tmean_edit_dist\tmean_edit_dist_total\t"
        tsv_header += "mean_edit_ratio\tmean_edit_ratio_total\trf_dist\tmax_rf_dist\trf_ratio\t"
        tsv_header += "#_orig_edges\t#_gen_edges\t#_common_edges\t#_missing_edges\tcorrect_edge_ratio\t"
        if info_header:
            tsv_header += info_header
        tsv_header += "\n"
        return tsv_header

    def get_tsv_entry(self, param_entry, taxa_comp, dist_comp, topo_comp):
        
        console_logger.info(f"topo-only: {self.topo_only}")
        tsv_entry = f"{get_filename(self.original_newick_path)}\t{get_filename(self.generated_newick_path)}\t"
        tsv_entry+=f"{self.taxa_only}\t{self.topo_only}\t"
        if taxa_comp:
            tsv_entry+=f"{taxa_comp["taxa_count_original"]}\t{taxa_comp["taxa_count_generated"]}\t"
            tsv_entry+=f"{taxa_comp["correct_taxa_ratio"]}\t{taxa_comp["count_equal_length"]}\t{taxa_comp["count_unequal_length"]}\t"
            tsv_entry+=f"{taxa_comp["mean_hamming_distance"]}\t{taxa_comp["mean_hamming_distance_ratio"]}\t"
            tsv_entry+=f"{taxa_comp["mean_edit_distance"]}\t{taxa_comp["mean_edit_distance_total"]}\t{taxa_comp["mean_edit_ratio"]}\t"
            tsv_entry+=f"{taxa_comp["mean_edit_ratio_total"]}\t"
        else:
            tsv_entry+="None\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\t"
        if topo_comp:
            tsv_entry+=f"{topo_comp["rf"]}\t{topo_comp["max_rf"]}\t{topo_comp["rf_ratio"]}\t{topo_comp["count_original_edges"]}\t"
            tsv_entry+=f"{topo_comp["count_generated_edges"]}\t{topo_comp["count_common_edges"]}\t"
            tsv_entry+=f"{topo_comp["count_missing_original_edges"]}\t{topo_comp["correct_edges_ratio"]}\t"
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
        original_taxa = get_taxa(self.original_newick)
        generated_taxa = get_taxa(self.generated_newick)
        # count taxa in both newicks
        taxa_dict["taxa_count_original"] = len(original_taxa)
        taxa_dict["taxa_count_generated"] = len(generated_taxa)
        # percentage of correct taxa, average hamming distance, average levenshtein distance
        match_counter = 0
        taxon_pairs = get_taxon_pairs_greedy(original_taxa, generated_taxa)
        for original, generated in taxon_pairs.items():
            if original == generated:
                match_counter += 1
        taxa_dict["correct_taxa_ratio"] = round(match_counter/len(original_taxa), 4) 
        # calculate average hamming distance and ratio for all taxa pairs whose strings have equal length
        accu_hamming_distances = 0
        accu_hamming_ratios = 0
        # count number of taxa with equal length (implies substitutions), unequal length (implies insertions/deletions)
        equal_length_counter = 0
        unequal_length_counter = 0
        # save which original taxa dont have a matching generated taxon
        mismatches = []
        for original, generated in taxon_pairs.items():
            if generated: 
                if len(original) == len(generated):
                    hdist = hamming_distance(original, generated)
                    accu_hamming_distances += hdist
                    accu_hamming_ratios += 1 - (hdist/len(original))
                    equal_length_counter += 1
                else:
                    unequal_length_counter += 1
            else:
                mismatches.append(original)
        taxa_dict["count_equal_length"] = equal_length_counter
        taxa_dict["count_unequal_length"] = unequal_length_counter
        taxa_dict["mismatches"] = mismatches
        taxa_dict["mean_hamming_distance"] = round(accu_hamming_distances/equal_length_counter, 4)
        taxa_dict["mean_hamming_distance_ratio"] = round(accu_hamming_ratios/equal_length_counter, 4)
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
        To compare the topologies we need to change the taxa of the generated newick to their counterpart in the original newick
        i.e. the actual taxa the AI just read wrong (correcting spelling mistakes made by the AI). 

        Args:
            original (str): original newick
            generated (str): generated newick
        """
        original = self.original_newick
        generated = self.generated_newick
        taxa_pairs = get_taxon_pairs_greedy(get_taxa(original),get_taxa(generated))
        # correct spelling mistakes by AI
        for original_taxon, generated_taxon in taxa_pairs.items():
            generated.replace(generated_taxon, original_taxon)
        original_tree = Tree(original, format=1)
        generated_tree = Tree(generated, format=1)
        # TODO can midpoint outgrouping change the results? Why do i have to midpoint outgroup?
        original_tree.set_outgroup(original_tree.get_midpoint_outgroup())
        generated_tree.set_outgroup(generated_tree.get_midpoint_outgroup())
        # dictionary for the comparison of taxa
        topo_dict = dict()
        # before calculating RF distance check if both trees have the same taxa
        rf, max_rf, _, original_edges, generated_edges, _, _ = original_tree.robinson_foulds(generated_tree)
        topo_dict["rf"] = rf
        topo_dict["max_rf"] = max_rf
        # modified normalized rf distance = 1 - rf / max_rf 
        topo_dict["rf_ratio"] = 1 - rf / max_rf
        # compare edges
        topo_dict["count_original_edges"] = len(original_edges)
        topo_dict["count_generated_edges"] = len(generated_edges)
        edge_intersection = list(set(original_edges).intersection(set(generated_edges)))
        topo_dict["count_common_edges"] = len(edge_intersection)
        topo_dict["common_edges"] = edge_intersection
        missing_original_edges = list(set(original_edges).difference(set(generated_edges)))
        topo_dict["missing_original_edges"] = missing_original_edges
        topo_dict["count_missing_original_edges"] = len(missing_original_edges)
        topo_dict["correct_edges_ratio"] = round(len(edge_intersection)/len(original_edges), 4)
        return topo_dict
    
    # TODO: for those common edges how different from each other are the distances, 
    def compare_distances(self):    
        distances_dict = dict()
    


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
    console_logger.info("Starting")
    ########## CREATE COMPARISON_JOB OBJECT ##########
    # get newicks from their files
    original_newick = get_newick_from_file(path_original_newick)
    generated_newick = get_newick_from_file(path_generated_newick)
    # check if newicks are valid and set their format
    format_original = ut.get_newick_format(original_newick)
    format_generated = ut.get_newick_format(generated_newick)
    # set flags, if one of the two newicks e.g. doesnt have distances, their distances wont be compared
    if format_original == 100 or format_generated == 100:
        topo_only = True
    else:
        topo_only = False
    if format_original == 9 or format_generated == 9:
        taxa_only = True
    else: 
        taxa_only = False
    console_logger.info(
        f"Newicks {"are topology-only" if topo_only else "are taxa-only" if taxa_only else "have taxa and branch lengths"}."
    )
    # create job object
    comparison_job = Comparison_Job(
        original_newick_path=path_original_newick,
        generated_newick_path=path_generated_newick,
        original_newick=original_newick,
        generated_newick=generated_newick,
        taxa_only=taxa_only,
        topo_only=topo_only,
    )
    ########## CHECKS ##########
    # if both newicks have different formats warn user and only compare attributes both newicks have
    if not format_original == format_generated:
        raise ValueError(f"Newicks have different formats. Format original newick: {format_original}. Format "
                         f"generated newick: {format_generated}.")
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
    if comparison_job.taxa_only:
        taxa_comp_dict = comparison_job.compare_taxa()
        dist_comp_dict = None
        topo_comp_dict = comparison_job.compare_topology()
    elif comparison_job.topo_only:
        # TODO add comparison for topology-only trees
        taxa_comp_dict = None
        dist_comp_dict = None
        topo_comp_dict = None
    else:
        # TODO add comparison of distances
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
    console_logger.info("Finished")
# execute main method
if __name__ == "__main__":
    main()
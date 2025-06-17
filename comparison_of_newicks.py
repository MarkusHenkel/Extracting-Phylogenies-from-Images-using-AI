from data_collection import phylogenetic_trees_and_newicks as ptan # writing the phylogenetic tree of the generated newick
import argparse # argument parsing 
from argparse import RawTextHelpFormatter
import datetime # for error messages and logging
import re # parsing newick string
import nltk # Levenshtein
from ete3 import Tree # turn newick into tree structure
import os # for checking outfile paths of the tsv

# Method for retrieving the Newick string from file
def get_newick(path):
    if not path:
        exit(f"[{datetime.datetime.now()}] Error in get_newick(): Path cannot be None.")
    with open(path, "r") as nwk_file:
        newick = nwk_file.read()
        return newick
    
def get_taxa(newick):
    """
    Returns all taxa from the given string in left to right order

    Args:
        newick (str): Newick string

    Returns:
        List(str): List of all taxa found 
    """
    return re.findall(r"(?<=[,\(])[A-Za-z0-9_.]+(?=:)", newick)

def get_filename(path):
    """
    Given a filepath returns the name of the file the path points to. 
    Returns an error if the path isnt valid or if the path points to a directory.

    Args:
        path (str): path to the file

    Returns:
        str: filename
    """
    filename_match = re.search(r"[\w\-\ \.]+(?=\.)", path)
    if filename_match:
        return filename_match.group()
    else:
        exit(f"[{datetime.datetime.now()}] Error in get_filename(): Path doesn't contain a file name.")
    
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
        exit(f"[{datetime.datetime.now()}] Error in hamming_distance(): Strings have to be of equal length.")
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
    if len(original_taxa) < len(generated_taxa):
        print(f"[{datetime.datetime.now()}] Warning in get_taxon_pairs_greedily: More taxa in generated Newick than in original one.")
    elif len(original_taxa) > len(generated_taxa):
        print(f"[{datetime.datetime.now()}] Warning in get_taxon_pairs_greedily: More taxa in original Newick than in generated one.")
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

# TODO: count insertions, deletions, substitutions
def compare_taxa(original, generated):
    """
    Given the strings of two newicks returns a dictionary of statistics of the taxa of both

    Args:
        original (str): original newick that the generated one is compared to
        generated (str): (AI-generated) newick whose quality is calculated

    Returns:
        dict: statistics of how similar the taxa of both newicks are
    """
    # check if newicks are valid
    if not ptan.is_newick(original, 1): 
        exit(f"[{datetime.datetime.now()}] Error in compare_taxa(): Newick 1 is not formatted correctly.")
    elif not ptan.is_newick(generated, 1):
        exit(f"[{datetime.datetime.now()}] Error in compare_taxa(): Newick 2 is not formatted correctly.")
    elif not get_taxa(original) or not get_taxa(generated):
        exit(f"[{datetime.datetime.now()}] Error in compare_taxa(): Lists of taxa can't be empty.")
    comparison_dict = dict()
    original_taxa = get_taxa(original)
    generated_taxa = get_taxa(generated)
    # count taxa in both newicks
    comparison_dict["taxa_count_original"] = len(original_taxa)
    comparison_dict["taxa_count_generated"] = len(generated_taxa)
    #
    # percentage of correct taxa, average hamming distance, average levenshtein distance
    # 
    match_counter = 0
    taxon_pairs = get_taxon_pairs_greedy(original_taxa, generated_taxa)
    for original, generated in taxon_pairs.items():
        if original == generated:
            match_counter += 1
    comparison_dict["correct_taxa_ratio"] = round(match_counter/len(original_taxa), 4) 
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
    comparison_dict["count_equal_length"] = equal_length_counter
    comparison_dict["count_unequal_length"] = unequal_length_counter
    if mismatches:
        comparison_dict["mismatches"] = mismatches
    comparison_dict["mean_hamming_distance"] = round(accu_hamming_distances/equal_length_counter, 4)
    comparison_dict["mean_hamming_distance_ratio"] = round(accu_hamming_ratios/equal_length_counter, 4)
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
    comparison_dict["mean_edit_distance"] = round(accu_edit_distances/len(taxon_pairs), 4)
    # mean distance calculated over all pairs (even those where the generated taxon is "" => max distance)
    # => should generally be lower than mean_edit_distance
    comparison_dict["mean_edit_distance_total"] = round(accu_edit_distances/(len(taxon_pairs)- len(mismatches)), 4)
    # also calculate a ratio over pairs with generated newick and over all pairs
    comparison_dict["mean_edit_ratio"] = round(accu_edit_ratios/(len(taxon_pairs)- len(mismatches)), 4)
    comparison_dict["mean_edit_ratio_total"] = round(accu_edit_ratios/len(taxon_pairs), 4)
    return comparison_dict

# TODO: count the number of edges in both trees, count the number of edges shared by both, ratio, mean of the average
def compare_topology(original, generated):
    """
    Given 2 Newick strings compares the topologies of their resulting trees.
    To compare the topologies we need to change the taxa of the generated newick to their counterpart in the original newick
    i.e. the actual taxa the AI just read wrong (correcting spelling mistakes made by the AI). 

    Args:
        original (str): original newick
        generated (str): generated newick
    """
    if not ptan.is_newick(original, 1): 
        exit(f"[{datetime.datetime.now()}] Error in compare_topology(): Original Newick is not formatted correctly.")
    elif not ptan.is_newick(generated, 1):
        exit(f"[{datetime.datetime.now()}] Error in compare_topology(): Generated Newick is not formatted correctly.")
    taxa_pairs = get_taxon_pairs_greedy(get_taxa(original),get_taxa(generated))
    # correct spelling mistakes by AI
    for original_taxon, generated_taxon in taxa_pairs.items():
        generated.replace(generated_taxon, original_taxon)
    original_tree = Tree(original, format=1)
    generated_tree = Tree(generated, format=1)
    original_tree.set_outgroup(original_tree.get_midpoint_outgroup())
    generated_tree.set_outgroup(generated_tree.get_midpoint_outgroup())
    comparison_dict = dict()
    # before calculating RF distance check if both trees have the same taxa
    rf, max_rf, common_leaves, original_edges, generated_edges, original_discarded, generated_discarded = original_tree.robinson_foulds(generated_tree)
    comparison_dict["rf"] = rf
    comparison_dict["max_rf"] = max_rf
    # modified normalized rf distance = 1 - rf / max_rf 
    comparison_dict["rf_ratio"] = 1 - rf / max_rf
    # compare edges
    comparison_dict["count_original_edges"] = len(original_edges)
    comparison_dict["count_generated_edges"] = len(generated_edges)
    edge_intersection = list(set(original_edges).intersection(set(generated_edges)))
    comparison_dict["count_common_edges"] = len(edge_intersection)
    comparison_dict["common_edges"] = edge_intersection
    missing_original_edges = list(set(original_edges).difference(set(generated_edges)))
    comparison_dict["missing_original_edges"] = missing_original_edges
    comparison_dict["count_missing_original_edges"] = len(missing_original_edges)
    comparison_dict["correct_edges_ratio"] = round(len(edge_intersection)/len(original_edges), 4)
    return comparison_dict
    
    
# TODO: for those common edges how different from each other are the distances, 
# def compare_distances(original, generated):    

# TODO: implement tsv writer class: 
# TODO: write the tsv header
def write_tsv_header(outfile_path):
    if os.path.exists(outfile_path):
        exit(f"[{datetime.datetime.now()}] Error in write_tsv_header: File already exists.")
    tsv_header = "newick1\tnewick2\t#_taxa1\t#_taxa2\tcorrect_taxa_ratio\t#_equal_length\t"
    tsv_header += "#_unequal_length\tmean_ham_dist\tmean_ham_ratio\tmean_edit_dist\tmean_edit_dist_total\t"
    tsv_header += "mean_edit_ratio\tmean_edit_ratio_total\trf_dist\tmax_rf_dist\trf_ratio\t"
    tsv_header += "#_orig_edges\t#_gen_edges\t#_common_edges\t#_missing_edges\tcorrect_edge_ratio"
    tsv_header += "\n"
    with open(outfile_path, "w") as tsv_file:
        tsv_file.write(tsv_header)
    print(f"[{datetime.datetime.now()}] write_tsv_header: File with header was created.")

def write_tsv_entry(outfile_path, path_original_newick, path_generated_newick):
    if not os.path.exists(outfile_path):
        exit(f"""[{datetime.datetime.now()}] Error in write_tsv_entry: 
            File hasn't been initialized.""")
    original_newick = get_newick(path_original_newick)
    generated_newick = get_newick(path_generated_newick)
    original_newick_name = get_filename(path_original_newick)
    generated_newick_name = get_filename(path_generated_newick)
    taxa_comp = compare_taxa(original_newick, generated_newick)
    top_comp = compare_topology(original_newick, generated_newick)
    tsv_entry = f"{original_newick_name}\t{generated_newick_name}\t"
    tsv_entry += f"{taxa_comp["taxa_count_original"]}\t{taxa_comp["taxa_count_generated"]}\t"
    tsv_entry += f"{taxa_comp["correct_taxa_ratio"]}\t{taxa_comp["count_equal_length"]}\t{taxa_comp["count_unequal_length"]}\t"
    tsv_entry += f"{taxa_comp["mean_hamming_distance"]}\t{taxa_comp["mean_hamming_distance_ratio"]}\t"
    tsv_entry += f"{taxa_comp["mean_edit_distance"]}\t{taxa_comp["mean_edit_distance_total"]}\t{taxa_comp["mean_edit_ratio"]}\t"
    tsv_entry += f"{taxa_comp["mean_edit_ratio_total"]}\t{top_comp["rf"]}\t{top_comp["max_rf"]}\t{top_comp["rf_ratio"]}\t"
    tsv_entry += f"{top_comp["count_original_edges"]}\t{top_comp["count_generated_edges"]}\t{top_comp["count_common_edges"]}\t"
    tsv_entry += f"{top_comp["count_missing_original_edges"]}\t{top_comp["correct_edges_ratio"]}"
    with open(outfile_path, "a") as tsv_file:
        tsv_file.write(tsv_entry)
    print(f"[{datetime.datetime.now()}] write_tsv_entry: Entry was written to file.")
        
def main():
    argument_parser = argparse.ArgumentParser(
        description="""Comparison of an AI-generated Newick and the original Newick.
        Modes:\n
        \tgenerate_image: Use the phylogenetic_trees_and_newicks script to generate an image of a given Newick.
        \tcompare_newicks: Compares two given Newicks by their taxa, distances and general topology. Results are saved to .tsv.""",
        formatter_class=RawTextHelpFormatter
    )
    argument_parser.add_argument("--mode", required=True, choices=["generate_image", "compare_newicks"], 
                                 help="Specify which functionality you want to use.")
    argument_parser.add_argument("-n", "--original_newick", required=False, type=str,
                                 help="Path to the .nwk of the original Newick.")
    argument_parser.add_argument("-g", "--generated_newick", required=False, type=str,
                                 help="Path to the .nwk of the AI-generated Newick.")
    argument_parser.add_argument("-o", "--outfile", required=False, type=str,
                                 help="Path where output is created at.")
    args = argument_parser.parse_args()
    path_original_newick = args.original_newick
    path_generated_newick = args.generated_newick
    outfile_path = args.outfile
    mode = args.mode
    if mode == "generate_image":
        if not outfile_path or (not path_generated_newick and not path_original_newick):
            exit(f"""[{datetime.datetime.now()}] Error in generate_image mode: When using '--mode generate_image' please
                provide the path where the image is supposed to be saved at and either an original newick or a generated newick.""")
        elif path_generated_newick:
            newick = get_newick(path_generated_newick)
            # Generate the image of the phylogenetic tree that represents the AI generated Newick
            ptan.save_newick_image(newick=newick, outfile_path=outfile_path)
            print(f"[{datetime.datetime.now()}] generate_image mode: Image was saved.")
        elif path_original_newick:
            newick = get_newick(path_original_newick)
            # Generate the image of the phylogenetic tree that represents the AI generated Newick
            ptan.save_newick_image(newick=newick, outfile_path=outfile_path)
            print(f"[{datetime.datetime.now()}] generate_image mode: Image was saved.")
    elif mode == "compare_newicks":
        if not outfile_path or not path_generated_newick or not path_original_newick:
            exit(f"[{datetime.datetime.now()}] Error in compare_newicks mode: When using '--mode compare_newicks' "
                 "please provide two newicks and the path where the .tsv containing the comparison is supposed to be saved at.")
        write_tsv_header(outfile_path)
        write_tsv_entry(outfile_path, path_original_newick, path_generated_newick)
        
# execute main method
main()
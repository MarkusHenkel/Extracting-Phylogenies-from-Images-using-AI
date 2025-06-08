from data_collection import phylogenetic_trees_and_newicks
import argparse
import datetime
import re
import nltk


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

def get_taxon_pairs_greedy(original_taxa, generated_taxa):
    """
    Given a dictionary with two lists: mismatched original taxa and mismatched generated taxa,
    returns a dict with each key being an original taxon and the value being the corresponding generated taxon.
    Iterates through the original taxa and greedily assigns the best generated taxon to the current original taxon.

    Args:
        original_taxa (List(str)): List of original taxa that couldnt be matched
        generated_taxa (List(str)): List of generated taxa that couldnt be matched
    Returns:
        dict: dict with original_taxon/generated_taxa pairs
    """
    pair_dict = dict()
    for original_taxon in original_taxa:
        score_dict = dict()
        for generated_taxon in generated_taxa:
            score_dict[generated_taxon] = nltk.edit_distance(
                original_taxon, generated_taxon, substitution_cost=1, transpositions=False
            )
        pair_dict[original_taxon] = min(score_dict, key=score_dict.get)
        # remove the generated taxon that was just assigned from the list of generated taxa
        generated_taxa.remove(pair_dict[original_taxon])
    return pair_dict 

def compare_taxa(original, generated):
    """
    Given the strings of two newicks returns a dictionary of statistics of the taxa of both

    Args:
        original (str): original newick that the generated one is compared to
        generated (str): (AI-generated) newick whose quality is calculated

    Returns:
        dict: statistics of how similar the taxa of both newicks are
    """
    if not get_taxa(original) or not get_taxa(generated):
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
    # calculate number of pairs with strings of equal length
    equal_length_counter = 0
    for original, generated in taxon_pairs.items():
        if len(original) == len(generated):
            hdist = hamming_distance(original, generated)
            accu_hamming_distances += hdist
            accu_hamming_ratios += 1 - (hdist/len(original))
            equal_length_counter += 1
    comparison_dict["mean_hamming_distance"] = round(accu_hamming_distances/equal_length_counter, 4)
    comparison_dict["mean_hamming_distance_ratio"] = round(accu_hamming_ratios/equal_length_counter, 4)
    # calculate average edit distance and ratio for all taxa pairs
    accu_edit_distances = 0
    accu_edit_ratios = 0
    for original, generated in taxon_pairs.items():
        edit_distance = nltk.edit_distance(original, generated, substitution_cost=1, transpositions=False)
        accu_edit_distances += edit_distance
        accu_edit_ratios += 1 - (edit_distance/max(len(original), len(generated)))
    comparison_dict["mean_edit_distance"] = round(accu_edit_distances/len(taxon_pairs), 4)
    comparison_dict["mean_edit_ratio"] = round(accu_edit_ratios/len(taxon_pairs), 4)
    # 
    # Testing
    # 
    # print("---Mismatched original taxa---")
    # for taxon in mismatched_taxa["original_taxa"]:
    #     print(taxon)
    # print("---Mismatched generated taxa---")
    # for taxon in mismatched_taxa["generated_taxa"]:
    #     print(taxon)
        
    return comparison_dict

# def compare_distances(original, generated):

# def compare_topology(original, generated):
    
# TODO: method that extracts all subtrees
# TODO: implement tsv writer class: 
# TODO: implement comparison class?

def main():
    argument_parser = argparse.ArgumentParser(
        description="""Comparison of an AI-generated Newick and the original Newick.
        Modes:\n
        \tgenerate_image: Use the phylogenetic_trees_and_newicks script to generate an image of a given Newick.
        \tcompare_newicks: Compares two given Newicks by their taxa, distances and general topology."""
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
    path_outfile = args.outfile
    mode = args.mode
    original_newick = get_newick(path_original_newick)
    generated_newick = get_newick(path_generated_newick)
    if mode == "generate_image":
        if not path_outfile or (not path_generated_newick and not path_original_newick):
            exit(f"""[{datetime.datetime.now}] Error in generate_image mode: When using '--mode generate_image' please provide the path where the image 
                 is supposed to be saved at and either an original newick or a generated newick.""")
        elif path_generated_newick:
            newick = get_newick(path_generated_newick)
            # Generate the image of the phylogenetic tree that represents the AI generated Newick
            phylogenetic_trees_and_newicks.save_newick_image(newick=newick, outfile_path=path_outfile)
            print(f"[{datetime.datetime.now()}] generate_image mode: Image was saved.")
        elif path_original_newick:
            newick = get_newick(path_original_newick)
            # Generate the image of the phylogenetic tree that represents the AI generated Newick
            phylogenetic_trees_and_newicks.save_newick_image(newick=newick, outfile_path=path_outfile)
            print(f"[{datetime.datetime.now()}] generate_image mode: Image was saved.")
    elif mode == "compare_newicks":
        comparison_dict = compare_taxa(original_newick, generated_newick)
        print(comparison_dict)
# execute main method
main()
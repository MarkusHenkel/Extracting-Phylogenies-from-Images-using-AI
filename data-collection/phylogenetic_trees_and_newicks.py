from Bio import Phylo
from ete3 import NCBITaxa
from io import StringIO 
from numpy import random
import datetime
import argparse
import matplotlib.pyplot as plt
import re

# TODO: replace comments with pydocs

# ncbi taxonomy
ncbi = NCBITaxa()

# Check if a given taxid is valid 
def is_valid(taxid):
    return True if ncbi.get_taxid_translator([taxid]) else False 
    
# Counts amount of invalid taxids in the first 10000
# Function for testing
def amount_invalid():
    counter = 0
    for taxid in range(10000):
        if not is_valid(taxid):
            counter += 1
    return counter

# Given a preferred number of taxa (amount as int) and a boolean generates either a list with fixed or randomized amount [2, amount] of random valid taxa
#
# default amount of taxa is 10
# not all numbers from 1-1000 are valid entries => check if taxon ID is valid (dictionary of translator is not empty)
def generate_random_taxids(amount = 10, randomize = False):
    valid_rand_taxids = []
    # if randomize = True, randomize the amount of taxa (0-amount)
    if randomize:
        amount = random.randint(amount)
    # in any case let the minimum amount of taxa be 2
    if (amount <= 1):
        amount = 2
    # Add a valid taxid <amount> times
    i = 0
    while i < amount:
        rand_taxid = random.randint(10000)
        if is_valid(rand_taxid):
            valid_rand_taxids.append(rand_taxid)
            i += 1
    return valid_rand_taxids

# Given a list of random valid taxids returns a list of the corresponding taxa names 
# 
# output of get_taxid_translator is a dictionary where the keys are the taxon IDs and the values the taxa
def taxids_to_taxa(valid_rand_taxids):  
    valid_rand_taxa_dict = ncbi.get_taxid_translator(valid_rand_taxids)
    valid_rand_taxa = []
    for taxon in valid_rand_taxa_dict.values():
        valid_rand_taxa.append(taxon)
    return valid_rand_taxa

# Given a list of random valid taxa generates its topology and returns its newick tree as a string
def generate_newick_tree(valid_rand_taxa):
    return ncbi.get_topology(valid_rand_taxa).write()

# Given a newick tree as a string, generates an image of the newick tree 
#
# draws the tree using bio.phylo and matplotlib therefore must first be made into a newick file
def tree_to_image_file(newick):
    # make the the newick tree (string) into a file so that it can be drawn by bio.phylo (not elegant?)
    newick_file = StringIO(newick)
    newick_tree = Phylo.read(newick_file, "newick")
    # create a matplotlib figure and 
    newick_tree.rooted = True
    fig = plt.figure(figsize=(10, 20), dpi=150)
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(newick_tree, axes=axes, do_show=False)
    plt.savefig(f"tree_{str(datetime.datetime.now().strftime(r"%Y-%m-%d-%H-%S-%f"))}.jpg")
    plt.close()

# Given a newick tree return a list of all taxids inside the string using regex
# CAUTION: This regex only works with newicks with distances of 1 so before distances are added
def taxids_from_newick(newick_string):
    return re.findall(r"[^)(,:][\d]+[^:]", newick_string)
    
def translate_newick_tree(newick_string):
    """
    Given a newick tree with taxids returns the newick tree with the corresponding taxa

    Args:
        newick_string (str): The tree in newick format

    Returns:
        newick_string (str): The same tree but taxon IDs are replaced with the corresponding taxa
    """
    taxids = taxids_from_newick(newick_string)
    taxa = taxids_to_taxa(taxids)
    for i in range(len(taxids)):
        newick_string = newick_string.replace(taxids[i], taxa[i])
    return newick_string
    
# TODO: add function that edits newick file so that the distances of the newick tree are randomized before its drawn
# TODO: add function that puts the image and a text file with the newick inside a folder and 

def main():
    argument_parser = argparse.ArgumentParser(
        description="Data Collection for Extracting phylogenies from images using AI"
    )
    #
    # Arguments
    #
    
    # TODO: add parameter to make random distances optional (0 for standard distance of 1, 1 for randomized distances)
    # TODO: change name of image file to something that increments
    
    # parameter for the preferred number of taxa generated, if not specified defaults to 10 taxa
    argument_parser.add_argument('-a', '--amount_taxa', required=False, type=int,
                                 help='Type=Int. Choose the preferred amount of generated taxa. If not specified defaults to 10 taxa.')
    argument_parser.add_argument('-ra', '--randomize_amount', required=False, action='store_true',
                                 help='On/Off flag. Randomizes the amount of taxa generated. Range: 2 to <amount_taxa>. Per default the amount of taxa is not randomized.')
    # argument_parser.add_argument('-rd', '--randomize_distances', required=False, type=int,
    #                              help='On/Off flag. ')
    
    # Specified parameters
    args = argument_parser.parse_args()
    
    #
    # Main functionality
    #
    
    # if the amount of taxa is not specified it defaults to 10 taxa
    amount_taxa = args.amount_taxa
    randomize_amount = args.randomize_amount
    
    # if amount of taxa and randomize_amount are specified then call generate_random_taxids() with the amount of taxa and randomize set to True
    random_taxids = []
    # -a and -r specified
    if amount_taxa != None and randomize_amount == True:
        random_taxids = generate_random_taxids(amount=amount_taxa, randomize=True)
        newick_taxids = generate_newick_tree(random_taxids)
        newick_taxa = translate_newick_tree(newick_taxids)
    # -a specified, -r=False
    elif amount_taxa != None and randomize_amount == False:
        random_taxids = generate_random_taxids(amount=amount_taxa)
        newick_taxids = generate_newick_tree(random_taxids)
        newick_taxa = translate_newick_tree(newick_taxids)
    # -r specified, -a=None
    elif amount_taxa == None and randomize_amount == True:
        random_taxids = generate_random_taxids(randomize=True)
        newick_taxids = generate_newick_tree(random_taxids)
        newick_taxa = translate_newick_tree(newick_taxids)
    # neither specified, -a=None, -r=False
    else:
        random_taxids = generate_random_taxids()
        newick_taxids = generate_newick_tree(random_taxids)
        newick_taxa = translate_newick_tree(newick_taxids)
    # Generates the newick tree
    print(f"Randomize amount: {randomize_amount}. Specified amount of taxa: {amount_taxa}. Actual amount of taxa: {len(random_taxids)}.")
    print("Newick tree:", newick_taxa)
    tree_to_image_file(newick_taxa)
    
    
    
    #
    # Testing main functionality
    #
    # print()
    # print("With taxids:","(((((((8610:1,8834:1)1:1,9881:1)1:1,8343:1)1:1,7755:1)1:1,6957:1)1:1,(4466:1,4151:1)1:1)1:1,(1381:1,2366:1)1:1);")
    # print("Translation:", translate_newick_tree("(((((((8610:1,8834:1)1:1,9881:1)1:1,8343:1)1:1,7755:1)1:1,6957:1)1:1,(4466:1,4151:1)1:1)1:1,(1381:1,2366:1)1:1);"))
    # taxid_list = generate_random_taxids(amount=10, randomize=False)
    # newick = generate_newick_tree(taxid_list)
    # print(len(taxid_list))
    
    # print("Amount taxa:", amount_taxa)
    # print("Randomize amount:", randomize_amount)
    
    # there are 1149 invalid taxids in the first 10000 taxids
    # print(amount_invalid())
    
    
# execute the main method
main()

#
# Testing the ETE Toolkit functionality:
#

# 2890 is no a valid taxid, so the returned dictionary of translations is empty and it doesnt result in an error 
# => if dictionary is empty then the ID is not valid
# taxid_list = [2890]
# taxid_name = ncbi.get_taxid_translator(taxid_list)
# print(taxid_name)

# when translating the taxids into their respective names get_taxid_translator() will return a dictionary where the keys are the taxids 
# and the values are the taxon names
# test_taxids = [23, 24, 25]
# print(ncbi.get_taxid_translator(test_taxids))

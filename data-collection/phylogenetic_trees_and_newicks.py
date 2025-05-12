from Bio import Phylo
from ete3 import NCBITaxa
from io import StringIO 
from numpy import random
import argparse
import matplotlib.pyplot as plt

# ncbi taxonomy
ncbi = NCBITaxa()

# Given a preferred number of taxa (amount as int) generates a list of a random amount (between 2 and the chosen amount+2) of random valid taxa
# + 2 to ensure that theres always at least 2 taxon and the get_topology function doesnt return an error
# not all numbers from 1-1000 are valid entries => check if taxon ID is valid
def generate_random_taxids(amount):
    valid_rand_taxids = []
    for i in range(random.randint(amount) + 2):
        rand_taxid = [random.randint(10000)]
        if ncbi.get_taxid_translator(rand_taxid):
            valid_rand_taxids.append(rand_taxid[0])
    return valid_rand_taxids

# given a list of random valid taxids returns a list of the corresponding taxa names 
# output of get_taxid_translator is a dictionary where the keys are the taxon IDs and the values the 
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
# draws the tree using bio.phylo and matplotlib therefore must first be made into a newick file
def tree_to_image_file(newick):
    # make the the newick tree (string) into a file so that it can be drawn by bio.phylo (not elegant?)
    newick_file = StringIO(newick)
    newick_tree = Phylo.read(newick_file, "newick")
    newick_tree.rooted = True
    Phylo.draw(newick_tree)
    plt.savefig("phylo_tree.jpg")

# TODO: add function that edits newick file so that the generated images contain the taxa names not their IDs
# TODO: add function that edits newick file so that the distances of the newick tree are randomized before its drawn

def main():
    argument_parser = argparse.ArgumentParser(
        description="Data Collection for Extracting phylogenies from images using AI"
    )
    
    # arguments
    
    # TODO: add parameter for choosing if the amount of taxa is randomized or not
    
    # TODO: add parameter to make random distances optional (0 for standard distance of 1, 1 for randomized distances)
    
    # parameter for the preferred number of taxa generated, if not specified defaults to 10 taxa
    argument_parser.add_argument('-m', '--max_amount_taxa', required=False,
                                 help='The amount of taxa generated is randomized. Choose the preferred upper limit of generated taxa. If not specified defaults to a maximum 10 taxa.')
    # Specified parameters
    args = argument_parser.parse_args()
    
    #
    # Main functionality:
    #
    # if the amount of taxa is not specified it defaults to 10 taxa
    max_amount_taxa = args.max_amount_taxa
    if(max_amount_taxa == None):
        max_amount_taxa = 10
        
    # Generates the newick tree
    newick = generate_newick_tree(generate_random_taxids(max_amount_taxa))
    print(newick)
    tree_to_image_file(newick)
    
    
    
    
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

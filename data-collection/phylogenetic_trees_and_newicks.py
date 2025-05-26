from Bio import Phylo
from ete3 import NCBITaxa
from io import StringIO 
from numpy import random
import datetime
import argparse
import matplotlib.pyplot as plt
import re
import os
import pathlib


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

def tree_to_image_file(newick, path, file_id):
    """
    Given a newick tree as a string, generates an image of the newick tree and saves it to the specified path.
    Draws the tree using bio.phylo and matplotlib therefore must first be made into a newick file.

    Args:
        newick (str): newick string of the phylogenetic tree to be created
        path (str): path where the file is saved  
    """
    # make the the newick tree (string) into a file so that it can be drawn by bio.phylo (not elegant?)
    newick_file = StringIO(newick)
    newick_tree = Phylo.read(newick_file, "newick")
    # create a matplotlib figure and 
    fig = plt.figure(figsize=(10, 20), dpi=150)
    axes = fig.add_subplot(1, 1, 1)
    newick_tree.rooted = True
    Phylo.draw(newick_tree, axes=axes, do_show=False)
    image_path = path + f"\\tree_{file_id}.jpg" 
    # print("Test path: " + image_path)
    plt.savefig(image_path)
    plt.close()

def taxids_from_newick(newick_string):
    """
    Given a newick tree containing taxon IDs return a list of all taxon IDs inside the string using regex.
    CAUTION: This regex only works with newicks with distances of 1 so before distances are added.

    Args:
        newick_string (str): The tree in newick format

    Returns:
        List[str]: List containing all taxon IDs in order
    """
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
    # print("Taxa Test:")
    # print(taxa)
    for i in range(len(taxids)):
        newick_string = newick_string.replace(taxids[i], taxa[i].replace(" ","_"))
    return newick_string

def randomize_distances(newick_string, max_distance):
    """
    Given a newick tree and a max distance returns the newick tree with randomized distances. Also works on 
    Randomized distances are rounded to 2 decimal and are in the range from 0.00 to the specified max amount

    Args:
        newick_string (str): The tree in newick format
        max_distance (float): upper limit of randomized distances
    """
    # if re.sub() is given a lambda it tries to give the match object to the lambda but the lambda doesnt take any inputs
    # and instead generates a new random float each time a match is found
    newick_string = re.sub(
        r"(?<=[:)])\d+",
        lambda _: str(round(random.uniform(max_distance), 2)),
        newick_string)
    return newick_string
    
def create_output_directory(newick_taxa, path, file_id):
    """
    Given the completely processed newick string creates a directory with at specified path where a folder with the 
    image and the newick string and a specified file ID that is appended to the end of the directory name

    Args:
        newick_taxa (_type_): _description_
        path (_type_): _description_
        file_id (_type_): _description_
    """
    if not os.path.exists(path):
        os.makedirs(path)
    tree_to_image_file(newick_taxa, path, file_id)
    newick_path = path + f"\\newick_{file_id}.nwk"
    nwk_file = open(newick_path, "w")
    nwk_file.write(newick_taxa)
    nwk_file.close() 

def main():
    argument_parser = argparse.ArgumentParser(
        description="Data Collection for Extracting phylogenies from images using AI"
    )
    #
    # Arguments
    #
    
    # parameter for the preferred number of taxa generated, if not specified defaults to 10 taxa
    argument_parser.add_argument('-a', '--amount_taxa', required=False, type=int,
                                 help='Type=Int. Choose the preferred amount of generated taxa. If not specified defaults to 10 taxa.')
    argument_parser.add_argument('-ra', '--randomize_amount', required=False, action='store_true',
                                 help='On/Off flag. Randomizes the amount of taxa generated. Range: 2 to <amount_taxa>. Per default the amount of taxa is not randomized.')
    argument_parser.add_argument('-rd', '--randomize_distances', required=False, type=int,
                                 help='Type=Int. If an Integer is specified the distances between taxa are randomized between 0.00 and the specified number. Per default all distances are set to 1.')
    
    # Specified parameters
    args = argument_parser.parse_args()
    
    #
    # Main functionality
    #
    
    # if the amount of taxa is not specified it defaults to 10 taxa
    amount_taxa = args.amount_taxa
    randomize_amount = args.randomize_amount
    max_distance = args.randomize_distances
    
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
    
    # randomize the distances if a max distance is specified
    if max_distance != None:
        newick_taxa = randomize_distances(newick_taxa, max_distance)
    
    # 
    # Generation of output
    #
    time = str(datetime.datetime.now().strftime(r"%Y-%m-%d-%H-%M-%S-%f"))
    print("Data generation for extracting phylogenies from images using AI.")
    print("Run at: "+ str(datetime.datetime.now()))
    print("File ID: "+ time)
    print(f"Parameters:\n  Randomize distances: {str(False if max_distance == None else True)}\n  Max distance: {max_distance if max_distance != None else 1}")
    print(f"  Specified amount of taxa: {amount_taxa}\n  Actual amount of taxa: {len(random_taxids)}")
    print(f"Newick tree:\n{newick_taxa}")
    file_id = time
    path = str(pathlib.Path(__file__).parent.resolve()) + f"\\generated_data\\data_{time}"
    create_output_directory(newick_taxa, path, file_id)
    
# execute the main method
main()


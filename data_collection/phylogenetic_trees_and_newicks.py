from Bio import Phylo
from ete3 import NCBITaxa
from ete3 import Tree
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
# Formatting issue? Each clade receives a 1
def generate_newick_tree(valid_rand_taxa):
    newick = ncbi.get_topology(valid_rand_taxa).write()
    newick = re.sub(r"(?<=\))\d", "", newick)
    return newick

def save_newick_image(newick, file_name = None, outfile_path = None, display_branch_lengths = True):
    """
    Given a newick tree as a string, generates an image of the newick tree and saves it to the specified path. If no path is specified it will instead 
    save the image to the current directory. If no file name is provided the image in the current directory will be called "tree.jpg".
    If the path contains directories that dont exist, they will be created. If no valid file extension is used, the image will be a .jpg.

    Args:
        newick (str): newick string of the phylogenetic tree to be created
        outfile_path (str): path to where the image is saved at
        file_name (str): Name of the file if no outfile_path was specified
    """
    # make the the newick tree (string) into a file so that it can be drawn by bio.phylo 
    newick_tree = Phylo.read(StringIO(newick), "newick")
    # create a matplotlib figure and 
    fig = plt.figure(figsize=(30, 20), dpi=150)
    axes = fig.add_subplot(1, 1, 1)
    newick_tree.rooted = True
    plt.axis("off")
    if display_branch_lengths:
        Phylo.draw(newick_tree, axes=axes, do_show=False, branch_labels=lambda c: c.branch_length)
    else: 
        Phylo.draw(newick_tree, axes=axes, do_show=False, branch_labels=lambda c: None)
    if outfile_path:
        if not os.path.exists(os.path.dirname(outfile_path)):
            os.makedirs(os.path.dirname(outfile_path), exist_ok=True)
            plt.axis("off")
            plt.savefig(outfile_path, bbox_inches='tight') 
            print(f"[{datetime.datetime.now()}] save_newick_image: Image was saved to newly created directory.")
        else: 
            plt.savefig(outfile_path, bbox_inches='tight')
            print(f"[{datetime.datetime.now()}] save_newick_image: Image was saved to specified directory.")
    elif file_name: 
        if not file_name.endswith(".jpg") or file_name.endswith(".png") or file_name.endswith(".jpeg") or file_name.endswith(".pdf") or file_name.endswith(".svg"):
            plt.axis("off")
            plt.savefig(file_name + ".jpg", bbox_inches='tight')
            print(f"[{datetime.datetime.now()}] save_newick_image: Image was saved to specified file.")
        else: 
            plt.axis("off")
            plt.savefig(file_name, bbox_inches='tight')
            print(f"[{datetime.datetime.now()}] save_newick_image: Image was saved to specified file.")
    else:
        plt.axis("off")
        plt.savefig(f"tree.jpg", bbox_inches='tight')
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
    Given a newick tree and a max distance returns the newick tree with randomized distances. 
    Randomized distances are rounded to 2 decimal and are in the range from 0.00 to the specified max amount.

    Args:
        newick_string (str): The tree in newick format
        max_distance (float): upper limit of randomized distances
    """
    # if re.sub() is given a lambda it tries to give the match object to the lambda but the lambda doesnt take any inputs
    # and instead generates a new random float each time a match is found
    
    newick_string = re.sub(
        r"(?<=:)\d+",
        lambda _: str(round(random.uniform(max_distance), 2)),
        newick_string)
    return newick_string
    
def create_output_directory(newick_taxa, outdirectory_path, file_id, display_branch_lengths):
    """
    Given the newick with taxa (not taxIDs) string creates a directory with at specified path where a folder with the 
    image and the newick string and a specified file ID that is appended to the end of the directory name

    Args:
        newick_taxa (str): newick string with taxa and distances 
        outdirectory_path (str): path to the directory to be created
        file_id (str): chosen file ID 
    """
    if not os.path.exists(outdirectory_path):
        os.makedirs(outdirectory_path)
        print(f"[{datetime.datetime.now()}] create_output_directory: Directory wasn't found. New directory was created.")
    newick_path = outdirectory_path + f"\\newick_{file_id}.nwk"
    # save the newick image
    if display_branch_lengths:
        save_newick_image(newick=newick_taxa, outfile_path=outdirectory_path + f"\\image_tree_{file_id}.jpg")
    else:
        save_newick_image(newick=newick_taxa, outfile_path=outdirectory_path + f"\\image_tree_{file_id}.jpg", display_branch_lengths=False)
    # save the .nwk
    with open(newick_path, "w") as nwk_file:
        nwk_file.write(newick_taxa)
    print(f"[{datetime.datetime.now()}] Newick string was saved into .nwk.")

def main():
    if __name__ == '__main__': 
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
        argument_parser.add_argument("-db", "--display_branch_lengths", required=False, action="store_true",
                                     help="On/Off flag. If specified then all branch lengths are added to the corresponding branch in the image. Per default all branch lengths are added.")
        argument_parser.add_argument("-o", "--outfile_path", required=False, type=str,
                                     help="""If specified the directory containing the Newick and the corresponding image will be saved to this path. If the path points to a directory that doesn't exist
                                     it will be created.""")
        
        # Specified parameters
        args = argument_parser.parse_args()
        
        #
        # Main functionality
        #
        
        # if the amount of taxa is not specified it defaults to 10 taxa
        amount_taxa = args.amount_taxa
        randomize_amount = args.randomize_amount
        max_distance = args.randomize_distances
        display_branch_lengths= args.display_branch_lengths
        outfile_path = args.outfile_path
        
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
        print(f"""Data generation for extracting phylogenies from images using AI.\n
              Default file ID: {time}\n
              Parameters:\n
              \tRandomize distances: {str(False if max_distance == None else True)}\n
              \tMax distance: {max_distance if max_distance != None else 1}\n
              \tSpecified amount of taxa: {amount_taxa}\n
              \tActual amount of taxa: {len(random_taxids)}
              Newick tree:\n
              \t{newick_taxa}
              """)
        file_id = time
        # save directory to outfile path if it was specified. If not create the generated data directory with the data directory inside
        if outfile_path:
            if display_branch_lengths:
                create_output_directory(newick_taxa, outfile_path, file_id, display_branch_lengths=True)
            else:
                create_output_directory(newick_taxa, outfile_path, file_id, display_branch_lengths=False)
        else:
            path = str(pathlib.Path(__file__).parent.resolve()) + f"\\generated_data\\data_{time}"
            if display_branch_lengths:
                create_output_directory(newick_taxa, path, file_id, display_branch_lengths=True)
            else:
                create_output_directory(newick_taxa, path, file_id, display_branch_lengths=False)
    
# execute the main method
main()


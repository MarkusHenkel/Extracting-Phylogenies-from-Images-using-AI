from Bio import Phylo
from ete3 import NCBITaxa
from ete3 import Tree, TreeStyle, NodeStyle
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

def is_newick(newick, format=1):
    """
    Given a newick string checks if it has valid formatting.

    Args:
        newick (str): Newick string
        format (int): format of the newick according to the ETE3 toolki

    Returns:
        bool: True only if the formatting is valid
    """
    try:
        if format:
            tree = Tree(newick, format=format)
        else:
            tree = Tree(newick)
    except:
        return False
    return True

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
    newick = ncbi.get_topology(valid_rand_taxa).write()
    newick = re.sub(r"(?<=\))\d", "", newick)
    return newick

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
    newick_string = re.sub(
        r"(?<=:)\d+",
        lambda _: str(round(random.uniform(max_distance), 2)),
        newick_string)
    return newick_string

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

class TreeRender:
    """
    Class for creating Newick + image and writing a directory containing both into a specified directory.
    """
    def __init__(
        self, 
        newick, 
        package,
        file_id = None,
        outdir_path = None, # Path to the directory where output directory is created at
        display_branch_lengths = True, 
        circular_tree = False, 
        right_to_left_orientation = False, # if False then default orientation left to right is applied
        dont_allow_multifurcations = True,
        branch_vertical_margin = 10 # number of pixels between adjacent branches
        ):
        # Attributes
        self.newick = newick
        self.package = package
        self.file_id = file_id
        self.outdir_path = outdir_path
        self.display_branch_lengths = display_branch_lengths
        self.circular_tree = circular_tree
        self.right_to_left_orientation = right_to_left_orientation
        self.dont_allow_multifurcations = dont_allow_multifurcations
        self.branch_vertical_margin = branch_vertical_margin
        
    def save_newick_image(self, outfile_path):
        """
        Given a newick tree as a string, generates an image of the newick tree and saves it to the specified path. 
        If no path is specified it will instead save the image to the current directory. 
        If no file name is provided the image in the current directory will be called "{package}_tree_{id}.jpg".
        If the path contains directories that dont exist, they will be created. If no valid file extension is used, the image will be a .jpg.
        """
        if not self.package in ["phylo", "ete3"]:
            exit(f"[{datetime.datetime.now()}] save_newick_image: Package {self.package} not valid.")
        elif self.package == "phylo":
            self.save_newick_image_phylo(outfile_path)
        elif self.package == "ete3":
            self.save_newick_image_ete3(outfile_path)
            
    def save_newick_image_phylo(self):
        """
        Used on a TreeRender object generates and saves an image of the TreeRender objects Newick using Biopython.Phylo

        Args:
            outfile_path (str): path where the image is saved at 
        """
        # check if the newick has correct formatting
        if not is_newick(self.newick):
            print(f"[{datetime.datetime.now()}] Error in save_newick_image(): Newick string doesn't have valid formatting.")
        # make the the newick tree (string) into a file so that it can be drawn by bio.phylo 
        newick_tree = Phylo.read(StringIO(self.newick), "newick")
        # create a matplotlib figure and 
        fig = plt.figure(figsize=(30, 20), dpi=150)
        axes = fig.add_subplot(1, 1, 1)
        newick_tree.rooted = True
        plt.axis("off")
        if self.display_branch_lengths:
            Phylo.draw(newick_tree, axes=axes, do_show=False, branch_labels=lambda c: c.branch_length)
        else: 
            Phylo.draw(newick_tree, axes=axes, do_show=False, branch_labels=lambda c: None)
        if self.file_id:
            outfile_path = self.outdir_path + f"\\phylo_tree_{self.file_id}.jpg"
        else:
            outfile_path = self.outdir_path + "\\phylo_tree.jpg"
        if outfile_path:
            if not os.path.exists(os.path.dirname(outfile_path)):
                os.makedirs(os.path.dirname(outfile_path), exist_ok=True)
                plt.axis("off")
                plt.savefig(outfile_path, bbox_inches='tight') 
                print(f"[{datetime.datetime.now()}] save_newick_image_phylo: Image was saved to newly created directory.")
            else: 
                plt.savefig(outfile_path, bbox_inches='tight')
                print(f"[{datetime.datetime.now()}] save_newick_image_phylo: Image was saved to specified directory.")
        else:
            plt.axis("off")
            plt.savefig(f"phylo_tree_{self.file_id}.jpg", bbox_inches='tight')
        plt.close()
        
    def save_newick_image_ete3(self, outfile_path):
        """
        Used on a TreeRender object generates and saves an image of the TreeRender objects Newick using the ETE3 toolkit 

        Args:
            outfile_path (str): path where the image is saved at
        """
        if not is_newick(self.newick):
            exit(f"[{datetime.datetime.now()}] Error in save_newick_image_ete3: Given newick isn't valid.")
        newick_tree = Tree(self.newick)
        treestyle = TreeStyle()
        # TODO: would it be smart to test if the AI is able to give leaves trivial names from top to bottom if the given 
        # tree doesnt have taxa? 
        # currently all leaf names are shown
        treestyle.show_leaf_name = True
        # TODO: depending on how common this is it might be interesting to test if the AI can differentiate both of the values
        treestyle.show_branch_support = False # branch support currently isnt shown 
        treestyle.show_branch_length = self.display_branch_lengths
        treestyle.scale = 50 # 50 pixels per branch length unit
        # default margin is 10 pixels
        if self.branch_vertical_margin:
            if self.branch_vertical_margin < 10:
                treestyle.branch_vertical_margin = 10
                print(f"[{datetime.datetime.now()}] Warning in save_newick_image_ete3: branch_vertical_margin should not deceed 5 pixels. Defaulting to 10 pixels.")
            else:
                treestyle.branch_vertical_margin = self.branch_vertical_margin  
        else:
            treestyle.branch_vertical_margin = 10
        if self.circular_tree:
            treestyle.mode = "c"
        else:     
            treestyle.mode = "r"
        # in ete3 0 is left to right, 1 is right to left 
        treestyle.orientation = int(self.right_to_left_orientation)
        # solve multifurcations if they aren't allowed
        if self.dont_allow_multifurcations == True:
            newick_tree.resolve_polytomy(recursive=True)    
        # save the image to a specified path or into the current working directory with a specified name
        if outfile_path:
            if not os.path.exists(os.path.dirname(outfile_path)):
                os.makedirs(os.path.dirname(outfile_path), exist_ok=True)
                newick_tree.render(file_name=outfile_path, tree_style=treestyle)
                print(f"[{datetime.datetime.now()}] save_newick_image_ete3: Image was saved to newly created directory.")
            else: 
                newick_tree.render(file_name=outfile_path, tree_style=treestyle)
                print(f"[{datetime.datetime.now()}] save_newick_image_ete3: Image was saved to specified directory.")
        else:
            newick_tree.render(file_name=f"ete3_tree_{self.file_id}", tree_style=treestyle)
    
    def create_output_directory(self):
        """
        Combines Newick and image generation to create the output directory.
        """
        # check if output directory needs to be created first
        if not self.outdir_path:
            path = str(pathlib.Path(__file__).parent.resolve()) + f"\\generated_data\\data_{self.file_id}\\"
            self.save_newick_image(path + f"{self.package}_tree_{self.file_id}.jpg")
            newick_path = path + f"newick_{self.file_id}.nwk"
            print(f"[{datetime.datetime.now()}] create_output_directory: Image was saved to generated_data directory in working directory.")
        elif os.path.exists(self.outdir_path):
            self.save_newick_image(self.outdir_path + f"\\generated_data\\data_{self.file_id}\\newick_{self.file_id}.nwk")
            newick_path = self.outdir_path + f"\\generated_data\\data_{self.file_id}\\newick_{self.file_id}.nwk"
            print(f"[{datetime.datetime.now()}] create_output_directory: Image was saved to specified directory.")
        else:
            os.makedirs(self.outdir_path)
            newick_path = self.outdir_path + f"\\generated_data\\data_{self.file_id}\\newick_{self.file_id}.nwk"
            print(f"[{datetime.datetime.now()}] create_output_directory: Directory doesn't exist. Specified directory was created.")
        # save the .nwk
        with open(newick_path, "w") as nwk_file:
            nwk_file.write(self.newick)
        print(f"[{datetime.datetime.now()}] create_output_directory: Newick string was saved into .nwk.")

def main():
    # prevent execution of the whole module when calling it in other modules
    if __name__ == '__main__': 
        argument_parser = argparse.ArgumentParser(
            description="Data Collection for Extracting phylogenies from images using AI"
        )
        #
        # Arguments
        #
        
        # parameter for the preferred number of taxa generated, if not specified defaults to 10 taxa
        argument_parser.add_argument("-p", "--package", required=True, choices=["phylo", "ete3"], 
                                 help="Specify which package is used in the creation of the image. Choose between Biopython.Phylo and ETE3 Toolkit.")
        argument_parser.add_argument('-a', '--amount_taxa', required=False, type=int, default=10,
                                    help='Type=Int. Choose the preferred amount of generated taxa. Default: 10.')
        argument_parser.add_argument('-ra', '--randomize_amount', required=False, action='store_true', default=False,
                                    help='On/Off flag. Randomizes the amount of taxa generated. Range: 2 to <amount_taxa>. Default: False.')
        argument_parser.add_argument('-rd', '--randomize_distances', required=False, type=int, default=None,
                                    help='Type=Int. If an Integer is specified the distances between taxa are randomized between 0.00 and the specified number. Default: 1, not randomized.')
        argument_parser.add_argument("-db", "--display_branch_lengths", required=False, action="store_true", default=True,
                                     help="On/Off flag. If specified then all branch lengths are added to the corresponding branch in the image. Default: True")
        argument_parser.add_argument("-o", "--outdir_path", required=False, type=str,
                                     help="""If specified the directory containing the Newick and the corresponding image will be saved to this path. If the path points to a directory that doesn't exist
                                     it will be created, path can't point to a file. If not specified then a generated_data directory will be created in the current working directory.""")
        argument_parser.add_argument("-c", "--circular_tree", required=False, action="store_true", default=False,
                                     help="On/Off flag. ETE3 only. If True the tree in the image will be in circular format. Default: False.")
        argument_parser.add_argument("-rl", "--right_to_left_orientation", required=False, action="store_true", default=False,
                                     help="On/Off flag. ETE3 only. Specify if tree is oriented from right to left (taxa on the left). Default: False (left to right orientation)")
        argument_parser.add_argument("-da", "--dont_allow_multifurcations", required=False, action="store_true", default=False,
                                     help="On/Off flag. ETE3 only. If True then the resulting tree will be binary. Default: False.")
        argument_parser.add_argument("-vm", "--branch_vertical_margin", required=False, type=int, default=10,
                                     help="ETE3 only. Amount of pixels between two adjacent branches. Should not be smaller than 5. Default: 10 pixels.")
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
        outdir_path = args.outdir_path
        package = args.package
        circular_tree = args.circular_tree
        dont_allow_multifurcations = args.dont_allow_multifurcations
        right_to_left_orientation = args.right_to_left_orientation 
        branch_vertical_margin = args.branch_vertical_margin
                
        # if amount of taxa and randomize_amount are specified then call generate_random_taxids() with the amount of taxa and randomize set to True
        # TODO: remove this line?
        random_taxids = []
        # -a and -r specified
        if amount_taxa != None and randomize_amount == True:
            random_taxids = generate_random_taxids(amount=amount_taxa, randomize=True)
            newick_with_taxids = generate_newick_tree(random_taxids)
            newick_with_taxa = translate_newick_tree(newick_with_taxids)
        # -a specified, -r=False
        elif amount_taxa != None and randomize_amount == False:
            random_taxids = generate_random_taxids(amount=amount_taxa)
            newick_with_taxids = generate_newick_tree(random_taxids)
            newick_with_taxa = translate_newick_tree(newick_with_taxids)
        # -r specified, -a=None
        elif amount_taxa == None and randomize_amount == True:
            random_taxids = generate_random_taxids(randomize=True)
            newick_with_taxids = generate_newick_tree(random_taxids)
            newick_with_taxa = translate_newick_tree(newick_with_taxids)
        # neither specified, -a=None, -r=False
        else:
            random_taxids = generate_random_taxids()
            newick_with_taxids = generate_newick_tree(random_taxids)
            newick_with_taxa = translate_newick_tree(newick_with_taxids)
        # randomize the distances if a max distance is specified
        if max_distance != None:
            newick_with_taxa = randomize_distances(newick_with_taxa, max_distance)
        # 
        # Generation of output
        #
        file_id = str(datetime.datetime.now().strftime(r"%Y-%m-%d-%H-%M-%S-%f"))
        print("Data generation for extracting phylogenies from images using AI.")
        print(f"Default file ID: {file_id}")
        print("Parameters:")
        print(f"  Randomize distances: {str(False if max_distance == None else True)}")
        print(f"  Max distance: {max_distance if max_distance != None else 1}")
        print(f"  Specified amount of taxa: {amount_taxa}")
        print(f"  Actual amount of taxa: {len(random_taxids)}")
        print(f"  Circular tree: {circular_tree}")
        print(f"  Used package: {package}")
        print(f"  Orientation: {"left to right" if not right_to_left_orientation else "right to left"}")
        print(f"  Don't allow multifurcations: {dont_allow_multifurcations}")
        print(f"  Vertical margin for adjacent branches: {branch_vertical_margin}")
        print("Newick tree:")
        print(f"  {newick_with_taxa}")
        # instantiate treerender object
        tree_render = TreeRender(
            newick=newick_with_taxa,
            outdir_path=outdir_path,
            file_id=file_id,
            package=package,
            display_branch_lengths=display_branch_lengths,
            circular_tree=circular_tree,
            right_to_left_orientation=right_to_left_orientation,
            dont_allow_multifurcations=dont_allow_multifurcations,
            branch_vertical_margin=branch_vertical_margin
        )
        # save directory to outfile path if it was specified. If not create the generated data directory with the data directory inside
        tree_render.create_output_directory()
    
# execute the main method
main()


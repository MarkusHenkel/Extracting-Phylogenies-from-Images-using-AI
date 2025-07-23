from Bio import Phylo
from ete3 import NCBITaxa
from ete3 import Tree, TreeStyle, NodeStyle, AttrFace, faces
from io import StringIO 
from numpy import random
import random
import datetime
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import re
import os
import pathlib

# ncbi taxonomy
ncbi = NCBITaxa()

# time for logs
time = datetime.datetime.now().strftime("%Y-%b-%d %H:%M:%S")

def is_taxid_valid(taxid):
    """
    Check if a given taxon ID is valid (present in the NCBI taxonomy and not a digit string) 

    Args:
        taxid (int): taxon ID from the NCBI taxonomy

    Returns:
        bool: True if taxon ID can be found in the NCBI taxonomy
    """
    print(f"is_taxid_valid: taxid: {taxid}")
    if ncbi.get_taxid_translator([taxid]) and not ncbi.get_taxid_translator([taxid])[taxid].isdigit():
        print(f"is_taxid_valid: ncbi: {ncbi.get_taxid_translator([taxid])}")
        return True
    else:
        return False
# Counts amount of invalid taxids in the first 10000
# Function for testing
# def amount_invalid():
#     counter = 0
#     for taxid in range(10000):
#         if not is_taxid_valid(taxid):
#             counter += 1
#     return counter

def is_newick(newick, format=None):
    """
    Given a newick string checks if it has valid formatting.

    Args:
        newick (str): Newick string
        format (int): format of the newick according to the ETE3 toolkit

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

def remove_support_vals(newick):
    """
    Given a newick removes the support values after closing parentheses.
    i.e. ((A:1,B:1)1:1,(C:1,D:1)1:1); => ((A:1,B:1):1,(C:1,D:1):1);
    
    Args:
        newick (str): newick 

    Returns:
        str: newick without support values
    """
    return re.sub(r"(?<=\))\d", "", newick)

def remove_special_chars(taxon):
    r"""
    Deletes special characters that cause issues with the newick format or the packages from a given newick.
    Warning: Remove special characters after getting the topology otherwise the taxon might not be found in the 
    NCBI taxonomy 

    Args:
        taxon (str): taxon that may contain special character 

    Returns:
        str: taxon without special characters [](),.;:'"\
    """
    return re.sub(r"[\[\]\,\.\;\:\'\"\\\(\)]", "", taxon) 

# Debug
test = "(test1)_[test2]_\"ahaha\"_\'ahaha;:.,\'"
print(remove_special_chars(test))

def generate_newick(amount_taxa=10):
    """
    Generates a newick tree with randomized taxa and the specified amount of taxa using the NCBI taxonomy.
    The topology of the newick corresponds to that of the minimal pruned NCBI taxonomy tree i.e. a taxonomic tree 
    based on hierarchical classification.  

    Args:
        amount_taxa (int): amount of taxa in the output newick

    Returns:
        str: newick string with taxa and without NCBIs support values and without taxa containing special characters
    """
    valid_taxids = []
    valid_taxa = []
    # in any case let the minimum amount of taxa be 2
    if (amount_taxa <= 1):
        amount_taxa = 2
    # Add a valid taxid <amount> times
    i = 0
    while i < amount_taxa:
        # not every taxid in the range 1 to 10000 is valid, we have to manually check each time a taxid is randomized
        taxid = random.randint(1,10000)
        # if the taxid is valid and not already in the taxid list save taxid and the corresponding taxon
        if is_taxid_valid(taxid) and not taxid in valid_taxids:
            print(f"taxid {i}: {taxid}")
            valid_taxids.append(taxid)
            taxon = ncbi.get_taxid_translator([taxid])[taxid]
            print(f"taxon {i}: {taxon}")
            valid_taxa.append(taxon)
            i += 1
    # get the topology of the taxids
    newick_taxids = ncbi.get_topology(valid_taxids).write()
    print("with support"+newick_taxids)
    # remove support values that get_topology adds
    newick = remove_support_vals(newick_taxids) 
    print("amount taxa "+ str(amount_taxa))
    print("amount taxids "+ str(len(valid_taxids)))
    print("amount taxa "+ str(len(valid_taxa)))
    print("taxids:" + str(valid_taxids))
    print("taxids:" + str(valid_taxa))
    print("without support"+newick)
    # special characters that are removed because they could interfere with newick parsing or packages 
    for i in range(len(valid_taxids)):
        # translate taxids to their corresponding taxa, replace spaces with _ for newick parsing reasons
        # newick = newick.replace(str(valid_taxids[i]), valid_taxa[i].replace(" ", "_"))
        newick = re.sub(fr"(?<!\d){valid_taxids[i]}(?!\d)", valid_taxa[i].replace(" ", "_"))
        print(newick) 
        # remove special chars 
        newick = newick.replace(valid_taxa[i], remove_special_chars(valid_taxa[i]))
    print(newick)
    return newick
# generate_newick()
# exit()
def randomize_distances_func(newick_string, max_distance):
    """
    Given a newick tree and a max distance returns the newick tree with randomized distances. 
    Randomized distances are rounded to 2 decimal and are in the range from 0.00 to the specified max amount.

    Args:
        newick_string (str): The tree in newick format
        max_distance (float): upper limit of randomized distances
    """
    newick_string = re.sub(
        r"(?<=:)\d+",
        lambda _: str(round(random.uniform(0.00, max_distance), 2)),
        newick_string)
    return newick_string

def taxids_from_newick(newick_string):
    """
    Given a newick tree containing taxon IDs return a list of all taxon IDs inside the string using regex.

    Args:
        newick_string (str): The tree in newick format

    Returns:
        List[str]: List containing all taxon IDs in order
    """
    return re.findall(r"[^)(,:][\d]+[^:]", newick_string)

class TreeRender:
    """
    An instance of the TreeRender class holds all values needed for creating the newick, image and the tsv
    containing the used parameters. 
    """
    def __init__(
        self, 
        newick, 
        fontsize,
        linewidth,
        randomize_distances, # attributes of newick string itself not the render
        max_distance, # attributes of newick string itself not the render
        amount_taxa, # attributes of newick string itself not the render
        package,
        file_id = None,
        outdir_path = None, # Path to the directory where output directory is created at
        dont_display_lengths= False, 
        circular_tree = False, 
        right_to_left_orientation = False, # if False then default orientation left to right is applied
        dont_allow_multifurcations = True,
        branch_vertical_margin = 10, # number of pixels between adjacent branches
        hierarchy_only = False # equivalent to removing distances in newick + dont_display_lengths = True
        ):
        # Attributes
        self.newick = newick
        self.randomize_distances = randomize_distances
        self.max_distance = max_distance
        self.amount_taxa = amount_taxa
        self.package = package
        self.file_id = file_id
        self.outdir_path = outdir_path
        self.dont_display_lengths = dont_display_lengths
        self.circular_tree = circular_tree
        self.right_to_left_orientation = right_to_left_orientation
        self.dont_allow_multifurcations = dont_allow_multifurcations
        self.branch_vertical_margin = branch_vertical_margin
        self.fontsize = fontsize
        self.linewidth = linewidth
        self.hierarchy_only = hierarchy_only 
    
    def is_multifurcating(self):
        """
        Checks if the newick of a treerender object contains a multifurcation

        Returns:
            bool: True if the newick contains a multifurcation
        """
        return self.is_multifurcating_worker(Tree(self.newick))
    
    def is_multifurcating_worker(self, newick_tree):
        if len(newick_tree.get_children()) > 2:
            return True
        for node in newick_tree.get_children():
            if self.is_multifurcating_worker(node):
                return True
        return False   
    
    def solve_multifurcations(self):
        """
        Solves multifurcations in the newick and the tree object using ete3's resolve_polytomy function
        """
        if not self.is_multifurcating:
            # Debugging
            # print("No multifurcations found.")
            return
        # print(f"Current newick: {self.newick}")
        newick_tree = Tree(self.newick)
        newick_tree.resolve_polytomy(recursive=True)    
        # update the newick, so that the tree in the image and the newick match and remove support values by ete3
        self.newick = remove_support_vals(newick_tree.write())
        # print(f"Updated newick: {self.newick}")
        # ete solves polytomies by adding branches with length zero, replace them with the specified branch lengths
        if self.randomize_distances == True:
            self.newick = re.sub(r"(?<=\):)0(?=[,\)])",
                                    lambda m: str(round(random.uniform(0.00, self.max_distance), 2)), 
                                    self.newick)
        else:
            self.newick = re.sub(r"(?<=\):)0(?=[,\)])", "1", self.newick)
        
    def save_newick_image(self, outfile_path):
        """
        Given a newick tree as a string, generates an image of the newick tree and saves it to the specified path. 
        If no path is specified it will instead save the image to the current directory. 
        If no file name is provided the image in the current directory will be called "{package}_tree_{id}.jpg".
        If the path contains directories that dont exist, they will be created. If no valid file extension is used, 
        the image will be a .jpg.
        """
        if not self.package in ["phylo", "ete3"]:
            exit(f"[{time}] save_newick_image: Package {self.package} not valid.")
        elif self.package == "phylo":
            self.save_newick_image_phylo(outfile_path)
        elif self.package == "ete3":
            self.save_newick_image_ete3(outfile_path)
            
    def save_newick_image_phylo(self, outfile_path):
        """
        Used on a TreeRender object generates and saves an image of the TreeRender objects Newick using Biopython.Phylo

        Args:
            outfile_path (str): path where the image is saved at 
        """
        # check if the newick has correct formatting
        if not is_newick(self.newick):
            print(f"[{time}] Error in save_newick_image(): Newick string doesn't have valid formatting.")
        # make the the newick tree (string) into a file so that it can be drawn by bio.phylo 
        newick_tree = Phylo.read(StringIO(self.newick), "newick")
        # create a matplotlib figure 
        fig = plt.figure(figsize=(30, 20), dpi=150)
        axes = fig.add_subplot(1, 1, 1)
        # set line width and fontsize
        mpl.rcParams['lines.linewidth'] = self.linewidth # reasonable range: [1,10]
        mpl.rcParams['font.size'] = self.fontsize # reasonable range: [8,30]
        newick_tree.rooted = True
        plt.axis("off")
        if self.dont_display_lengths or self.hierarchy_only:
            Phylo.draw(newick_tree, axes=axes, do_show=False, branch_labels=lambda c: None)
        else: 
            Phylo.draw(newick_tree, axes=axes, do_show=False, branch_labels=lambda c: c.branch_length)
        ###### DEBUGGING
        # plt.show() # show the tree instead of writing a file each time
        if outfile_path:
            if not os.path.exists(os.path.dirname(outfile_path)):
                os.makedirs(os.path.dirname(outfile_path), exist_ok=True)
                plt.savefig(outfile_path, bbox_inches='tight') 
                print(f"[{time}] save_newick_image_phylo: Image was saved to newly created directory.")
            else: 
                plt.savefig(outfile_path, bbox_inches='tight')
                print(f"[{time}] save_newick_image_phylo: Image was saved to specified directory.")
        else:
            plt.savefig(f"phylo_tree_{self.file_id}.jpg", bbox_inches='tight')
        plt.close()
        
    def save_newick_image_ete3(self, outfile_path):
        """
        Used on a TreeRender object generates and saves an image of the TreeRender objects Newick using the ETE3 toolkit 

        Args:
            outfile_path (str): path where the image is saved at
        """
        if not is_newick(self.newick):
            exit(f"[{time}] Error in save_newick_image_ete3: Given newick isn't valid.")
        # create ete3 Tree object 
        newick_tree = Tree(self.newick)
        # create treestyle object to adjust tree properties of Tree object
        treestyle = TreeStyle()
        # create nodestyle object to adjust node and edge properties of Tree object
        nodestyle = NodeStyle()
        # solve multifurcations if they aren't allowed and update the newick if there is a multifurcation
        if self.dont_allow_multifurcations == True and self.is_multifurcating():
            self.solve_multifurcations()
            print(f"[{time}] Solved multifurcation(s).")
        # adjust line width
        nodestyle["vt_line_width"] = self.linewidth # reasonable range: [1,10]
        nodestyle["hz_line_width"] = self.linewidth # reasonable range: [1,10]
        # apply nodestyle adjustments to the whole tree
        for node in newick_tree.traverse():
            node.set_style(node_style=nodestyle)
        # adjust fontsize 
        fontsize = self.fontsize # reasonable range: [8,30]
        # create layout for adjusting the fontsize of taxa and distances
        def layout(node):
            if node.is_leaf():
                # create face for taxon with chosen fontsize
                taxa_face = faces.AttrFace("name", fsize=fontsize)
                # apply face to the current leaf
                faces.add_face_to_node(taxa_face, node, column=0)
            # if dont_display_lengths or hierarchy_only are True, then dont create faces for distances
            if not (self.dont_display_lengths or self.hierarchy_only):
                if not node.is_root():
                    # create face for the distance 
                    dist_face = faces.AttrFace(f"dist", fsize=fontsize)
                    # apply face to the current node
                    faces.add_face_to_node(dist_face, node, column=0, position="branch-top")
        # set the layout
        treestyle.layout_fn = layout
        # TODO: would it be smart to test if the AI is able to give leaves trivial names from top to bottom if the given 
        # tree doesnt have taxa? 
        treestyle.show_leaf_name = True # currently all leaf names are shown
        treestyle.show_branch_support = False # branch support currently isnt shown 
        # treestyle.scale = 100 # 100 pixels per branch length unit
        # set default margin is 10 pixels
        if self.branch_vertical_margin:
            if self.branch_vertical_margin < 5:
                treestyle.branch_vertical_margin = 10
                print(f"[{time}] Warning in save_newick_image_ete3: branch_vertical_margin should not deceed 10 pixels. \
                    Defaulting to 10 pixels.")
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
        # creating text faces for adjusting fontsize does not replace the default font, turn off the default font 
        treestyle.show_leaf_name = False
        treestyle.show_branch_length = False
        ###### DEBUGGING
        # newick_tree.show(tree_style=treestyle) # show the tree instead of writing a file each time
        # save the image to a specified path or into the current working directory with a specified name
        if outfile_path:
            if not os.path.exists(os.path.dirname(outfile_path)):
                os.makedirs(os.path.dirname(outfile_path), exist_ok=True)
                newick_tree.render(file_name=outfile_path, tree_style=treestyle)
            else: 
                newick_tree.render(file_name=outfile_path, tree_style=treestyle)
        else:
            newick_tree.render(file_name=f"ete3_tree_{self.file_id}", tree_style=treestyle)
            
    def write_params_to_tsv(self, outfile_path):
        """
        Writes the parameters used in the creation of the image into a tsv so that we can compare the performance of 
        AI on different parameters.
        This is done for each image seperately.
        """
        tsv_header = "random_distances\tmax_distance\tamount_taxa\tpackage\tbranch_lengths\tcircular_tree\t"
        tsv_header += "right_to_left_orientation\tmultifurcations\tbranch_vertical_margin[px]\t"
        tsv_header += "fontsize\tlinewidth[px]\thierarchy_only\t"
        tsv_header += "\n"
        params = f"{self.randomize_distances}\t{self.max_distance}\t{self.amount_taxa}\t"
        params += f"{self.package}\t{self.dont_display_lengths}\t{self.circular_tree}\t"
        params += f"{self.right_to_left_orientation}\t{not self.dont_allow_multifurcations}\t"
        params += f"{self.branch_vertical_margin}\t{self.fontsize}\t{self.linewidth}\t{self.hierarchy_only}\t" 
        with open(outfile_path, "w") as tsv_file:
            tsv_file.write(tsv_header)
            tsv_file.write(params)       
        
    def create_output_directory(self):
        """
        Combines Newick and image generation to create the output directory.
        """
        # check if a path to the output directory was specified and if it needs to be created first   
        if self.outdir_path:
            if not os.path.exists(self.outdir_path):
                os.makedirs(self.outdir_path)
                print(f"[{time}] create_output_directory: Specified directory was created")
            image_path = self.outdir_path + f"\\data{self.file_id}\\tree{self.file_id}.jpg"
            newick_path = self.outdir_path + f"\\data{self.file_id}\\newick{self.file_id}.nwk"
            tsv_path = self.outdir_path + f"\\data{self.file_id}\\params{self.file_id}.tsv"
            print(f"[{time}] create_output_directory: Image is saved to specified directory.")
        else:
            path = str(pathlib.Path(__file__).parent.resolve()) + f"\\data{self.file_id}\\"
            image_path = path + f"tree{self.file_id}.jpg"
            newick_path = path + f"newick{self.file_id}.nwk"
            tsv_path = path + f"params{self.file_id}.tsv"
            print(f"[{time}] create_output_directory: Image is saved to generated_data directory in working directory.")
        # save the image
        self.save_newick_image(image_path)
        # update the newick before writing it to file if --hierarchy_only was specified
        if self.hierarchy_only:
            self.remove_distances_from_newick()
        # save the .nwk
        with open(newick_path, "w") as nwk_file:
            nwk_file.write(self.newick)
        print(f"[{time}] create_output_directory: Newick string was saved into .nwk.")
        # save the tsv
        self.write_params_to_tsv(tsv_path)
        print(f"[{time}] create_output_directory: Used parameters were saved into .tsv.")
        
    # TODO: implement create_rand_treerender()
    def randomize_treerender(self):
        """
        When applied to a TreeRender object, randomizes the following parameters:
        package, amount_taxa, randomize_distances, max_distance, dont_allow_multifurcations, branch_vertical_margin,
        fontsize, linewidth
        This way a user that does not want to understand all parameters can create a simple but 
        diverse dataset quickly.

        Returns:
            TreeRender: TreeRender object with user parameters randomized 
        """
        self.package = random.choice(["ete3","phylo"])
        if self.package == "ete3":
            self.dont_allow_multifurcations = random.choice([True,False])
            self.branch_vertical_margin = random.randint(10, 100)
            self.fontsize = random.randint(8,50)
            self.linewidth = random.randint(1,20)
        if self.package == "phylo":
            self.fontsize = random.randint(8,50)
            # TODO: until no solution for increasing the branch length offset is found the linewidth cant be increased
            # for phylo because it omits the branch lengths
            self.linewidth = 1
        self.amount_taxa = random.randint(5,20) # range: [5, 20]
        self.randomize_distances = random.choice([True,False])
        if self.randomize_distances:
            self.max_distance = random.randint(1, 10)
        return self
    
    def remove_distances_from_newick(self):
        """
        Function for creating hierarchy-only tree renders. Removes distances from the treerenders newick using regex. 
        """
        self.newick = re.sub(r":\d+(\.\d+){0,1}", "", self.newick)
        
def ask_user_to_continue():
    while True:
        print("Yes[y]/No[n]?")
        yes_or_no = input()
        if yes_or_no == "n":
            print("Cancelled.")
            exit()
        elif yes_or_no == "y":
            print("Continuing.")
            break  
      
def main():
    # prevent execution of the whole module when calling it in other modules
    if __name__ == '__main__': 
        argument_parser = argparse.ArgumentParser(
            description="Data Collection for Extracting phylogenies from images using AI"
        )
        ########## ARGUMENTS ##########
        argument_parser.add_argument(
            "--create_rand_dataset", 
            required=False, 
            action="store_true", 
            default=False,
            help="""On/Off flag. Quickly create a diverse dataset instead of one type of image by specifying 
            --create_dataset. This will randomize the following parameters and flags within reasonable ranges: 
            package, amount_taxa, randomize_distances, max_distance, dont_allow_multifurcations, branch_vertical_margin, 
            fontsize, linewidth. 
            --number_directories sets the size of the dataset. --hierarchy_only removes distances from the newick and 
            the image of each directory created with --create_rand_dataset."""
        )
        argument_parser.add_argument(
            "--hierarchy_only",
            required=False,
            action="store_true",
            default=False,
            help="""On/Off flag. If specified removes distances from newick and image of the created directory. This is 
            helpful for creating datasets used for an initial training step in which models are just trained to recognize
            the tree's hierarchy and not yet the branch lenghts. Can be combined with --create_rand_dataset. Default:
            False."""
        )
        argument_parser.add_argument("-o", "--outdir_path", required=False, type=str,
                                     help="""If specified the directory containing the Newick and the corresponding 
                                     image will be saved to this path. If the path points to a directory that
                                     doesn't exist it will be created, path can't point to a file. If not specified
                                     then a generated_data directory will be created in the current working 
                                     directory.""")
        argument_parser.add_argument("-n", "--number_directories", type=int, required=False, default=1,
                                     help="""Choose the number of directories created with the chosen parameters. 
                                     Default: 1.""")
        argument_parser.add_argument("-f", "--fontsize", type=int, required=False, 
                                     help="""Choose the size of the font. Also applies to branch lengths if applicable. 
                                     Default: 8 (ETE3), 16 (Bio.Phylo)""")
        argument_parser.add_argument("-l", "--linewidth", type=int, required=False, 
                                     help="Choose the width of branches in pixels. Default: 1 (ETE3), 2 (Bio.Phylo)")
        argument_parser.add_argument("-p", "--package", required=False, choices=["phylo", "ete3"], default="ete3",
                                     help="""Specify which package is used in the creation of the image. 
                                     Choose between Bio.Phylo and ETE3 Toolkit. Default: ete3.""")
        argument_parser.add_argument('-a', '--amount_taxa', required=False, type=int, default=10,
                                     help='Type=Int. Choose the preferred amount of generated taxa. Default: 10.')
        argument_parser.add_argument('-rd', '--randomize_distances', required=False, action="store_true", default=False,
                                     help="""On/Off flag. If specified the distances (branch lengths) will be 
                                     randomized. Default: False.""")
        argument_parser.add_argument('-md', '--max_distance', required=False, type=int, default=1,
                                     help="""Type=Int. If distances are randomized this will be the maximum distance. 
                                     Default: 1.""")
        argument_parser.add_argument("-db", "--dont_display_lengths", required=False, action="store_true", default=False,
                                     help="""On/Off flag. If specified then no branch lengths will be displayed. 
                                     Default: False.""")
        argument_parser.add_argument("-c", "--circular_tree", required=False, action="store_true", default=False,
                                     help="""On/Off flag. ETE3 only. If True the tree in the image will be in 
                                     circular format. Default: False.""")
        argument_parser.add_argument("-rl", "--right_to_left_orientation", required=False, action="store_true", 
                                     default=False,help="""On/Off flag. ETE only. Specify if tree is oriented from 
                                     right to left (taxa on the left). Default: False (left to right orientation)""")
        argument_parser.add_argument("-da", "--dont_allow_multifurcations", required=False, action="store_true", 
                                     default=False,help="""On/Off flag. ETE3 only. If specified then the resulting 
                                     tree will be binary. Default: False.""")
        argument_parser.add_argument("-vm", "--branch_vertical_margin", required=False, type=int, 
                                     help="""ETE3 only. Amount of pixels between two adjacent branches. 
                                     Should not be smaller than 5. Default: 10 pixels.""")
        # Specified parameters
        args = argument_parser.parse_args()
        # if the amount of taxa is not specified it defaults to 10 taxa
        number_directories = args.number_directories
        amount_taxa = args.amount_taxa
        randomize_distances = args.randomize_distances
        max_distance = args.max_distance
        dont_display_lengths= args.dont_display_lengths
        outdir_path = args.outdir_path
        package = args.package
        circular_tree = args.circular_tree
        dont_allow_multifurcations = args.dont_allow_multifurcations
        right_to_left_orientation = args.right_to_left_orientation 
        branch_vertical_margin = args.branch_vertical_margin
        create_rand_dataset = args.create_rand_dataset
        fontsize = args.fontsize
        linewidth = args.linewidth
        hierarchy_only = args.hierarchy_only
        # set the default values of fontsize and linewidth
        if not fontsize:
            if package == "ete3":
                fontsize = 8
            elif package == "phylo":
                fontsize = 16
        if not linewidth:
            if package == "ete3":
                linewidth = 1
            elif package == "phylo":
                linewidth = 2
        ########## CHECKS ##########
        # negative or zero linewidth is set to 1
        if package == "ete3" and linewidth <= 0:
            linewidth = 1
        # warn user if they use ete3 parameters with a module other than ete3
        if (not package == "ete3") and (
            branch_vertical_margin or 
            dont_allow_multifurcations or 
            circular_tree or 
            right_to_left_orientation
        ):
            # warn user which parameters cant be applied and reset them to default values
            print(f"[{time}] Warning: You are using parameters that are ete3 only: ")
            if branch_vertical_margin:
                print("branch_vertical")
                branch_vertical_margin = None
            if dont_allow_multifurcations:
                print("dont_allow_multifurcations") 
                dont_allow_multifurcations = False 
            if circular_tree:
                print("circular_tree") 
                circular_tree = False
            if right_to_left_orientation:
                print("right_to_left_orientation") 
                right_to_left_orientation = False
            print(f"[{time}] Resetting ete3 parameters. Do you want to proceed?")
            ask_user_to_continue()
        if number_directories < 1:
            exit(f"[{time}] --number_directories {number_directories} not valid. Has to be positive integer > 1.")
        # warn user if he wants to create more than one image
        if number_directories > 1:
            print(f"[{time}] Warning: You are about to create {number_directories} directories. Are you sure you want to proceed?")
            ask_user_to_continue()
        ########## MAIN LOOP ##########
        print("Data generation for extracting phylogenies from images using AI.")
        # execute module <number_directories> times
        for i in range(number_directories):
            # unique file ID is just the current time (hour, minute, second, microsecond)
            # could be shortened to second and microsecond maybe
            file_id = str(datetime.datetime.now().strftime(r"%H%M%S%f"))
            ########## INSTANTIATING TREERENDER OBJECT ##########
            tree_render = TreeRender(
                newick="", # initialize with empty string, finished newick is added at a later point
                randomize_distances = randomize_distances,
                max_distance=max_distance,
                amount_taxa=amount_taxa,
                outdir_path=outdir_path,
                file_id=file_id,
                package=package,
                dont_display_lengths=dont_display_lengths,
                circular_tree=circular_tree,
                right_to_left_orientation=right_to_left_orientation,
                dont_allow_multifurcations=dont_allow_multifurcations,
                branch_vertical_margin=branch_vertical_margin,
                fontsize=fontsize,
                linewidth=linewidth,
                hierarchy_only=hierarchy_only,
            )
            # if the user wants to generate a dataset with randomized parameters
            if create_rand_dataset:
                tree_render.randomize_treerender()
            ########## CREATE THE NEWICK ##########
            # TODO: fix double _ issue
            newick = generate_newick(tree_render.amount_taxa)
            # randomize the distances if a max distance is specified
            if tree_render.randomize_distances:
                newick = randomize_distances_func(newick, tree_render.max_distance)
            # set the finished newick as the treerenders newick
            tree_render.newick = newick
            # save directory to outfile path if it was specified
            # If not create the generated data directory with the data directory inside
            tree_render.create_output_directory()
            print(f"Default file ID: {file_id}")
            print("Parameters:")
            print(f"  Fontsize: {tree_render.fontsize}")
            print(f"  Linewidth: {tree_render.linewidth}")
            print(f"  Hierarchy only: {tree_render.hierarchy_only}")
            print(f"  Randomize distances: {tree_render.randomize_distances}")
            print(f"  {f"Max distance: {tree_render.max_distance}" if tree_render.randomize_distances \
                else "Distances: all exactly 1"}")
            print(f"  Amount of taxa: {tree_render.amount_taxa}")
            print(f"  Circular tree: {tree_render.circular_tree}")
            print(f"  Used package: {tree_render.package}")
            print(f"  Orientation: {"left to right" if not tree_render.right_to_left_orientation else "right to left"}")
            print(f"  Don't allow multifurcations: {tree_render.dont_allow_multifurcations}")
            print(f"  Vertical margin for adjacent branches: {tree_render.branch_vertical_margin}") \
                if tree_render.package == "ete3" else None
            print(f"  Branch lengths displayed: {not tree_render.dont_display_lengths}")
            print("Newick:")
            print(f"  {tree_render.newick}")
# execute the main method
main()


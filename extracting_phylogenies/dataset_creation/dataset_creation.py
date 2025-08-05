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
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar # size bar phylo
import matplotlib.font_manager as fm # size bar phylo
from matplotlib.transforms import Bbox # size bar phylo
from extracting_phylogenies.utilities import utilities as ut
import sys
import logging

# create logger
console_logger = logging.getLogger(__name__)
# ncbi taxonomy
ncbi = NCBITaxa()
# ncbi.update_taxonomy_database()

def is_taxid_valid(taxid, max_length=35):
    """
    Check if a given taxon ID is valid (present in the NCBI taxonomy and not a digit string) 

    Args:
        taxid (int): taxon ID from the NCBI taxonomy

    Returns:
        bool: True if taxon ID can be found in the NCBI taxonomy
    """
    if ncbi.get_taxid_translator([taxid]) and not ncbi.get_taxid_translator([taxid])[taxid].isdigit():
        # dont allow absurdly long taxa that could malform the image
        if len(ncbi.get_taxid_translator([taxid])[taxid]) > max_length:
            return False
        return True
    else:
        return False

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
    taxids = []
    # in any case let the minimum amount of taxa be 2
    if (amount_taxa <= 1):
        amount_taxa = 2
    # Add a valid taxid <amount> times
    i = 0
    while i < amount_taxa:
        # not every taxid in the range 1 to 10000 is valid, we have to manually check each time a taxid is randomized
        taxid = random.randint(1,10000)
        # if the taxid is valid and not already in the taxid list save it
        if is_taxid_valid(taxid) and taxid not in taxids:
            taxids.append(taxid)
            i += 1
    # remove support values that get_topology adds
    newick = ut.remove_support_vals(ncbi.get_topology(taxids).write()) 
    # some taxon IDs are replaced because of ncbi.get_topology, update taxids 
    taxids = re.findall(r"(?<=[,\(])\d+(?=:)", newick)
    # substitute taxids with their corresponding taxa, replace spaces with _ for newick parsing reasons
    # dont accidentally sub taxids that are within other taxids e.g. target taxid 71, actual taxid 3712
    for taxid in taxids:
        taxon = ncbi.get_taxid_translator([taxid])[int(taxid)]
        newick = re.sub(fr"(?<!\d){taxid}(?!\d)", ut.remove_special_chars(taxon).replace(" ", "_"), newick) 
    return newick

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
        taxa_only = False, # equivalent to removing distances in newick + dont_display_lengths = True
        topology_only = False,
        align_taxa = False,
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
        self.taxa_only = taxa_only
        self.topology_only = topology_only
        self.align_taxa = align_taxa
    
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
            return
        newick_tree = Tree(self.newick)
        newick_tree.resolve_polytomy(recursive=True)    
        # update the newick, so that the tree in the image and the newick match and remove support values by ete3
        self.newick = ut.remove_support_vals(newick_tree.write())
        # ete solves polytomies by adding branches with length zero, replace them with the specified branch lengths
        if self.randomize_distances == True:
            self.newick = re.sub(r"(?<=\):)0(?=[,\)])",
                                    lambda m: str(round(random.uniform(1.00, self.max_distance), 2)), 
                                    self.newick)
        else:
            self.newick = re.sub(r"(?<=\):)0(?=[,\)])", "1.0", self.newick)
        
    def save_newick_image(self, outfile_path):
        """
        Given a newick tree as a string, generates an image of the newick tree and saves it to the specified path. 
        If no path is specified it will instead save the image to the current directory. 
        If no file name is provided the image in the current directory will be called "{package}_tree_{id}.jpg".
        If the path contains directories that dont exist, they will be created. If no valid file extension is used, 
        the image will be a .jpg.
        """
        if not self.package in ["phylo", "ete3"]:
            console_logger.warning(f"Package {self.package} not valid.")
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
        newick_format = 100 if self.topology_only else 0
        # check if the newick has correct formatting 
        if not ut.is_newick(self.newick, format=newick_format):
            console_logger.warning(f"Newick string doesn't have valid formatting.")
            raise ValueError(f"Newick {self.newick} is not valid.")
        # make the the newick into a file so that it can be drawn by bio.phylo 
        newick_tree = Phylo.read(StringIO(self.newick), "newick")
        # create a matplotlib figure 
        fig = plt.figure(figsize=(30, 20), dpi=150)
        fontprops = fm.FontProperties(size=self.fontsize)
        axes = fig.add_subplot(1, 1, 1)
        # if topo or taxa only is specified dont add the sizebar
        if not (self.topology_only or self.taxa_only):
            # add scalebar to plot manually (bio.phylo does not provide it)
            scalebar = AnchoredSizeBar(
                axes.transData,
                1, '1.0', 'lower left', 
                pad=2,
                color='black',
                frameon=False,
                # size_vertical=1,
                fontproperties=fontprops,
                # bbox_to_anchor=Bbox.from_bounds(0,0,1,1),
                bbox_to_anchor=(0.15, 0),
                bbox_transform=axes.figure.transFigure
            )
            axes.add_artist(scalebar)
        # set line width and fontsize
        mpl.rcParams['lines.linewidth'] = self.linewidth # reasonable range: [1,10] (if no branch lengths)
        mpl.rcParams['font.size'] = self.fontsize # reasonable range: [8,16]
        newick_tree.rooted = True
        plt.axis("off")
        # remove taxa from newick if topology_only
        if self.topology_only:
            self.newick = ut.remove_taxa_from_newick(self.newick)
        if self.dont_display_lengths or self.taxa_only or self.topology_only:
            Phylo.draw(newick_tree, axes=axes, do_show=False, branch_labels=lambda c: None)
        else: 
            Phylo.draw(newick_tree, axes=axes, do_show=False, branch_labels=lambda c: c.branch_length)
        ###### DEBUGGING
        # plt.show() # show the tree instead of writing a file each time
        if not os.path.exists(os.path.dirname(outfile_path)):
            os.makedirs(os.path.dirname(outfile_path), exist_ok=True)
            plt.savefig(outfile_path, bbox_inches='tight') 
            console_logger.info(f"Image was saved to newly created directory.")
        else: 
            plt.savefig(outfile_path, bbox_inches='tight')
            console_logger.info(f"Image was saved to specified directory.")
        plt.close()
        
    def save_newick_image_ete3(self, outfile_path):
        """
        Used on a TreeRender object generates and saves an image of the TreeRender objects Newick using the ETE3 toolkit 

        Args:
            outfile_path (str): path where the image is saved at
        """
        # solve multifurcations if they aren't allowed and update the newick if there is a multifurcation
        if self.dont_allow_multifurcations == True and self.is_multifurcating():
            self.solve_multifurcations()
            console_logger.info(f"Solved multifurcation(s).")
        # create ete3 Tree object of the finished newick, all further modifications are taken on the Tree object not the
        # newick itself
        newick_tree = Tree(self.newick)
        # create treestyle object to adjust tree properties of Tree object
        treestyle = TreeStyle()
        # create nodestyle object to adjust node and edge properties of Tree object
        nodestyle = NodeStyle()
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
            # dont put leaf labels if topology_only is True
            if not self.topology_only:
                if node.is_leaf():
                    # create face for taxon with chosen fontsize
                    taxa_face = faces.AttrFace("name", fsize=fontsize)
                    # apply face to the current leaf, align faces if specified
                    if self.align_taxa:
                        faces.add_face_to_node(taxa_face, node, column=0, position="aligned")
                    else:
                        faces.add_face_to_node(taxa_face, node, column=0)
            # if dont_display_lengths, taxa_only or topology_only are True, then dont create faces for distances
            if not (self.dont_display_lengths or self.taxa_only or self.topology_only):
                if not node.is_root():
                    # create face for the distance 
                    dist_face = faces.AttrFace(f"dist", fsize=fontsize)
                    # apply face to the current node
                    faces.add_face_to_node(dist_face, node, column=0, position="branch-top")
        # set the layout
        treestyle.layout_fn = layout
        # treestyle.scale = 100 # 100 pixels per branch length unit
        # set default margin is 10 pixels
        if self.branch_vertical_margin:
            if self.branch_vertical_margin < 5:
                treestyle.branch_vertical_margin = 10
                console_logger.info(f"branch_vertical_margin should not deceed 10 pixels. Defaulting to 10 pixels.")
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
        # if topo or taxa only is specified remove the sizebar
        if self.topology_only or self.taxa_only:
            treestyle.show_scale = False
        # add dotted lines that ete3 adds to edges when edge is too short for branch length label
        treestyle.complete_branch_lines_when_necessary = True 
        ###### DEBUGGING
        # newick_tree.show(tree_style=treestyle) # show the tree instead of writing a file each time
        # save the rendered image to file
        if not os.path.exists(os.path.dirname(outfile_path)):
            os.makedirs(os.path.dirname(outfile_path), exist_ok=True)
            newick_tree.render(file_name=outfile_path, tree_style=treestyle)
        else: 
            newick_tree.render(file_name=outfile_path, tree_style=treestyle)
            
    def write_params_to_tsv(self, outfile_path):
        """
        Writes the parameters used in the creation of the image into a tsv so that we can compare the performance of 
        AI on different parameters.
        This is done for each image seperately.
        """
        tsv_header = "random_distances\tmax_distance\tamount_taxa\tpackage\tbranch_lengths\tcircular_tree\t"
        tsv_header += "right_to_left_orientation\tmultifurcations\tbranch_vertical_margin[px]\t"
        tsv_header += "fontsize\tlinewidth[px]\ttaxa_only\ttopo_only\talign_taxa\t"
        tsv_header += "\n"
        params = f"{self.randomize_distances}\t{self.max_distance}\t{self.amount_taxa}\t"
        params += f"{self.package}\t{not self.dont_display_lengths}\t{self.circular_tree}\t"
        params += f"{self.right_to_left_orientation}\t{not self.dont_allow_multifurcations}\t"
        params += f"{self.branch_vertical_margin}\t{self.fontsize}\t{self.linewidth}\t{self.taxa_only}\t"
        params += f"{self.topology_only}\t{self.align_taxa}\t" 
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
                console_logger.info(f"create_output_directory: Specified directory was created")
            image_path = self.outdir_path + f"\\data{self.file_id}\\tree{self.file_id}.jpg"
            newick_path = self.outdir_path + f"\\data{self.file_id}\\newick{self.file_id}.nwk"
            tsv_path = self.outdir_path + f"\\data{self.file_id}\\params{self.file_id}.tsv"
            console_logger.info(f"Image is saved to specified directory.")
        else:
            path = str(pathlib.Path(__file__).parent.resolve()) + f"\\data{self.file_id}\\"
            image_path = path + f"tree{self.file_id}.jpg"
            newick_path = path + f"newick{self.file_id}.nwk"
            tsv_path = path + f"params{self.file_id}.tsv"
            console_logger.info(f"Image is saved to generated_data directory in working directory.")
        # remove taxa before removing distances and saving the image 
        if self.package == "phylo" and self.topology_only:
            self.newick = ut.remove_taxa_from_newick(self.newick)
        # save the image
        self.save_newick_image(image_path)
        # remove taxa from topology_only trees before removing distances but after saving img because ete3 doesnt allow
        # empty leaf nodes
        if self.package == "ete3" and self.topology_only:
            self.newick = ut.remove_taxa_from_newick(self.newick)
        # remove distances from taxa and topology_only trees
        if self.taxa_only or self.topology_only:
            self.newick = ut.remove_distances(self.newick)
        # save the .nwk
        with open(newick_path, "w") as nwk_file:
            nwk_file.write(self.newick)
        console_logger.info(f"Newick string was saved into .nwk.")
        # save the tsv
        self.write_params_to_tsv(tsv_path)
        console_logger.info(f"Used parameters were saved into .tsv.")
        
    def randomize_treerender(self, excluded_params):
        """
        When applied to a TreeRender object, randomizes the following parameters:
        package, amount_taxa, randomize_distances, max_distance, dont_allow_multifurcations, branch_vertical_margin,
        fontsize, linewidth
        This way a user that does not want to understand all parameters can create a simple but 
        diverse dataset quickly.

        Args:
            excluded_params (List(str)): List of parameters that arent randomized

        """
        if "package" not in excluded_params:
            self.package = random.choice(["ete3","phylo"])
        if "dont_display_lengths" not in excluded_params:
            self.dont_display_lengths = random.choice([True,False])
        if self.package == "ete3":
            if "dont_allow_multifurcations" not in excluded_params:
                self.dont_allow_multifurcations = random.choice([True,False])
            if "branch_vertical_margin" not in excluded_params:
                self.branch_vertical_margin = random.randint(10, 100)
            if "fontize" not in excluded_params:
                self.fontsize = random.randint(8,20)
            if "linewidth" not in excluded_params: 
                self.linewidth = random.randint(1,15)
            if "align_taxa" not in excluded_params:
                self.align_taxa = random.choice([True,False])
        if self.package == "phylo":
            if "fontsize" not in excluded_params:
                self.fontsize = random.randint(12,20)
            if "linewidth" not in excluded_params:
                # increase linewidth only if no edge labels are present
                if self.topology_only or self.taxa_only or self.dont_display_lengths:
                    self.linewidth = random.randint(1,15)
                else: 
                    # in phylo branch lengths arent legible if linewidth is increased
                    self.linewidth = 1 
        if "amount_taxa" not in excluded_params:
            self.amount_taxa = random.randint(5,20) # range: [5, 20]
        if "randomize_distances" not in excluded_params:
            self.randomize_distances = random.choice([True,False])
        if "max_distance" not in excluded_params:
            if self.randomize_distances:
                self.max_distance = random.randint(1, 10)
        ##### DEBUGGING
        # return self
    
    def randomize_distances_func(self):
        """
        Given a newick tree and a max distance returns the newick tree with randomized distances. 
        Randomized distances are rounded to 2 decimal and are in the range from 0.00 to the specified max amount.

        Args:
            newick_string (str): The tree in newick format
            max_distance (float): upper limit of randomized distances
        """
        # only randomize between 1 and max distance to stop distances from being drawn over edges and being illegible 
        # because their corresponding branch is too short
        if self.max_distance >= 1:
            self.newick = re.sub(
                r"(?<=:)\d+",
                lambda _: str(round(random.uniform(1.00, self.max_distance), 2)),
                self.newick)
    
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
    argument_parser = argparse.ArgumentParser(
        description="Data Collection for Extracting phylogenies from images using AI"
    )
    ########## ARGUMENTS ##########
    argument_parser.add_argument(
        "-r",
        "--create_rand_tree", 
        required=False, 
        action="store_true", 
        default=False,
        help="""On/Off flag. Randomizes the following parameters and flags within reasonable ranges if they haven't
        been specified: 
        package, amount_taxa, randomize_distances, max_distance, dont_allow_multifurcations, branch_vertical_margin, 
        fontsize, linewidth, dont_display_lengths, align_taxa. If --create_rand_tree and e.g. package are specified
        then package won't be randomized. Quickly create a diverse dataset instead of one type of image by 
        specifying --create_rand_tree in combination with --number_directories. --topology_only removes distances 
        and taxa from the newick and the image of each directory created with --create_rand_tree. --taxa_only
        removes just the distances in newick and image."""
    )
    argument_parser.add_argument(
        "-x",
        "--taxa_only",
        required=False,
        action="store_true",
        default=False,
        help="""On/Off flag. If specified removes distances from newick and image of the created directory. This is 
        helpful for creating datasets used for an initial training step in which models are just trained to recognize
        the tree's topology and taxa but not yet the branch lenghts. Can be combined with --create_rand_tree. 
        Default: False."""
    )
    argument_parser.add_argument(
        "-t",
        "--topology_only",
        required=False,
        action="store_true",
        default=False,
        help="""On/Off flag. If specified removes distances and taxa from newick and image of the created directory. 
        This is helpful for creating datasets used for an initial training step in which models are just trained to recognize
        the tree's topology. Can be combined with --create_rand_tree. 
        Default: False."""
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
                                    help="""Choose the width of branches in pixels. Default: 1. Should not exceed 1 if
                                    used with package phylo unless branch lengths are not displayed or they don't need
                                    to be fully visible.""")
    argument_parser.add_argument("-p", "--package", required=False, choices=["phylo", "ete3"],
                                    help="""Specify which package is used in the creation of the image. 
                                    Choose between Bio.Phylo and ETE3 Toolkit. Default: ete3.""")
    argument_parser.add_argument('-a', '--amount_taxa', required=False, type=int,
                                    help='Type=Int. Choose the preferred amount of generated taxa. Default: 10.')
    argument_parser.add_argument('-rd', '--randomize_distances', required=False, action="store_true", 
                                    help="""On/Off flag. If specified the distances (branch lengths) will be 
                                    randomized. Default: False.""")
    argument_parser.add_argument('-md', '--max_distance', required=False, type=int, 
                                    help="""Type=Int. If distances are randomized this will be the maximum distance. 
                                    Default: 1.""")
    argument_parser.add_argument("-db", "--dont_display_lengths", required=False, action="store_true", 
                                    help="""On/Off flag. If specified then no branch lengths will be displayed in the 
                                    image. Distances are still going to be present in the newick. 
                                    Default: False.""")
    argument_parser.add_argument("-c", "--circular_tree", required=False, action="store_true", default=False,
                                    help="""On/Off flag. ETE3 only. If True the tree in the image will be in 
                                    circular format. Default: False.""")
    argument_parser.add_argument("-rl", "--right_to_left_orientation", required=False, action="store_true", 
                                    default=False,help="""On/Off flag. ETE only. Specify if tree is oriented from 
                                    right to left (taxa on the left). Default: False (left to right orientation)""")
    argument_parser.add_argument("-da", "--dont_allow_multifurcations", required=False, action="store_true", 
                                    help="""On/Off flag. ETE3 only. If specified then the resulting 
                                    tree will be binary. Default: False.""")
    argument_parser.add_argument("-vm", "--branch_vertical_margin", required=False, type=int, 
                                    help="""ETE3 only. Amount of pixels between two adjacent branches. 
                                    Should not be smaller than 5. Default: 10 pixels.""")
    argument_parser.add_argument("-at", "--align_taxa", required=False, action="store_true", 
                                    help="""On/Off flag. ETE3 only. If specified taxa are aligned to the right instead 
                                    of written at the leaf node. Default: False.""")
    argument_parser.add_argument("--quiet", required=False, action="store_true", default=False,
                                    help="""On/Off flag. If --quiet is specified all console logs will be disabled and 
                                    only the output printed out. This is useful inside a pipeline where the newick is 
                                    piped into another application.""")
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
    create_rand_tree = args.create_rand_tree
    fontsize = args.fontsize
    linewidth = args.linewidth
    taxa_only = args.taxa_only
    topology_only = args.topology_only
    align_taxa = args.align_taxa
    quiet = args.quiet
    ########## CONFIGURE LOGGER ##########
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(logFormatter)
    console_logger.addHandler(stream_handler)
    if quiet:
        console_logger.setLevel(logging.WARNING)
    else:
        console_logger.setLevel(logging.INFO)
    print("\n\n\tData generation for extracting phylogenies from images using AI.\n\n")
    console_logger.info('Started')
    ########## SET DEFAULT VALUES ##########
    # remember if parameters relevant for randomize_treerender() were specified, those wont be randomized
    used_parameters = []
    # set default package
    if not package:
        package = "ete3"
    else:
        # package was specified, dont randomize the parameter later
        used_parameters.append("package")
    # set default fontsize and linewidth
    if not fontsize:
        if package == "ete3":
            fontsize = 8
        elif package == "phylo":
            fontsize = 16
    else:
        used_parameters.append("fontsize")
    if not linewidth:
        if package == "ete3":
            linewidth = 1
        elif package == "phylo":
            linewidth = 1
    else:
        used_parameters.append("linewidth")
    # set default ete3 params
    if package == "ete3":
        if not branch_vertical_margin:
            branch_vertical_margin = 10
        else:
            used_parameters.append("branch_vertical_margin")
        if not align_taxa:
            align_taxa = False
        else:
            used_parameters.append("align_taxa")
    if not amount_taxa:
        amount_taxa = 10
    else: 
        used_parameters.append("amount_taxa")
    if not max_distance:
        max_distance = 1.0
    else: 
        used_parameters.append("max_distance")
    if not dont_display_lengths:
        dont_display_lengths = False
    else: 
        used_parameters.append("dont_display_lengths")
    if not dont_allow_multifurcations:
        dont_allow_multifurcations = False
    else: 
        used_parameters.append("dont_allow_multifurcations")
    if not randomize_distances:
        randomize_distances = False
    else: 
        used_parameters.append("randomize_distances")
    
    ########## CHECKS ##########
    # TODO: remove unnecessary parameters from treerenders if topo or taxa_only is specified (not fontsize !!!)
    # TODO: warn user if he uses params compatible with topo or taxa only => font, linewidth, branch lengths etc
    # TODO: show default vertical margin for ete3 
    # TODO: implement max_taxa_length argument, is_taxid_valid 
    # TODO: remove node faces from inner nodes and leafs in ete3?
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
        param_warning = "You are using parameters that are ete3 only:\n"
        if branch_vertical_margin:
            param_warning += "branch_vertical\n"
            branch_vertical_margin = None
        if dont_allow_multifurcations:
            param_warning += "dont_allow_multifurcations\n" 
            dont_allow_multifurcations = False 
        if circular_tree:
            param_warning += "circular_tree\n"
            circular_tree = False
        if right_to_left_orientation:
            param_warning += "right_to_left_orientation\n"
            right_to_left_orientation = False
        console_logger.warning(param_warning)
        print(f"Resetting ete3 parameters. Do you want to proceed?")
        ask_user_to_continue()
    if number_directories < 1:
        raise ValueError(f"--number_directories {number_directories} not valid. Has to be "
                            "positive integer >= 1.")
    # warn user if he wants to create more than one image
    if number_directories > 1:
        print(f"Warning: You are about to create {number_directories} directories. Are you sure you want to proceed?")
        ask_user_to_continue()      
    ########## MAIN LOOP ##########
    # execute module <number_directories> times
    for i in range(number_directories):
        # unique file ID is just the current time (hour, minute, second, microsecond)
        # could be shortened to second and microsecond maybe
        file_id = str(datetime.datetime.now().strftime(r"%H%M%S%f"))
        console_logger.info(f"File ID: {file_id}")
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
            taxa_only=taxa_only,
            topology_only=topology_only,
            align_taxa=align_taxa,
        )
        # if the user wants to generate a tree with randomized parameters
        if create_rand_tree:
            tree_render.randomize_treerender(used_parameters)
        ########## CREATE THE NEWICK ##########
        tree_render.newick = generate_newick(tree_render.amount_taxa)
        # randomize the distances if a randomize_distances is specified
        if tree_render.randomize_distances:
            tree_render.randomize_distances_func()
        # save directory to outfile path if it was specified
        # If not create the generated data directory with the data directory inside
        tree_render.create_output_directory()
        iteration_info =f"""Default file ID: {file_id}
        Parameters:
        Fontsize: {tree_render.fontsize}
        Linewidth: {tree_render.linewidth}
        Taxa only: {tree_render.taxa_only}
        Topology only: {tree_render.topology_only}
        Randomize distances: {tree_render.randomize_distances}
        {"Max distance: "+{tree_render.max_distance} if tree_render.randomize_distances else "Distances: all exactly 1"}
        Amount of taxa: {tree_render.amount_taxa}
        Circular tree: {tree_render.circular_tree}
        Used package: {tree_render.package}
        Orientation: {"left to right" if not tree_render.right_to_left_orientation else "right to left"}
        Don't allow multifurcations: {tree_render.dont_allow_multifurcations}
        Vertical margin for adjacent branches: {tree_render.branch_vertical_margin if tree_render.package == "ete3" else None}
        Branch lengths displayed: {not tree_render.dont_display_lengths}
        Taxa aligned: {tree_render.align_taxa}
        """
        console_logger.info(iteration_info)
        
        print(f"Newick {i+1}:")
        print(f"  {tree_render.newick}")
# execute the main method
if __name__ == "__main__":
    main()


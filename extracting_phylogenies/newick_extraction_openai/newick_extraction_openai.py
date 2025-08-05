import sys
import os
import base64
import argparse
import re
from openai import OpenAI
import datetime
from argparse import RawTextHelpFormatter
from ete3 import Tree
from extracting_phylogenies.utilities import utilities as ut
import logging

# create logger
console_logger = logging.getLogger(__name__)

# environment variable for api key safety
client = OpenAI(api_key = os.getenv("BA_API_KEY"))

def get_image_from_directory(dir_path):
    """
    Given a path to a directory returns the image path. Supports the following formats: png, jpg, jpeg. 
    Iterates through all files in the given directory, checks if there is an image inside the directory and if it finds 
    the path to this one image is returned. 
    Multiple images inside one directory is not supported right now. 

    Args:
        dir_path (str): path to the directory containing the image of a phylogenetic tree and the corresponding newick tree
    """
    image_path = ""
    if not os.path.exists(dir_path):
        exit(f"[{ut.get_time()}] Error in get_image_from_directory: Path does not exist.") 
    elif not os.path.isdir(dir_path):
        exit(f"[{ut.get_time()}] Error in get_image_from_directory: Path is not a directory.")
    else:
        for image in os.listdir(dir_path):
            if image.endswith(".png") or image.endswith(".jpg") or image.endswith(".jpeg"):
                image_path = dir_path + str(f"\\{image}")
                break
    if image_path:
        return image_path
    else: 
        exit("No image found in the specified directory.") 

# TODO: add to utilities class
def get_file_id(dir_path):
    """
    Given the directory storing a file with an ID, returns its file ID.
    The ID is a number appended to the data directory and the image and newick file.
    If not file ID is found an empty string is returned

    Args:
        file_path (str): path to the file whose ID is returned 

    Returns:
        file_id (str): string of number at the end of the image filename 
    """
    file_path = get_image_from_directory(dir_path)
    if not (match_obj := re.search(r"\d+(?=\..{2,4})", file_path)) == None:
        file_id = match_obj.group()
        return file_id
    else: 
        return ""
        
# Function to encode the image
def encode_image(image_path):
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode("utf-8")
    
class Newick_Extraction_Job:
    """
    An instance of the Extractor class hold all values needed for extracting the newick from an image of a phylogenetic
    tree and saving it to file.
    """
    def __init__(
        self,
        infile_path, # path to image or dir containing image
        model = "gpt-4.1",
        newick = None,
        topology = None,
        taxa_only = False,
        topo_only = False,
        outfile_path = None,
        file_id = None,
        approach = "extract_nwk",
        extract_taxa = False,
    ):
        self.infile_path = infile_path
        self.model = model
        self.newick = newick
        self.topology = topology
        self.taxa_only = taxa_only
        self.topo_only = topo_only
        self.outfile_path = outfile_path
        self.file_id = file_id
        self.approach = approach
        self.extract_taxa = extract_taxa 
        
    def write_extracted_newick_to_file(self):
        """
        Called on a Newick_Extraction_Job object saves its newick at the specified path or in the current working 
        directory if no outfile path was specified. If outfile_path points to a file the newick is written into the file 
        if the file doesnt exist it is created. If outfile_path points to a directory the newick is written into the 
        directory. 

        Args:
            outfile_path (str): path where newick is supposed to be saved at
        """
        if self.outfile_path:
            app = "nwk_" if self.approach == "extract_nwk" else "topo_" if self.approach == "extract_topo" else ""
            flag = "taxa_" if self.extract_taxa else ""
            # if outfile path points to dir save it in the dir, if it points to a file then create the file
            if os.path.isdir(self.outfile_path):
                # create predictions directory for the extracted newick if it doesnt already exist in the specified dir
                predictions_path = "predictions"
                parent_path = self.outfile_path
                newick_path = f"{self.model}_{app}{flag}{self.file_id}.nwk"
                os.makedirs(dir_path := os.path.join(parent_path, predictions_path), exist_ok=True)  
                # # if the file already exists the current file content will be overwritten
                with open(os.path.join(dir_path, newick_path), "w") as nwk_file:
                    nwk_file.write(self.newick)
            elif os.path.isfile(self.outfile_path):
                # if the file already exists the current file content will be overwritten
                with open(self.outfile_path, "w") as nwk_file:
                    nwk_file.write(self.newick)
            else: 
                raise FileNotFoundError(
                    f"[{ut.get_time()}] Error in write_extracted_newick_to_file: outfile_path {self.outfile_path} is "
                    "not valid."
                )
        # if no outfile path was given, save the newick in the current working directory
        else:
            if os.path.isfile(f".\\{self.model}_{app}{flag}{self.file_id}.nwk"):
                raise FileExistsError(f"[{ut.get_time()}] Error in write_newick_to_file: File already exists.") 
            else:
                with open(f".\\{self.model}_{app}{flag}{self.file_id}.nwk", "w") as nwk_file:
                    nwk_file.write(self.newick)
                    
    # TODO utilities class?
    def get_image_path(self):
        """
        Function for returning the image_path from a infile_path. If infile_path points to a dir the path to the first
        image inside is returned. If it points to a file then self.infile_path is returned. If the path doesnt exist an
        exception is raised.

        Returns:
            str: path to image in Newick_Extraction_Job object's infile_path
        """
        # if the infile path is a dir get the img from the dir
        if os.path.isdir(self.infile_path):
            return get_image_from_directory(self.infile_path)
        # if infile path is a file then get the img directly
        elif os.path.isfile(self.infile_path):
            return self.infile_path
        else:
            raise FileNotFoundError(
                f"[{ut.get_time()}] Error in get_image_path(): Infile path {self.infile_path} does not exist."
            )
            
    def extract_newick_directly(self):
        """
        Given a path to an image (png, jpg or jpeg) and a OpenAI model returns the newick.
        "Directly" refers to the fact that the AI is not given the taxa in a first step and tasked with generating
        the Newick format directly instead of letting the AI write out the basic hierarchical structure of the image and 
        then extracting the Newick from that manually.

        Args:
            image_path (str): path to the image
            model (str, optional): OpenAI model that is used. Defaults to "gpt-4.1".

        Returns:
            str: response text by the model
        """
        if not self.model in ["gpt-4.1", "o4-mini", "gpt-4o"]:
            raise ValueError(f"""[{ut.get_time()}] Error in generate_newick_from_image_directly: Given model is not 
                             supported. Please choose from gpt-4.1, o4-mini or gpt-4o""")
        img_path = self.get_image_path()
        prompt = """Give me the Newick format of this phylogenetic tree, preserving all taxa (if visible), all branch
        lengths (if visible) and topology."""
        instructions =  """You are given an image of a phylogenetic tree that may have multifurcations. You task is to 
        output only the tree in valid Newick format. In case there are no branch lengths then infere the branch 
        lengths from the scale bar.
        
        Guidelines:
        - Reply with nothing but the newick, no explanation, no prefix, no suffix
        - Dont put backticks around the newick 
        - dont but new lines in the newick
        - The Newick string must include the taxon names exactly as shown in the image
        - The Newick string must correspond to the topology exactly as shown in the image
        - The Newick string must include branch lengths exactly as shown in the image
        - In case there is a scale bar but no branch lengths in the whole image then infere the branch lengths from the 
        scale bar
        - In case there are no branch lengths in the whole image and no scale bar then dont include branch lengths in 
        the Newick string
        - In case there are no branch lengths, no scale bar and no taxa in the whole tree then dont include branch 
        lengths in the newick and put empty strings for taxon names e.g. ((,),(,(,)));
        
        Examples:
        
        Example with taxa and branch lengths: 
        ((<taxon1>:<branch_length1>,<taxon2>:<branch_length2>):<branch_length4>,<taxon3>:<branch_length3>);
        ((A:2.37,B:1.55):4.58,((C:1.43,D:3.63):0.27,E:1.66):4.07);
        Example without branch lengths.
        ((<taxon1>,<taxon2>),<taxon3>);
        ((A,B),((C,D),E));
        Example without branch lengths and taxa:
        ((,),((,),));
        Example with multifurcations:
        ((A:4.24,B:2.21,C:9.11):3.31,(D:2.22,E:1.02):1.21,((F:4.26,G:6.66):5.01,H:1.32):3.21);
        ((A,B,C),(D,E),((F,G),H));
        ((,,),(,),((,),));
        """
        b64_image = encode_image(img_path)
        # conditionally include temperature and top_p in the create() functions args because o4-mini does not support it
        args = {
            "model":self.model,
            "instructions":instructions,
            "input":[
                {
                    "role": "user", 
                    "content": [
                        { "type": "input_text", "text": prompt},
                        { "type": "input_image", "image_url": 
                            f"data:image/{ut.get_image_format(img_path)};base64,{b64_image}"},
                    ],
                }
            ],
        }
        if self.model == "gpt-4.1" or self.model == "gpt-4o":
            args["temperature"] = 0.0
            args["top_p"] = 0.0
        response = client.responses.create(**args)
        console_logger.info(f"Model response: {response.output_text}")
        # parse the newick from the reponse
        if re.search(r"\(.+;", response.output_text):
            newick = re.search(r"\(.+;", response.output_text).group()
        else:
            raise ValueError(f"No Newick found in model response.")
        console_logger.info(f"Extracted newick: {newick}")
        # set taxa_only and topo_only to true if topology has no distances or no distances and no taxon names
        if not (branch_lengths := re.search(r":\d+(\.\d+){0,1}", newick)) and not \
            re.search(r"(?<=[,(])[\d\w\.\-\"\'#]+(?=[\:\)\,])", newick):
                self.topo_only = True
        elif not branch_lengths:
            self.taxa_only = True
        return newick

    def extract_newick_from_topology(self):
        """
        Used on a Newick_Extraction_Job object with a topology similar to that of Bio.Phylo created by a model  returns 
        the newick of the topology

        Args:
            topology (str): string in hierarchical text format 

        Returns:
            str: newick corresponding to the hierarchical text format
        """
        def get_name(line):
            return match.group() if (match := re.search(r"(?<=name=[\'\"]).+(?=[\'\"]\))", line)) else None
        def get_dist(line):
            return match.group() if (match := re.search(r"(?<=branch_length=)\d+(\.\d+){0,1}", line)) else None
        def get_indentation_level(line):
            # one tab is 4 spaces
            return (len(line) - len(line.lstrip())) / 4
        # get each line of the topology string
        lines = self.topology.split("\n")
        # create ete3 tree object with root
        tree = Tree()
        current_parent = tree
        # remember the order of the internal nodes with a stack
        internal_nodes_stack = [tree]
        # loop starts with first child of the root
        for i in range(1, len(lines)):
            current_indentation = get_indentation_level(lines[i]) 
            next_indentation = get_indentation_level(lines[i+1]) if i+1 < len(lines) else -1 
            if i+1 < len(lines) and next_indentation > current_indentation:
                # if next line is indented more then the current line is an internal node
                internal_node = current_parent.add_child(dist=get_dist(lines[i]))
                internal_nodes_stack.append(internal_node)
                # update the current node to add children to current line
                current_parent = internal_node
            elif i+1 < len(lines) and next_indentation == current_indentation:
                # if next line has the same indentation then the current line is a leaf 
                current_parent.add_child(name=get_name(lines[i]), dist=get_dist(lines[i]))
            elif i+1 < len(lines) and next_indentation < current_indentation:
                # if next line has less indentation then the current line is a leaf and the next line belongs to another 
                # parent node and the current one is finished
                current_parent.add_child(name=get_name(lines[i]), dist=get_dist(lines[i]))
                # pop all internal nodes that are now finished
                for j in range(int(current_indentation - next_indentation)):
                    internal_nodes_stack.pop()
                # assign current_parent the new topmost node of the stack
                current_parent = internal_nodes_stack[len(internal_nodes_stack)-1]
            else:
                # there is no next line, add leaf to current internal node
                current_parent.add_child(name=get_name(lines[i]), dist=get_dist(lines[i]))
        # remove support vals to not confuse postprocessing
        newick = ut.remove_support_vals(tree.write())
        # remove placeholder distances from newick
        newick = ut.remove_distances(newick)
        console_logger.info(f"Newick extracted from topology: {newick}")
        return newick

    def extract_topology(self):
        """
        Given a path to an image (png, jpg or jpeg) and a OpenAI model returns the phylogenetic tree in hierarchical 
        text format for manual extraction. 
        Example with taxa and branch lengths:
        Clade()
            Clade(branch_length=3.54)
                Clade(branch_length=3.42)
                    Clade(branch_length=4.88, name='Crassulaceae')
                    Clade(branch_length=3.53, name='Calycanthus_chinensis')
                Clade(branch_length=1.8, name='Verruciconidia_persicina')
            Clade(branch_length=3.27)
                Clade(branch_length=0.42, name='Aquaspirillum_serpens')
                Clade(branch_length=1.74, name='Wolbachia_pipientis')
                
        for:
                                        _____________________ Crassulaceae
                         ______________|
          ______________|              |_______________ Calycanthus_chinensis
         |              |
        _|              |_______ Verruciconidia_persicina
         |
         |              _ Aquaspirillum_serpens
         |_____________|
                       |_______ Wolbachia_pipientis
        
        Args:
            image_path (str): path to the image
            model (str, optional): OpenAI model that is used. Defaults to "gpt-4.1".

        Returns:
            str: response text by the model
        """
        if not self.model in ["gpt-4.1", "o4-mini", "gpt-4o"]:
            raise ValueError(f"""Error in generate_newick_from_image_directly: Given model is not supported. Please 
                             choose from gpt-4.1, o4-mini or gpt-4o""")
        img_path = self.get_image_path()
        prompt = """Give me the topology, taxon names (if visible) and branch lengths (if visible) of this phylogenetic 
        tree in hierarchical text format as if you created a Bio.Phylo Tree object and printed it."""
        instructions =  """You are given an image of a phylogenetic tree. Your task is to output the topology, taxon 
        names (if visible) and branch lengths (if visible) in a hierarchical text format similar to that of Bio.Phylo
        when a Tree object is printed.
        
        Example with taxon names and branch lengths: 
        Clade()
            Clade(branch_length=3.54)
                Clade(branch_length=3.42)
                    Clade(branch_length=4.88, name='Crassulaceae')
                    Clade(branch_length=3.53, name='Calycanthus_chinensis')
                Clade(branch_length=1.8, name='Verruciconidia_persicina')
            Clade(branch_length=3.27)
                Clade(branch_length=0.42, name='Aquaspirillum_serpens')
                Clade(branch_length=1.74, name='Wolbachia_pipientis')
                
        this corresponds to a tree like:
        
                                        _____________________ Crassulaceae
                         ______________|
          ______________|              |_______________ Calycanthus_chinensis
         |              |
        _|              |_______ Verruciconidia_persicina
         |
         |              _ Aquaspirillum_serpens
         |_____________|
                       |_______ Wolbachia_pipientis
        
        Guidelines:
        - Reply with nothing but the topology, no explanation, no prefix, no suffix 
        - Each indentation is 4 spaces
        - Each indentation level corresponds to one clade deeper in the tree's hierarchy
        - No indentation means it is the root
        - One tab of indentation means that the clade or leaf is the child of the root
        - Two tabs mean that the clade or leaf is the grandchild of the root and so on
        - The topology string must include all taxon names and branch lengths
        - In case there are no taxon names in the image then put empty strings in their place e.g. ..., name='')
        - In case there are no branch lengths in the whole image then infere the branch lengths from the scale bar
        - In case there are no branch lengths, no scale bar in the whole image then set all branch lengths to empty 
        string e.g. Clade(branch_lengths=), Clade(branch_lengths=, name='taxon2')
        """
        b64_image = encode_image(img_path)
        # conditionally include temperature and top_p in the create() functions args because o4-mini does not support it
        args = {
            "model":self.model,
            "instructions":instructions,
            "input":[
                {
                    "role": "user", 
                    "content": [
                        { "type": "input_text", "text": prompt},
                        { "type": "input_image", "image_url": 
                            f"data:image/{ut.get_image_format(img_path)};base64,{b64_image}"},
                    ],
                }
            ],
        }
        if self.model == "gpt-4.1" or self.model == "gpt-4o":
            args["temperature"] = 0.0
            args["top_p"] = 0.0
        response = client.responses.create(**args)
        topology = response.output_text
        # set taxa_only and topo_only to true if topology has no distances or no distances and no taxon names
        if not (branch_lengths := re.search(r"(?<=branch_length=)\d", topology)) and not \
            re.search(r"(?<=name=)[\']{0,1}\w", topology):
                self.topo_only = True
        elif not branch_lengths:
            self.taxa_only = True
        console_logger.info(f"Extracted topology: {topology}")
        console_logger.info(f"Topology is taxa-only: {self.taxa_only}")
        console_logger.info(f"Topology is topo-only: {self.topo_only}")
        return topology
    
    def postprocess_newick(self, newick):
        """
        Given a newick returns the postprocessed one. 
        Passes the newick to the Job object's model if necessary and does manual postprocessing if it fails.
        Returns empty tree (";") if the postprocessing fails i.e. the newick isn't parseable by the ete3 module.

        Args:
            newick (str): newick string that needs postprocessing

        Returns:
            str: postprocessed newick
        """
        # set default format 
        format = 0
        # set format for topology only newicks
        if self.topo_only:
            console_logger.info("Newick is topology-only.")
            self.topo_only = True
            format = 100
        # if newick's formatting is correct skip AI postprocessing otherwise do it
        if ut.is_newick(newick, format=format):
            console_logger.info(f"Skipping AI post-processing.")    
            # create tree and write it to let ete3 remove unnecessary parentheses and spaces
            newick = Tree(newick, format=format).write()
            # remove placeholder distances by ete3 if newick originally didnt have any
            if self.taxa_only or self.topo_only:
                newick = ut.remove_distances(newick)
            # remove support values
            newick = ut.remove_support_vals(newick)
            console_logger.info(f"Post-processing was successful. Newick: {newick}")
            return newick
        else:
            console_logger.info(f"Starting AI-postprocessing.")
            newick = self.correct_newick(newick)
            if not ut.is_newick(newick):
                console_logger.info(f"AI post-processing failed. Continuing with manual postprocessing.")
            # remove special chars from taxa by removing them from each taxon seperately
            newick = re.sub(
                r"(?<=[,(])[\d\w\.\-\"\'#]+(?=[\:\)\,])", # matches taxa
                lambda t: ut.remove_special_chars(t.group()), 
                newick)
            # balance parentheses if necessary
            if not ut.is_balanced(newick):
                console_logger.info("Balancing out parentheses.")
                newick = ut.balance_parentheses(newick)
            # create tree and write it to let ete3 remove unnecessary parentheses and spaces
            if ut.is_newick(newick, format=format):
                newick = Tree(newick, format=format).write()
                # remove placeholder distances by ete3 if newick originally didnt have any
                if self.taxa_only or self.topo_only:
                    newick = ut.remove_distances(newick)
                # remove support values
                newick = ut.remove_support_vals(newick)
                console_logger.info(f"Post-processing was successful. Newick: {newick}")
                return newick
            # if formatting is still wrong, warn user and return empty tre
            else:
                console_logger.warning(f"Post-processing failed. Returning empty tree.")
                return ";"
        
    def correct_newick(self, erroneous_newick):
        """
        Given a newick uses the model to correct the given erroneous newick and returns the corrected newick
        
        Args
        """
        b64_image = encode_image(self.get_image_path())
        prompt=f"""This is a phylogenetic tree and this the corresponding erroneous Newick: {erroneous_newick}. The 
        Newick might have spelling mistakes, missing taxa and wrong formatting. Give me the correct Newick."""
        instructions = """You are given an image of a phylogenetic tree that may have multifurcations and a string in 
        Newick format that may have false formatting and spelling mistakes. Your task is to fix errors in the string in 
        Newick format and output a string in Newick format with valid formatting and no spelling mistakes.  
        
        Example of correct Newick format: 
        ((<taxon1>:<branch_length1>,<taxon2>:<branch_length2>):<branch_length4>,<taxon3>:<branch_length3>);
        
        Guidelines:
        - Reply with nothing but the newick, no explanation, no prefix, no suffix
        - Dont put backticks around the newick and dont but new lines in the newick
        - Correct parentheses with respect to the topology of the phylogenetic tree in the given image 
        - Make sure every opening parenthesis has a closing parenthesis 
        - Make sure there aren't any unnecessary parentheses
        - If the newick doesn't contain branch lengths and/or taxa, this is not part of the error 
        """
        # conditionally include temperature and top_p in the create() functions args because o4-mini does not support it
        args = {
            "model":self.model,
            "instructions":instructions,
            "input":[
                {
                    "role": "user", 
                    "content": [
                        { "type": "input_text", "text": prompt},
                        { "type": "input_image", "image_url": 
                            f"data:image/{ut.get_image_format(self.get_image_path())};base64,{b64_image}"},
                    ],
                }
            ],
        }
        if self.model == "gpt-4.1" or self.model == "gpt-4o":
            args["temperature"] = 0.0
            args["top_p"] = 0.0
        response = client.responses.create(**args)
        if re.search(r"\(.+;", response.output_text):
            updated_newick = re.search(r"\(.+;", response.output_text).group()
        else:
            console_logger.warning(f"No Newick found in the model's response. Returning empty newick.")
            return ";"
        return updated_newick

def main():
    argument_parser = argparse.ArgumentParser(
        description="""Newick encoding of images containing phylogenetic trees using AI.\n  Modes:
        generate_newick_directly: \n
          You can either provide a single directory with an image or the path to your image 
          + "the path where the image is saved at. If you provide just the directory then the first image found in the
          directory will be used and the resulting newick saved into this directory. 'directly' refers to the AI being
          told to generate the Newick format of an image directly without an intermediate manual processing step which
          tests the AIs Newick formatting ablities.""",
    )
    # argument for passing the path where the newick/image pair is saved at
    argument_parser.add_argument(
        "-a",
        "--approach", 
        required=True, 
        choices=["extract_nwk", "extract_topo"],
        help="""Specify which functionality you want to use. extract_nwk extracts the newick directly i.e the model 
        produces text in newick format itself. extract_topo extracts the topology of the image and then the newick
        is extracted algorithmically.""")
    argument_parser.add_argument('-i', '--infile_path', required=True, type=str,
                                 help="""Path to a directory containing an image of a phylogenetic tree or path to the 
                                 image itself. Be aware that the first image in a specified directory will be 
                                 chosen.""")
    argument_parser.add_argument('-o', '--outfile_path', required=False, type=str, 
                                 help="""Path where the newick is saved at. If the path points to a directory the 
                                 newick will be saved there inside a predictions directory. If no path is provided the 
                                 output is printed to console to allow for use in pipelines.""")
    argument_parser.add_argument("-m", "--model", required=False, choices=["gpt-4o", "gpt-4.1", "o4-mini"], type=str, 
                                 default="gpt-4.1",
                                 help="""Choose which OpenAI model is used in the generation of the newick string. 
                                 Choose from GPT-4o, GPT-4.1 or o4-mini. Default model used is GPT-4.1.""")
    # argument_parser.add_argument("-t", "--extract_taxa", required=False, action="store_true", default=False,
    #                              help="""On/Off flag. If specified taxa will be extracted by the chosen model in a 
    #                              first step and then passed to the model together with the image.""")
    argument_parser.add_argument("--quiet", required=False, action="store_true", default=False,
                                 help="""On/Off flag. If --quiet is specified all console logs will be disabled and 
                                 only the output printed out. This is useful inside a pipeline where the newick is 
                                 piped into another application.""")
    
    # Specified parameters
    args = argument_parser.parse_args()
    infile_path = args.infile_path
    outfile_path = args.outfile_path
    model = args.model 
    approach = args.approach
    quiet = args.quiet
    # extract_taxa = args.extract_taxa
    
    # file ID for the current job 
    file_id = ut.get_file_id()
    
    # configure logger
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(logFormatter)
    console_logger.addHandler(stream_handler)
    if quiet:
        console_logger.setLevel(logging.WARNING)
    else:
        console_logger.setLevel(logging.INFO)
    console_logger.info('Started')
    console_logger.info(f"File ID: {file_id}")
    
    # create Newick_Extraction_Job object
    extraction_job = Newick_Extraction_Job(
        infile_path=infile_path,
        model=model,
        outfile_path=outfile_path,
        file_id = file_id,
        approach=approach,
    )
        
    if approach == "extract_nwk":
        # extract newick
        newick = extraction_job.extract_newick_directly()
        # set the postprocessed newick as the job objects newick
        extraction_job.newick = extraction_job.postprocess_newick(newick)
        if extraction_job.outfile_path:
            extraction_job.write_extracted_newick_to_file()
        else:
            print(extraction_job.newick) 
    elif approach == "extract_topo":
        # extract topology
        extraction_job.topology = extraction_job.extract_topology() 
        # extract newick from topology
        newick = extraction_job.extract_newick_from_topology()
        # set the postprocessed newick as the job objects newick
        extraction_job.newick = extraction_job.postprocess_newick(newick)
        # if not outfile path is specified then print the newick to console for use in pipelines
        if extraction_job.outfile_path:
            extraction_job.write_extracted_newick_to_file()
        else:
            print(extraction_job.newick)
    else:
        raise ValueError(f"[{ut.get_time()}] Approach {approach} not valid.")
    
    console_logger.info('Finished')
            
# execute the main method
if __name__ == "__main__":
    main()
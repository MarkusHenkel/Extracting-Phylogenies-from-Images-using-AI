import os
import base64
import argparse
import re
from openai import OpenAI
import datetime
from argparse import RawTextHelpFormatter
from ete3 import Tree
import utilities as ut

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
    Given the directory storing the image to be encoded, returns its file ID, a number appended to the data directory 
    and the image and newick file 

    Args:
        image_path (str): Path to the image to be encoded 

    Returns:
        file_id (str): string of number at the end of the image filename 
    """
    image_path = get_image_from_directory(dir_path)
    if not (match_obj := re.search(r"\d+(?=\..{2,4})", image_path)) == None:
        file_id = match_obj.group()
        return file_id
    else: 
        return ""
        
# TODO: add to utilities class        
def get_image_format(image_path):
    """
    Using regex matches 2-4 chars at the end of of the path and returns the matched string if it is either jpeg, 
    png or jpg

    Args:
        image_path (str): path of the file

    Returns:
        str: image file format
    """
    if not re.search(r"(?<=\.)\w{3,4}$", image_path):
        exit(f"[{ut.get_time()}] Error in get_image_format: Given path doesn't have a file ending.")
    elif (file_ending := re.search(r"(?<=\.)\w{3,4}$", image_path).group()) in ["jpeg", "png", "jpg"]:
        return file_ending
    else:
        exit(f"[{ut.get_time()}] Error in get_image_format: Given file path does not end with .jpeg, .png or .jpg.")
        
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
        outfile_path = None,
        file_id = None,
        approach = "extract_nwk",
        extract_taxa = False,
    ):
        self.infile_path = infile_path
        self.model = model
        self.newick = newick
        self.topology = topology
        self.outfile_path = outfile_path
        self.file_id = file_id
        self.approach = approach
        self.extract_taxa = extract_taxa 
        
    def write_newick_to_file(self):
        """
        Called on a Newick_Extraction_Job object saves its newick at the specified path or in the current working 
        directory if no outfile path was specified. If outfile_path points to a file the newick is written into the file 
        if the file doesnt exist it is created. If outfile_path points to a directory the newick is written into the 
        directory. 

        Args:
            outfile_path (str): path where newick is supposed to be saved at
        """
        if self.outfile_path:
            app = "nwk" if self.approach == "extract_nwk" else "topo" if self.approach == "extract_topo" else None 
                
            # if outfile path points to dir save it in the dir, if it points to a file then create the file
            if os.path.isdir(self.outfile_path):
                newick_path = self.outfile_path + f"\\{self.model}_{app}_{self.file_id}.nwk"
                with open(newick_path, "w") as nwk_file:
                    nwk_file.write(self.newick)
            elif os.path.isfile(self.outfile_path):
                with open(self.outfile_path, "w") as nwk_file:
                    nwk_file.write(self.newick)
            else: 
                exit(f"[{ut.get_time()}] Error in write_newick_to_file: outfile_path {self.outfile_path} is not valid.")
        # if no outfile path was given, save the newick in the current working directory
        else:
            if os.path.isfile(f".\\{self.model}_{app}_{self.file_id}.nwk"):
                raise FileExistsError(f"[{ut.get_time()}] Error in write_newick_to_file: File already exists.") 
            else:
                with open(f".\\{self.model}_{app}_{self.file_id}.nwk", "w") as nwk_file:
                    nwk_file.write(self.newick)
                    
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
            raise FileNotFoundError(f"[{ut.get_time()}] Error in get_image_path(): Infile path {self.infile_path} does \
                not exist.")
            
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
        prompt = """Give me the newick of this phylogenetic tree with correct taxa and topology."""
        instructions =  """You are given images of phylogenetic trees and are tasked to respond with corresponding Newick 
        and just the Newick. 
        Respond with the corresponding Newick in correct format: 
        ((<taxon1>:<distance1>,<taxon2>:<distance2>):<clade_distance1>,<taxon3>=<distance3>); e.g.
        ((Sipunculidae:2.37,Tenebrio_molitor:1.55):4.58,((Siren:1.43,Catreus:3.63):0.27,Notothenioidei:1.66):4.07);.
        Make sure the Newick string includes the correct taxon names, all distances 
        (if scalebar or branch lenghts are present) and most importantly the correct tree topology.
        Make sure every opening parenthesis has a closing parenthesis. Make sure no parentheses pair is superfluous."""
        b64_image = encode_image(img_path)
        response = client.responses.create(
            # make model more deterministic
            temperature=0,
            model=self.model,
            instructions=instructions,
            input=[
                {
                    "role": "user", 
                    "content": [
                        { "type": "input_text", "text": prompt},
                        { "type": "input_image", "image_url": 
                            f"data:image/{get_image_format(img_path)};base64,{b64_image}"},
                    ],
                }
            ],
        )
        # parse the newick from the reponse
        if re.search(r"\(.+;", response.output_text):
            newick = re.search(r"\(.+;", response.output_text).group()
            print(f"[{ut.get_time()}] Generated Newick before postprocessing using AI:\n{newick}")
        else:
            exit(f"[{ut.get_time()}] Error in generate_newick_format(): No Newick found in the response.")
        # postprocess if newick formatting is wrong TODO: get results on if AI is capable of correcting the newick
        if not ut.is_newick(newick):
            print(f"[{ut.get_time()}] generate_newick_from_image_directly: Formatting is wrong. Beginning postprocessing.")
            newick = correct_newick(img_path, self.model)
            # check if formatting is still wrong after postprocessing TODO: do manual postprocessing
            if not ut.is_newick(newick):
                print(f"[{ut.get_time()}] generate_newick_from_image_directly: AI-Postprocessing failed.") 
        else:
            # TODO replace with logging
            print(f"""[{ut.get_time()}] generate_newick_from_image_directly: Formatting is correct. Skipping 
                  AI-postprocessing.""")
        print(f"[{ut.get_time()}] Generated Newick:\n{newick}")
        return newick

    def extract_newick_from_topology(topology):
        if not topology.split("\n")[0] == "Clade()":
            exit(f"[] Error in extract_newick_from_topology: Topology does not start with 'Clade()'")
        def get_name(line):
            return match.group() if (match := re.search(r"(?<=name=[\'\"]).+(?=[\'\"]\))", line)) else None
        def get_dist(line):
            return match.group() if (match := re.search(r"(?<=branch_length=)\d+(\.\d+){0,1}", line)) else None
        def get_indentation_level(line):
            # one tab is 4 spaces
            return (len(line) - len(line.lstrip())) / 4
        # get each line of the topology string
        lines = topology.split("\n")
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
        return ut.remove_support_vals(tree.write())

    def extract_topology(self):
        """
        Given a path to an image (png, jpg or jpeg) and a OpenAI model returns the phylogenetic tree in hierarchical text
        format for manual extraction. 
        E.g.
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
         |             |
        _|             |_______ Verruciconidia_persicina
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
            raise ValueError(f"""[{ut.get_time()}] Error in generate_newick_from_image_directly: Given model is not 
                             supported. Please choose from gpt-4.1, o4-mini or gpt-4o""")
        img_path = self.get_image_path()
        prompt = """Describe this phylogenetic tree and its topology with hierarchical text that includes all taxa and
        distances."""
        instructions =  """You encode images of phylogenetic trees into this hierarchical text format similar to that of 
        Bio.Phylo:
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
                    
        Each tab of indentation means one level further down in hierarchy of the phylogenetic tree. No indentation means 
        it is the root. One space of indentation means that the clade or leaf is the child of the root. Two tabs mean 
        that the clade or leaf is the grandchild of the root and so on. 
        Make sure all distances are included, that correct formatting is correct and that the taxon names are correct."""
        b64_image = encode_image(img_path)
        response = client.responses.create(
            temperature=0,
            model=self.model,
            instructions=instructions,
            input=[
                {
                    "role": "user", 
                    "content": [
                        { "type": "input_text", "text": prompt},
                        { "type": "input_image", "image_url": 
                            f"data:image/{get_image_format(img_path)};base64,{b64_image}"},
                    ],
                }
            ],
        )
        topology = response.output_text
        return topology
    
def correct_newick(image_path, generated_newick, model="gpt-4.1"):
    """
    Given a newick with incorrect formatting and the original image, does a postprocessing step with focus on correcting 
    spelling mistakes and mistakes in the newick formatting using the API once again

    Args:
        image_path (str): path to the orignal image
        generated_newick (str): generated newick that needs postprocessing

    Returns:
        updated_newick: newick with improved spelling and formatting 
    """
    b64_image = encode_image(image_path)
    prompt=f"""This is a phylogenetic tree and this the corresponding Newick: {generated_newick}. The Newick might have 
    spelling mistakes, missing taxa and wrong formatting. Give me the correct Newick."""
    instructions = """You encode images of phylogenetic trees into Newick format. Respond only with the newick.
    Taxa are from the NCBI taxonomy. Correct Newick format: 
    ((<taxon1>:<distance1>,<taxon2>:<distance2>):<clade_distance1>,<taxon3>=<distance3>);. 
    Make sure every opening parenthesis has a closing parenthesis. Make sure there aren't any unnecessary parentheses, 
    every parenthesis pair where the closing parentheses isn't followed by a distance i.e. ':0.1' or a semicolon is 
    unnecessary."""
    response = client.responses.create(
        model=model,
        instructions=instructions,
        input=[
            {
                "role": "user", 
                "content": [
                    { "type": "input_text", "text": prompt},
                    { "type": "input_image", "image_url": 
                        f"data:image/{get_image_format(image_path)};base64,{b64_image}"},
                ],
            }
        ],
    )
    if re.search(r"\(.+;", response.output_text):
        updated_newick = re.search(r"\(.+;", response.output_text).group()
    else:
        exit(f"[{ut.get_time()}] Error in generate_newick_format(): No Newick found in the response.")
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
        # formatter_class=RawTextHelpFormatter
    )
    #
    # Arguments
    #
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
                                 newick will be saved there. If no path is provided the newick is saved in the 
                                 current working directory.""")
    argument_parser.add_argument("-m", "--model", required=False, choices=["gpt-4o", "gpt-4.1", "o4-mini"], type=str, 
                                 default="gpt-4.1",
                                 help="""Choose which OpenAI model is used in the generation of the newick string. 
                                 Choose from GPT-4o, GPT-4.1 or o4-mini. Default model used is GPT-4.1.""")
    argument_parser.add_argument("-t", "--extract_taxa", required=False, action="store_true", default=False,
                                 help="""On/Off flag. If specified taxa will be extracted by the chosen model in a 
                                 first step and then passed to the model again """)
    
    # Specified parameters
    args = argument_parser.parse_args()
    infile_path = args.infile_path
    outfile_path = args.outfile_path
    model = args.model 
    approach = args.approach
    extract_taxa = args.extract_taxa
    
    # file ID for the current job 
    file_id = ut.get_file_id()
    
    # create Newick_Extraction_Job object
    extraction_job = Newick_Extraction_Job(
        infile_path=infile_path,
        model=model,
        outfile_path=outfile_path,
        file_id = file_id,
        approach=approach,
    )
        
    # Mode: generation of newick
    if approach == "extract_newick_directly":
        extraction_job.newick = extraction_job.extract_newick_directly()
        print(f"[{ut.get_time()}] Extracted Newick:")
        print(extraction_job.newick)
        extraction_job.write_newick_to_file()
    elif approach == "extract_topology":
        extraction_job.topology = extraction_job.extract_topology()
        print(f"[{ut.get_time()}] Extracted topology:")
        print(extraction_job.topology)
        extraction_job.newick = extraction_job.extract_newick_from_topology()
        print(f"[{ut.get_time()}] Extracted Newick:")
        print(extraction_job.newick)
        extraction_job.write_newick_to_file()
    else:
        raise ValueError(f"[{ut.get_time()}] Approach {approach} not valid.")
            
# execute the main method
main()

# path to test directory for testing
# Random distances (0.00 and 20.00), 20 taxa: 
# "C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\data-collection\20_taxa_rand_distances"
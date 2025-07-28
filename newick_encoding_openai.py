from data_collection import phylogenetic_trees_and_newicks as ptan # writing the phylogenetic tree of the generated newick
import os
import base64
import argparse
import re
from openai import OpenAI
import datetime
from argparse import RawTextHelpFormatter
from ete3 import Tree

# environment variable for api key safety
client = OpenAI(api_key = os.getenv("BA_API_KEY"))

# TODO: add to utilities class
# time for logs 
def get_time():
    return datetime.datetime.now().strftime("%Y-%b-%d %H:%M:%S")

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
        exit(f"[{get_time()}] Error in get_image_from_directory: Path does not exist.") 
    elif not os.path.isdir(dir_path):
        exit(f"[{get_time()}] Error in get_image_from_directory: Path is not a directory.")
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
        
# Function to encode the image
def encode_image(image_path):
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode("utf-8")

# TODO: add to utilities class
def write_newick_to_file(newick, model, file_id="", outfile_path=None):
    """
    Given a newick saves it at the specified path or in the current working directory if not outfile path was specified. 
    If outfile_path points to a file the newick is written into the file if the file doesnt exist it is created. If 
    outfile_path points to a directory the newick is written into the directory. 

    Args:
        outfile_path (str): path where newick is supposed to be saved at
    """
    if outfile_path:
        # if outfile path points to dir save it in the dir, if it points to a file then create the file
        if os.path.isdir(outfile_path):
            newick_path = outfile_path + f"\\{model}_{file_id}.nwk"
            with open(newick_path, "w") as nwk_file:
                nwk_file.write(newick)
        elif os.path.isfile(outfile_path):
            with open(outfile_path, "w") as nwk_file:
                nwk_file.write(newick)
        else: 
            exit(f"[{get_time()}] Error in write_newick_to_file: outfile_path {outfile_path} is not valid.")
    # if no outfile path was given, save the newick in the current working directory
    else:
        if os.path.isfile(f".\\{model}_{file_id}.nwk"):
            exit(f"[{get_time()}] Error in write_newick_to_file: File already exists.") 
        else:
            with open(f".\\{model}_{file_id}.nwk", "w") as nwk_file:
                nwk_file.write(newick)

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
        exit(f"[{get_time()}] Error in get_image_format: Given path doesn't have a file ending.")
    elif (file_ending := re.search(r"(?<=\.)\w{3,4}$", image_path).group()) in ["jpeg", "png", "jpg"]:
        return file_ending
    else:
        exit(f"[{get_time()}] Error in get_image_format: Given file path does not end with .jpeg, .png or .jpg.")
        
def extract_newick_directly(image_path, model="gpt-4.1"):
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
    if not model in ["gpt-4.1", "o4-mini", "gpt-4o"]:
        exit(f"""[{get_time()}] Error in generate_newick_from_image_directly: Given model is not supported. 
             Please choose from gpt-4.1, o4-mini or gpt-4o""")
    prompt = """Give me the newick of this phylogenetic tree with correct taxa and topology."""
    instructions =  """You are given images of phylogenetic trees and are tasked to respond with corresponding Newick 
    and just the Newick. 
    Respond with the corresponding Newick in correct format: 
    ((<taxon1>:<distance1>,<taxon2>:<distance2>):<clade_distance1>,<taxon3>=<distance3>); e.g.
    ((Sipunculidae:2.37,Tenebrio_molitor:1.55):4.58,((Siren:1.43,Catreus:3.63):0.27,Notothenioidei:1.66):4.07);.
    Make sure the Newick string includes the correct taxon names, all distances 
    (if scalebar or branch lenghts are present) and most importantly the correct tree topology.
    Make sure every opening parenthesis has a closing parenthesis. Make sure no parentheses pair is superfluous."""
    b64_image = encode_image(image_path)
    response = client.responses.create(
        # make model more deterministic
        temperature=0,
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
    # parse the newick from the reponse
    if re.search(r"\(.+;", response.output_text):
        newick = re.search(r"\(.+;", response.output_text).group()
        print(f"[{get_time()}] Generated Newick before postprocessing using AI:\n{newick}")
    else:
        exit(f"[{get_time()}] Error in generate_newick_format(): No Newick found in the response.")
    # postprocess if newick formatting is wrong
    if not ptan.is_newick(newick):
        print(f"[{get_time()}] generate_newick_from_image_directly: Formatting is wrong. Beginning postprocessing.")
        newick = correct_newick(image_path, model)
        # check if formatting is still wrong after postprocessing
        if not ptan.is_newick(newick):
            print(f"[{get_time()}] generate_newick_from_image_directly: AI-Postprocessing failed.") 
    else:
        print(f"[{get_time()}] generate_newick_from_image_directly: Formatting is correct. Skipping AI-postprocessing.")
    print(f"[{get_time()}] Generated Newick:\n{newick}")
    return newick

def extract_topology(image_path, model="gpt-4.1"):
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
    if not model in ["gpt-4.1", "o4-mini", "gpt-4o"]:
        exit(f"""[{get_time()}] Error in generate_newick_from_image_directly: Given model is not supported. Please choose 
             from gpt-4.1, o4-mini or gpt-4o""")
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
                   
    Each tab of indentation means one level further down in hierarchy of the phylogenetic tree. No indentation means it
    is the root. One space of indentation means that the clade or leaf is the child of the root. Two tabs mean that 
    the clade or leaf is the grandchild of the root and so on. 
    Make sure all distances are included, that correct formatting is correct and that the taxon names are correct."""
    b64_image = encode_image(image_path)
    response = client.responses.create(
        temperature=0,
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
    # TODO: extract newick manually from the topology
    topology = response.output_text
    print(topology)
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
        exit(f"[{get_time()}] Error in generate_newick_format(): No Newick found in the response.")
    return updated_newick

# def extract_newick_from_topology(topology):
    

def extract_newick_wrapper(extract_func, infile_path, outfile_path, model, file_id):
    """
    Wrapper for extracting newicks from images. Uses the given function on the given image using the given model 
    and saves the final newick at the specified outfile_path. If the outfile_path is a directory then a file ID is added
    for identification. If no outfile path is given then the newick is saved to the current working directory.

    Args:
        extract_func (function): function for extracting newick from image
        infile_path (str): path to the image or path to a directory that contains an image
        outfile_path (str): path where the newick is saved at 
        model (str): OpenAI model "gpt-4.1", "o4-mini", "gpt-4o"  
        file_id (str): string appended to the newick filename for identification
    """
    # if the infile path is a dir get the img from the dir
    if os.path.isdir(infile_path):
        image_from_dir = get_image_from_directory(infile_path)
        print(f"[{get_time()}] Given image: " + str(image_from_dir))
        print(f"[{get_time()}] Used model: " + model)
        generated_newick = extract_func(image_path=image_from_dir, model=model)
        print(f"[{get_time()}] Newick was generated.")
        write_newick_to_file(newick=generated_newick, model=model, file_id=file_id, outfile_path=outfile_path)
    # if infile path is a file then get the img directly
    elif os.path.isfile(infile_path):
        print(f"[{get_time()}] Given image: " + str(infile_path))
        print(f"[{get_time()}] Used model: " + model)
        generated_newick = extract_func(image_path=infile_path, model=model)
        write_newick_to_file(newick=generated_newick, model=model, outfile_path=outfile_path)
        print(f"[{get_time()}] Newick was saved to {outfile_path}")
    else:
        exit(f"[{get_time()}] Error in extract_newick_wrapper: infile path {infile_path} is not valid.")
                    
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
        choices=["extract_newick_directly", "extract_topology"],
        help="Specify which functionality you want to use.")
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
    # file ID 
    file_id = str(datetime.datetime.now().strftime(r"%H%M%S%f"))
    # Mode: generation of newick
    if approach == "extract_newick_directly":
        extract_newick_wrapper(
            extract_func=extract_newick_directly, 
            infile_path=infile_path, 
            outfile_path=outfile_path, 
            model=model, 
            file_id=file_id
        )
    elif approach == "extract_topology":
        extract_newick_wrapper(
            extract_func=extract_topology, 
            infile_path=infile_path, 
            outfile_path=outfile_path, 
            model=model, 
            file_id=file_id
        )
    else:
        exit(f"[{get_time()}] Error in main(): Approach {approach} not valid.")
            
# execute the main method
main()

# path to test directory for testing
# Random distances (0.00 and 20.00), 20 taxa: 
# "C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\data-collection\20_taxa_rand_distances"
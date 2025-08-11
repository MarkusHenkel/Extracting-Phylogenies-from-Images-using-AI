import datetime
import re
from ete3 import Tree
import logging
import os 

# logs
logger = logging.getLogger(__name__)

# time for logs
def get_time():
    return datetime.datetime.now().strftime("%Y-%b-%d %H:%M:%S")

# file ID 
def get_file_id():
    return str(datetime.datetime.now().strftime(r"%H%M%S%f"))

def get_image_format(image_path):
    """
    Using regex matches 2-4 chars at the end of of the path and returns the matched string if it is either jpeg, 
    png or jpg

    Args:
        image_path (str): path of the file

    Returns:
        str: image file format
    """
    if not os.path.isfile(image_path):
        raise ValueError("Filepath is not valid.")
    elif (file_ending := re.search(r"(?<=\.)\w{3,4}$", image_path).group()) in ["jpeg", "png", "jpg"]:
        return file_ending
    else:
        raise ValueError(f"Expected image file (jpeg, png, jpg) but got {file_ending}")
    
def taxids_from_newick(newick):
    """
    Given a newick tree containing taxon IDs return a list of all taxon IDs inside the string using regex.

    Args:
        newick_string (str): The tree in newick format

    Returns:
        List[str]: List containing all taxon IDs in order
    """
    return re.findall(r"[^)(,:][\d]+[^:]", newick)

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

def is_newick(newick, format=None):
    """
    Given a newick string checks if it has valid formatting. \n
    Formats:\n
    0:	flexible with support values.\n
    1:	flexible with internal node names.\n
    2:	all branches + leaf names + internal supports.\n
    3:	all branches + all names.\n
    4:	leaf branches + leaf names.\n
    5:	internal and leaf branches + leaf names.\n
    6:	internal branches + leaf names.\n
    7:	leaf branches + all names.\n
    8:  all names.\n
    9:  leaf names.\n
    100:    topology only.\n
    (according to ete3's Tree class)
    
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
    except Exception as e:
        logger.info(f"Newick not parseable: {e}")
        return False
    return True

def get_newick_format(newick):
    """
    Given a newick checks if it has format 100, 9, 5 or 0 in that order to accurately assign topology-only (100), 
    taxa-only (9), newicks with taxa and branch lenghts for internal nodes and leaf nodes (5) and newicks with support 
    values (0) their corresponding format

    Args:
        newick (str): newick with format 0, 5, 9 or 100

    Raises:
        ValueError: Newick isn't valid

    Returns:
        int: format
    """
    if is_newick(newick, format=100): # topology-only
        return 100
    elif is_newick(newick, format=9): # taxa-only
        return 9
    elif is_newick(newick, format=5): # taxa + branch lengths for leaves and internal nodes
        return 5
    elif is_newick(newick, format=0): # flexible
        return 0 
    else:
        raise ValueError("Newick is not valid.")
    
def remove_distances(newick):
    """
    Removes distances from the treerenders newick using regex.
    Regex expects them to be right after a ":".

    Args:
        newick (str): newick

    Returns:
        str: newick string without distances
    """
    return re.sub(r":\d+(\.\d+){0,1}", "", newick)
    
def remove_taxa_from_newick(newick):
    """
    Given a newick removes all taxa inside using regex.
    Regex expects them right after a '(' and right before ')', ',' and ':'.
    Special characters allowed inside the taxa: _ . - " ' / and #
    
    
    Args:
        newick (str): newick string with taxa, distances are optional

    Returns:
        str: newick string without taxa
    """
    return re.sub(r"(?<=[,(])[\d\w\.\-\"\'#\/]+(?=[\:\)\,])", "", newick)

def get_taxa(newick):
    """
    Returns all taxa from the given string in left to right order.
    Special characters allowed inside the taxa: _ . - " ' / and #

    Args:
        newick (str): Newick string

    Returns:
        List(str): List of all taxa found 
    """
    return re.findall(r"(?<=[,(])[\d\w\.\-\"\'#\/]+(?=[\:\)\,])", newick)

def remove_special_chars(taxon):
    r"""
    Deletes special characters that cause issues with the newick format or the packages from a given taxon.
    
    Warning for dataset module: Remove special characters after getting the topology otherwise the taxon might not be 
    found in the NCBI taxonomy 

    Args:
        taxon (str): taxon that may contain special character 

    Returns:
        str: taxon without special characters [](),;:'"\
    """
    return re.sub(r"[\[\]\,\;\:\'\"\\\(\)]", "", taxon) 

def is_balanced(newick):
    """
    Given a newick checks if its parentheses are balanced i.e. each closing parenthesis has a closing parenthesis 

    Args:
        newick (str): newick string

    Returns:
        bool: True if newicks parentheses are balanced 
    """
    opening_count = 0
    closing_count =  0
    for char in newick:
        if char == "(":
            opening_count += 1
        elif char == ")":
            closing_count += 1
    return True if opening_count == closing_count else False

def balance_parentheses(newick):
    """
    Naively solves unbalanced parentheses error in a given newick by placing opening or closing parentheses at the 
    beginning or end of the newick 

    Args:
        newick (str): newick string that may have unbalanced parentheses

    Returns:
        str: newick with balanced parentheses
    """
    newick.removesuffix(";")
    opening_count = 0
    closing_count =  0
    for char in newick:
        if char == "(":
            opening_count += 1
        elif char == ")":
            closing_count += 1
    if opening_count > closing_count:
        number_closing = opening_count - closing_count
        newick += ")" * number_closing
    if opening_count < closing_count:
        number_opening = closing_count - opening_count
        newick = ("(" * number_opening) + newick
    newick += ";"
    return newick

def is_multifurcating(node):
    if len(node.get_children()) > 2:
        return True
    for node in node.get_children():
        if is_multifurcating(node):
            return True
    return False

def count_multifurcations(newick):
    def count_multifurcations_worker(node):
        count = 0
        if len(node.get_children()) > 2:
            count += 1
        for child in node.get_children():
            count += count_multifurcations_worker(child) 
        return count
    return count_multifurcations_worker(Tree(newick))
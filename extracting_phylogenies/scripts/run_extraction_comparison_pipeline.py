import argparse
import subprocess
import os
import extracting_phylogenies.utilities.newick_util as ut
import sys
import logging
import re

# create logger
console_logger = logging.getLogger(__name__)

def main():
    arg_parser = argparse.ArgumentParser(
        description="Pipeline for iterating through dataset created by the dataset_creation module, extracting the "
        "newick using the newick_extraction_openai module and then directly creating a comparison directory with all "
        "comparison in a chosen location.")
    arg_parser.add_argument("-d","--dataset", required=True,help="Path to the dataset directory.")
    arg_parser.add_argument("-o","--outfile", required=True,
                            help="""Path to a .tsv where all comparisons are saved at. If file doesn't exist it is 
                            created.""")
    arg_parser.add_argument("-m", "--model", required=False,default="gpt-4.1",nargs="+",
                            help="""Choose the OpenAI model used for creating the newick. If multiple models are 
                            specified after --model e.g. '-m gpt-4.1 gpt-5' then for each model specified a newick is
                            generated and compared to the original. Options: gpt-4.1, gpt-4o, o4-mini, gpt-5 or 
                            gpt-4.1_finetuned""")
    arg_parser.add_argument(
        "-a",
        "--approach", 
        required=True, 
        choices=["extract_nwk", "extract_topo"],
        help="""Specify which functionality you want to use. extract_nwk extracts the newick directly i.e the model 
        produces text in newick format itself. extract_topo extracts the topology of the image and then the newick
        is extracted algorithmically.""")
    arg_parser.add_argument(
        "-f",
        "--format", 
        required=True, 
        choices=["taxa_only","topo_only","regular","no_format"],
        help="""Option to specify what the model should extract from the image. Depending on the format chosen different
        instructions are given to the model. Specify which format the output newick has: taxa_only: Output newick will    
        only have taxa and no branch lengths e.g. ((A,B,C),(D,E),((F,G),H)); topo_only: Output newick will have neither   
        taxa nor branch lengths, only topology e.g. ((,,),(,),((,),)); regular: Output newick will have taxa and branch   
        lengths e.g. ((A:2.37,B:1.55):4.58,((C:1.43,D:3.63):0.27,E:1.66):4.07); no_format: Based on the information in    
        the given image (taxa/no taxa, branch lengths/no branch lengths) the model decides on its own if the tree has     
        taxa and branch lengths or not and responds with a newick with a corresponding format. Default: no_format""")
    args = arg_parser.parse_args()
    dataset = args.dataset
    outfile = args.outfile
    format = args.format
    model_list = args.model
    approach = args.approach
    # check inputs
    if not os.path.exists(dataset):
        raise FileNotFoundError(f"--dataset: No directory found: {dataset}")
    if os.path.isdir(outfile):
        raise IsADirectoryError(f"--outfile: Expected filepath but got directory path: {outfile}")
    # configure logger
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(logFormatter)
    console_logger.addHandler(stream_handler)
    console_logger.setLevel(logging.INFO)
    console_logger.info('Starting run_extraction_comparison_pipeline')
    
    # get list of all subdirectories of the dataset
    data_dirs = [data for data in os.listdir(dataset) if os.path.isdir(os.path.join(dataset, data))]
    console_logger.info(f"Got {len(data_dirs)} sub-directories.")
    
    # iterate over all subdirectories of the given dataset
    for data in data_dirs:
        console_logger.info(f"Current subdirectory: {data}")
        # complete filepath of the current data dir
        data_path = os.path.join(dataset, data)
        # get files in data directory
        data_filenames = [file for file in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, file))]
        # get directories in data directory
        data_dirnames = [dir for dir in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, dir))]
        if (length:=len(data_dirnames)) > 1:
            raise ValueError(f"Expected 0 or 1 directories (predictions) in data directory but got {length}: {data}")
        # get the complete filepath of the current data dir
        root = os.path.join(dataset, data)
        # print("root" + str(root)) 
        newick_filename = None
        params_filename = None
        image_filename = None
        # get files for generating newick
        for file in data_filenames:
            # print("filename" + str(file))
            if file.endswith(".nwk"):
                newick_filename = file
            elif file.endswith(".tsv"):
                params_filename = file
            elif file.endswith(".jpg") or file.endswith(".jpeg") or file.endswith(".png"):  
                image_filename = file
        # check if all necessary files are present 
        if not newick_filename or not params_filename or not image_filename:
            raise ValueError(
                f"Expected image, newick and params in {root} but only got newick: {newick_filename}, params: "
                f"{params_filename}, image: {image_filename}"
            )
        
        # get the complete filepaths for newick, params and image
        newick_path = os.path.join(root, newick_filename)
        params_path = os.path.join(root, params_filename)
        image_path = os.path.join(root, image_filename)
        # get list of directories in the current data dir, should at most be 1 (predictions)
        subdirs = [dir for dir in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, dir))]
        if (len_subdirs := len(subdirs)) == 1:
            # get the name of the directory containing the predictions
            predictions = subdirs[0]
            # get the full directory path of the predictions directory
            predictions_path = os.path.join(data_path,predictions)
            # get all files in predictions directory and put their complete filepaths into a list
            filepaths_predictions = [os.path.join(predictions_path, file) for file in os.listdir(predictions_path) 
                                        if os.path.isfile(os.path.join(predictions_path, file))]
        elif len_subdirs > 1:
            raise ValueError(f"Expected one directory in {data} but got {len_subdirs}")
        else:
            filepaths_predictions = None
        
        # generate one newick for each model specified
        for model in model_list:
            # bool for potentially skipping subprocess
            skip_extraction = False
            # check if there already exists a prediction with the current parameters
            if filepaths_predictions:
                approach_abbr = "nwk" if approach == "extract_nwk" else "topo" if approach == "extract_topo" else None
                for prediction in filepaths_predictions:
                    pred_filename = os.path.split(prediction)[1]
                    if re.search(model, pred_filename) and re.search(approach_abbr,pred_filename):
                        skip_extraction = True
            # skip the extraction subprocess if the prediction already exists
            if not skip_extraction:
                console_logger.info(f"Extracting newick using {model}")
                subprocess.run(
                    [
                        "python", 
                        "-m", 
                        "extracting_phylogenies.newick_extraction_openai.newick_extraction_openai", 
                        "--approach",f"{approach}",
                        "--infile_path",f"{image_path}",
                        "--outfile_path",f"{root}",
                        "--model",f"{model}",      
                        "--format",f"{format}"           
                    ],
                    check=True   
                )
            else:
                console_logger.info(f"Skipping extraction for model {model} and approach {approach}")
        # get list of directories in the current data dir (should be one: predictions)
        subdirs = [dir for dir in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, dir))]
        # print(subdirs)
        if (len_subdirs := len(subdirs)) > 1:
            raise ValueError(f"Expected one directory in {data} but got {len_subdirs}")
        else:
            # get the name of the directory containing the predictions
            predictions = subdirs[0]
        # get the full directory path of the predictions directory
        predictions_path = os.path.join(data_path,predictions)
        # print(predictions_path)
        # get all files in predictions directory and put their complete filepaths into a list
        filepaths_predictions = [os.path.join(predictions_path, file) for file in os.listdir(predictions_path) 
                                 if os.path.isfile(os.path.join(predictions_path, file))]
        # print(filepaths_predictions)
        
        # create comparison for each model by iterating over all prediction filepaths
        for filepath in filepaths_predictions:
            # print(filepath)
            subprocess.run(
                [
                    "python", 
                    "-m", 
                    "extracting_phylogenies.newick_comparison.newick_comparison", 
                    "--original_newick", f"{newick_path}",
                    "--generated_newick", f"{filepath}",
                    "--outfile", f"{outfile}",
                    "--params", f"{params_path}"
                ],
                check=True
            )
    console_logger.info("Finished run_extraction_comparison_pipeline")
if __name__ == "__main__":
    main()
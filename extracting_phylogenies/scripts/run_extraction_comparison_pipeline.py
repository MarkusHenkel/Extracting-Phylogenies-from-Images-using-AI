import argparse
import subprocess
import os
import extracting_phylogenies.utilities.utilities as ut
import sys
import logging

# create logger
console_logger = logging.getLogger(__name__)

def main():
    arg_parser = argparse.ArgumentParser(
        description="Pipeline for iterating through dataset created by the dataset_creation module, extracting the "
        "newick using the newick_extraction_openai module and then directly creating a comparison directory with all "
        "comparison in a chosen location.")
    arg_parser.add_argument("-d","--dataset", required=True,help="Path to the dataset directory.")
    arg_parser.add_argument("-o","--outfile", required=True,help="Path to a directory where all comparison files are "
                            "saved at.")
    # arg_parser.add_argument("-m", "--model", required=False, choices=["gpt-4o", "gpt-4.1", "o4-mini"], type=str, 
    #                              default="gpt-4.1",
    #                              help="""Choose which OpenAI model is used in the generation of the newick string. 
    #                              Choose from GPT-4o, GPT-4.1 or o4-mini. Default model used is GPT-4.1.""")
    arg_parser.add_argument(
        "-a",
        "--approach", 
        required=True, 
        choices=["extract_nwk", "extract_topo"],
        help="""Specify which functionality you want to use. extract_nwk extracts the newick directly i.e the model 
        produces text in newick format itself. extract_topo extracts the topology of the image and then the newick
        is extracted algorithmically.""")
    arg_parser.add_argument("--quiet", action="store_true",help="On/Off flag. Option to turn off logging to console.")
    args = arg_parser.parse_args()
    dataset = args.dataset
    outfile = args.outfile
    # model = args.model
    approach = args.approach
    quiet = args.dataset
    
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
    
    # check if dataset path is valid 
    if not os.path.exists(dataset) or not os.path.exists(dataset):
        raise ValueError("Path is not valid.")
    
    # get list of all subdirectories of the dataset
    data_dirs = [data for data in os.listdir(dataset) if os.path.isdir(os.path.join(dataset, data))]
    print("data_dirs"+str(data_dirs))
    console_logger.info(f"Got {len(data_dirs)} sub-directories.")
    
    # iterate over all subdirectories of the given dataset
    for data in data_dirs:
        console_logger.info(f"Current subdirectory: {data}")
        # complete filepath of the current data dir
        data_path = os.path.join(dataset, data)
        # get files in data directory
        data_filenames = [file for file in os.listdir(data_path) if os.path.isfile(os.path.join(data_path, file))]
        print("data_filenames" + str(data_filenames))
        
        # get the complete filepath of the current data dir
        root = os.path.join(dataset, data)
        print("root" + str(root))
        
        newick_filename = None
        params_filename = None
        image_filename = None
        # get files for generating newick
        for file in data_filenames:
            print("filename" + str(file))
            if file.endswith(".nwk"):
                newick_filename = file
            elif file.endswith(".tsv"):
                params_filename = file
            elif file.endswith(".jpg") or file.endswith(".jpeg") or file.endswith(".png"):  
                image_filename = file
        print("newick_filename" + str(newick_filename))
        print("params_filename" + str(params_filename))
        print("image_filename" + str(image_filename))
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
        
        # generate one newick for each model
        for model in ["gpt-4.1", "gpt-4o", "o4-mini"]:
            console_logger.info(f"Extracting newick using {model}")
            subprocess.run(
                [
                    "python", 
                    "-m", 
                    "extracting_phylogenies.newick_extraction_openai.newick_extraction_openai", 
                    "--approach",
                    f"{approach}",
                    "--infile_path",
                    f"{image_path}",
                    "--outfile_path",
                    f"{root}",
                    "--model",
                    f"{model}",                 
                ]
            )
        # get predictions dir, if there is more than one dir in the current data dir, raise exception
        if (number_dirs := len(dirs := [dir for dir in os.listdir(data_path) 
                                        if os.path.isdir(os.path.join(data_path, dir))])) > 1:
            raise ValueError(f"Expected one directory in {data} but got {number_dirs}")
        else:
            predictions = dirs[0]
        # create comparison for each model by iterating over all predictions
        for prediction in predictions:
            prediction_path = os.path.join(data_path,"predictions",prediction)
            console_logger.info(f"Comparing {data}'s newick to {prediction}")
            subprocess.run(
                [
                    "python", 
                    "-m", 
                    "extracting_phylogenies.newick_comparison.newick_comparison", 
                    "--original_newick", 
                    f"{newick_path}",
                    "-generated_newick",
                    f"{prediction_path}"
                    "--outfile",
                    f"{outfile}",
                    "--params",
                    f"{params_path}"
                ]
            )
    console_logger.info("Finished")
if __name__ == "__main__":
    main()
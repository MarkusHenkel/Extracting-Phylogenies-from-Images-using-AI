from openai import OpenAI
import json
import os
from extracting_phylogenies.newick_extraction_openai import instructions as instr
from extracting_phylogenies.image_augmentation import image_augmentation as aug

import base64

# set client with API key
client = OpenAI(api_key = os.getenv("BA_API_KEY"))

def load_dataset(dataset):
    # create list for image newick tuples
    img_nwk_pairs = []
    # list with all data subdirectories
    data_dirs = [data for data in os.listdir(dataset) if os.path.isdir(os.path.join(dataset, data))]
    # iterate over subdiretories
    for data_dir in data_dirs:
        # complete filepath of the current data dir
        data_path = os.path.join(dataset, data_dir)
        # get files in data directory
        data_filepaths = [os.path.join(data_path, file) for file in os.listdir(data_path) 
                          if os.path.isfile(os.path.join(data_path, file))]
        nwk_path = [file for file in data_filepaths if os.path.split(file)[1].endswith("nwk")][0]
        img_path = [file for file in data_filepaths if os.path.split(file)[1].endswith("jpg")][0]
        with open(nwk_path, "r") as nwk_file: 
            nwk = nwk_file.read()
        img_nwk_pairs.append((img_path, nwk))
    return img_nwk_pairs

def create_jsonl_entry_base64(base64_string, prompt, instructions, truth):
    """
    Create entry for a JSONL for openai finetuning. The base64-encoded image is put into a data URL,
    truth is the original newick.

    Args:
        base64_string (str): base64-encoded image
        prompt (str): model prompt
        instructions (str): model instructions
        truth (str): original newick

    Returns:
        json: json entry for openai finetuning containing prompt, instructions, image, truth
    """
    entry={
        "messages":[
            {"role":"system","content":instructions},
            {"role":"user","content":prompt},
            {"role":"user","content":[
                {
                    "type":"image_url",
                    "image_url":{"url": f"data:image/jpeg;base64,{base64_string}"}
                }
            ]
            },
            {"role": "assistant", "content": truth}
        ]
    }
    return json.dumps(entry)


augment= aug.pad_resize_augment_wrapper(image_size=1024, normalize=False, to_tensor=False)

def create_jsonl_file_base64(dataset_path, jsonl_dirpath, augment=augment, chunk_size=100):
    """
    Takes a path to a dataset and instead of creating one huge jsonl, splits up the jsonl into smaller jsonls of at 
    most <chunk_size> image/nwk pairs to fly under the maximum jsonl size

    Args:
        dataset_path (str): path to dataset
        jsonl_dirpath (str): path to directory where jsonl chunks are saved
        chunk_size (int): number 
    """
    # load dataset 
    dataset = load_dataset(dataset_path)
    # set counter for jsonl entries
    entry_count = 0
    # set file counter
    file_count = 1
    # initialize the file contents of current 
    contents = "" 
    for i in range(len(dataset)):
        img_newick_pair = dataset[i]
        print(f"File counter: {file_count}")
        print(f"Entry counter: {entry_count}")
        # get prompt and instructions from instructions file i.e. let model produce newick directly
        prompt = instr.prompt
        instructions = instr.instr_nwk_regular
        # get image path
        img_path = img_newick_pair[0]
        # get newick
        newick = img_newick_pair[1]
        # augment image before uploading
        augmented_image = aug.get_augmented_image(augmentation=augment,image_path=img_path)
        # base64-encode image for the jsonl entry
        base64_string = base64.b64encode(augmented_image).decode("utf-8")
        # set json outfile path 
        jsonl_outfile = os.path.join(jsonl_dirpath, f"data{file_count}")
        # append the current entry
        contents+=create_jsonl_entry_base64(base64_string=base64_string,prompt=prompt, instructions=instructions, truth=newick)
        contents+="\n"
        entry_count += 1
        # make chunks of at most size 100
        if entry_count >= chunk_size or i == len(dataset)-1:
            # create the file if it doesnt exist, otherwise append entries  
            with open(jsonl_outfile, "w") as jsonl:
                jsonl.write(contents)        
            print(f"length contents: {len(contents)}")
            file_count += 1
            entry_count = 0
            contents = ""

    
def upload_jsonl_to_openai(jsonl_path):
    """
    Uploads the jsonl and returns the file ID for finetuning

    Args:
        jsonl_path (str): path to jsonl

    Returns:
        str: file ID
    """
    if not os.path.exists(jsonl_path):
        FileNotFoundError(f"File not found: {jsonl_path}")
    response = client.files.create(
        file=open(jsonl_path, "rb"),
        purpose="fine-tune"
    )
    return response.id

def upload_chunks(jsonl_dir):
    chunks = [os.path.join(jsonl_dir, file) for file in os.listdir(jsonl_dir) 
                                        if os.path.isfile(os.path.join(jsonl_dir, file))]
    if chunks:
        for chunk in chunks:
            print(f"Chunk: {chunk}")
            file_id = upload_jsonl_to_openai(chunk)
            print(file_id)
    else:
        raise ValueError(f"No chunks to upload, dir has no files: {jsonl_dir}")
        

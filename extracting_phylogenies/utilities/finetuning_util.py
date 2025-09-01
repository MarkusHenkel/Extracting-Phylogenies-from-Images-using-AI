from PIL import Image 
from datasets import Dataset
import re, os
from extracting_phylogenies.image_augmentation import image_augmentation as augm

# maximum amount of tokens in one newick
max_token_length = 512

def load_dataset(dataset_path, torch_dataset=False):
    """
    Given a path to a directory, iterates through all subdirs and creates a dict with three lists 
    IDs, paths to each image and the newicks. Ordering is strict, the first image in the image list corresponds 
    to the first newick in the newick list etc.
    If hf_dataset is True the dict is turned into a Dataset object

    Args:
        dataset_path (str): path to the dataset with subdirs that each have an id e.g. data202541627 and 
        contain an image_path/newick pair

    Returns:
        dict: dict that contains 3 lists: IDs, paths to each image, newicks
    """
    dataset = {
        "id": [],
        "image": [],
        "label": []
    }
    if not os.path.isdir(dataset_path):
        raise FileNotFoundError(f"Given path is not a directory: {dataset_path}")
    # iterate over each subdirectory in the dataset directory
    for subdir, _, files in os.walk(dataset_path):
        image = ""
        newick = ""
        # skip root of os.walk, the dataset itself which has no id etc because it only contains subdirs with pairs
        if subdir == dataset_path:
            continue
        # assign each entry in the dataset an ID; looks for a number at the end of the name of the current subdir
        if re.search(r"\d+$", subdir):
            id = re.search(r"\d+$", subdir).group()
        else: 
            # skip over subdirectories without an id => predictions
            continue
        for file in files:
            # 
            if file.endswith(".jpg"):
                image = os.path.join(subdir, file)
            # accept nwk or txt in case topology is provided for the ground truth
            elif file.endswith(".nwk") or file.endswith(".txt"):
                with open(os.path.join(subdir, file), "r") as nwk_file:
                    newick = nwk_file.read()
        if id and image and newick:
            dataset["id"].append(id)
            dataset["image"].append(image)
            dataset["label"].append(newick)
        else:
            print(f"Skipping subdir {subdir} because id, image or newick are missing.")
        # reset all variables so that missing values in the next iteration arent assigned previous values
        id = ""
        image = ""
        newick = ""
    if torch_dataset:
        return Dataset.from_dict(dataset)
    else:
        return dataset
    

from extracting_phylogenies.newick_extraction_openai import instructions as instr

def format_data(sample, instructions=None, prompt=None):
    if not instructions:
        instructions = instr.instr_nwk_regular
    if not prompt:
        prompt = instr.prompt
    return [
        {
            "role": "system",
            "content": [{"type": "text", "text": instructions}],
        },
        {
            "role": "user",
            "content": [
                {
                    "type": "image",
                    "image": sample["image"],
                },
                {
                    "type": "text",
                    "text": prompt,
                },
            ],
        },
        {
            "role": "assistant",
            "content": [{"type": "text", "text": sample["label"]}],
        },
    ]

def format_data_inference(image, instructions, prompt):
    return [
        {
            "role": "system",
            "content": [{"type": "text", "text": instructions}],
        },
        {
            "role": "user",
            "content": [
                {
                    "type": "image",
                    "image": image,
                },
                {
                    "type": "text",
                    "text": prompt,
                },
            ],
        }
    ]

class ChatbotDataset(Dataset):
    """
    Class for creating a chatbot dataset from a Dataset object that contains completely formatted messages as entries 
    containing instructions,prompt, image and label. Initialized with an augmentation this class provides on-the-fly 
    image augmentation to avoid having to augment all images of a dataset before fine-tuning.
    Inherits from Dataset. 

    Unfortunately SFT trainer by huggingface expects a huggingface dataset not a torch dataset
    """
    def __init__(
        self, 
        dataset, 
        instructions,
        prompt,
        augmentation=None,
        ):
        self.dataset = dataset
        self.instructions = instructions
        self.prompt = prompt
        self.augmentation = augmentation
    
    def __len__(self):
        return len(self.dataset)
    
    
    def __getitem__(self, index):
        # load entry at index
        entry = self.dataset[index]
        # get image for some reason sft trainer returns the image as a list with one element
        image = entry["image"][0]
        # get label
        label = entry["label"][0]
        if self.augmentation:
            # augment image
            image = augm.get_augmented_image(augmentation=self.augmentation, image_path=image)
        return format_data(image, label, self.instructions, self.prompt)
    
    
from extracting_phylogenies.utilities import finetuning_util as util
from datasets import Image


def upload_dataset_to_hf(dataset_path, hf_username, hf_dataset_name, get_images=False):
    dataset = util.load_dataset(dataset_path, torch_dataset = True)    
    if get_images:
        dataset = dataset.cast_column("image", Image(decode=True))
    dataset.push_to_hub(f"{hf_username}/{hf_dataset_name}")
    

import gc
import torch


def clear_memory():
    # Delete variables if they exist in the current global scope
    if "inputs" in globals():
        del globals()["inputs"]
    if "model" in globals():
        del globals()["model"]
    if "processor" in globals():
        del globals()["processor"]
    if "trainer" in globals():
        del globals()["trainer"]
    if "peft_model" in globals():
        del globals()["peft_model"]
    if "bnb_config" in globals():
        del globals()["bnb_config"]
    gc.collect()
    torch.cuda.empty_cache()
    torch.cuda.synchronize()
    gc.collect()
    print(f"GPU allocated memory: {torch.cuda.memory_allocated() / 1024**3:.2f} GB")
    print(f"GPU reserved memory: {torch.cuda.memory_reserved() / 1024**3:.2f} GB")
    
    
from qwen_vl_utils import process_vision_info

def generate_text_from_sample(model, processor, sample, max_new_tokens=1024):
    # Prepare the text input by applying the chat template
    text_input = processor.apply_chat_template(sample, tokenize=False, add_generation_prompt=True)
    # debug
    # print("############## text_input ############## ")
    # print(text_input)
    # Process the visual input from the sample
    image_inputs, _ = process_vision_info(sample)
    # print("############## image inputs ##############")
    # print(image_inputs)
    # Prepare the inputs for the model
    model_inputs = processor(
        text=[text_input],
        images=image_inputs,
        return_tensors="pt",
    ).to(model.device)
    # print("############## model inputs ##############")
    # print(model_inputs)
    # Generate text with the model
    generated_ids = model.generate(**model_inputs, min_new_tokens=64, max_new_tokens=max_new_tokens,)
    # Decode the output text
    output_text = processor.batch_decode(
        generated_ids, skip_special_tokens=True, clean_up_tokenization_spaces=False
    )
    return output_text[0] 
    
from PIL import Image 
from datasets import Dataset
import re, os

# maximum amount of tokens in one newick
max_token_length = 512

def load_dataset(dataset_path):
    """
    Given a path to a directory, iterates through all subdirs and creates a dict with three lists 
    IDs, paths to each image and the newicks. Ordering is strict, the first image in the image list corresponds 
    to the first newick in the newick list etc

    Args:
        dataset_path (str): path to the dataset with subdirs that each have an id e.g. data202541627 and 
        contain an image_path/newick pair

    Returns:
        dict: dict that contains 3 lists: IDs, paths to each image, newicks
    """
    dataset = {
        "id": [],
        "image": [],
        "newick": []
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
            # skip over subdirectories without an id (TODO: maybe subdirs shouldnt be forced to have an ID, only img, nwk etc)
            print(f"Warning: Skipping subdir {subdir} because it does not have an ID in its name.")
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
            dataset["newick"].append(newick)
            id = ""
            image = ""
            newick = ""
        else:
            print(f"Skipping subdir {subdir} because id, image or newick are missing.")
            id = ""
            image = ""
            newick = ""
    return dataset

def preprocess_dataset(dataset, tokenizer, image_processor):
    """
    Given a dict with 3 lists: ID, image path and newick returns a huggingface dataset object 
    (https://huggingface.co/docs/datasets/v1.2.1/loading_datasets.html) containing the pixel values of the images and
    tokenized newicks
    
    Args:
        dataset (dict): dict with 3 lists: IDs, image paths, newicks

    Returns:
        dataset object: huggingface dataset object containing pixel values and tokenized newicks
    """
    # rename keys image and newick in the dataset dict (model expects keys to be pixel_values and labels)
    dataset["pixel_values"] = dataset.pop("image")
    dataset["labels"] = dataset.pop("newick")
    # replace image and newick with the entries the model expects => pixel_values and labels
    for i in range(len(dataset["id"])):
        pil_image = Image.open(dataset["pixel_values"][i])
        # get the pixel values of the pil_image as a numpy array (for albumentations)
        pixel_values = image_processor(pil_image, return_tensors="pt").pixel_values[0].numpy()
        # Debugging
        # if i == 0:
        #     test = pixel_values
        # if not str(type(pixel_values)) == "<class 'torch.Tensor'>":
        #     print(100 * "#")
        #     print("Not a tensor")
        #     print(100 * "#")
        label = tokenizer(
            dataset["labels"][i], 
            padding="max_length", # pad up to max length if necessary
            max_length=max_token_length, # set maximum of tokens in each newick string
            return_tensors="pt" # tensor is a pytorch tensor
        ).input_ids[0].numpy()
        dataset["pixel_values"][i] = pixel_values
        dataset["labels"][i] = label
    # Debugging
    # print(100 * "#")
    # print("Pixel values before turning them into a numpy array:")
    # print(test)
    # print(type(test))
    # print(100 * "#")
    # print(100 * "#")
    # print("Pixel values after turning them into a numpy array:")
    # print(dataset["pixel_values"][0])
    # print(type(dataset["pixel_values"][0]))
    # print(100 * "#")
    return Dataset.from_dict(dataset)
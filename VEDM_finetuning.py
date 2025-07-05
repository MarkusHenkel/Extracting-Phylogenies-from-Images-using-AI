# model and model configs
from transformers import VisionEncoderDecoderModel, GPT2TokenizerFast, ViTImageProcessor
import albumentations as A # for image augmentation
import matplotlib.pyplot as plt
import cv2 # 
import os # for getting dataset from filesystem
import argparse
import requests
from PIL import Image # 
import torch
import datetime

# time for logging
time = datetime.datetime.now().strftime("%Y-%b-%d %H:%M:%S")

# image size in px
image_size = 2048

# load a fine-tuned image captioning model and corresponding tokenizer and image processor
model = VisionEncoderDecoderModel.from_pretrained("nlpconnect/vit-gpt2-image-captioning")
tokenizer = GPT2TokenizerFast.from_pretrained("nlpconnect/vit-gpt2-image-captioning")
image_processor = ViTImageProcessor.from_pretrained("nlpconnect/vit-gpt2-image-captioning")

"""
Basic transform for training the baseline. Uniform res of <image_size>x<image_size>. 
Breaks aspect ratio!
"""
transform_1 = A.Resize(image_size,image_size)
"""
Basic transform for training the baseline. Uniform res of <image_size>x<image_size>.
Keeps aspect ratio!
"""
transform_1_ar = A.Compose(
    [
        A.LongestMaxSize(max_size=image_size, interpolation=1), # resize while keeping aspect ratio
        A.PadIfNeeded(min_height=image_size, min_width=image_size, border_mode=0) # padding for square res
    ]
)

"""
Resize images to uniform resolution and augment them to get more samples and make model more robust
towards different conditions of the image. Only used after baseline is established i.e. second training round
"""
transform_2 = A.Compose(
    [
        A.LongestMaxSize(max_size=image_size, interpolation=1), # resize while keeping aspect ratio
        A.PadIfNeeded(min_height=image_size, min_width=image_size, border_mode=0), # padding for square res
        A.Rotate(limit=10), # rotation may cause taxa to be outside the visible frame!
        A.OneOf([
            # holes not squre, possibly overlapping, not in grid formation
            A.CoarseDropout(num_holes_range=(1, 8), hole_height_range=(0.01, 0.1),hole_width_range=(0.01, 0.1),  p=1.0),
            # square holes in grid formation, smaller, fewer holes, non-overlapping
            A.GridDropout(ratio=0.1, p=1.0),
            # more square holes in grid format, non-overlapping
            A.GridDropout(ratio=0.2, p=1.0),
            # fewer but bigger square holes, , non-overlapping
            A.GridDropout(ratio=0.1, unit_size_range=(500, 1000),p=1.0)
        ], p=0.1), # apply dropout only every 10th image
        
    ]
)

def load_dataset(dataset_path):
    """
    Given a path to a directory, iterates through all subdirs and add the pairs as tuples to a list

    Args:
        dataset_path (str): path to the dataset with subdirs each containing a image/newick pair

    Returns:
        List((str, str)): list of image/newick tuples
    """
    dataset = []
    for subdir, dirs, files in os.walk(dataset_path):
        image = ""
        newick = ""
        for file in files:
            # support other file formats as well in case ptan module uses other formats in the future
            if file.endswith(".jpg") or file.endswith(".png") or file.endswith(".jpeg") or file.endswith(".pdf"):
                image = file
            elif file.endswith(".nwk") or file.endswith(".txt"):
                with open(os.path.join(subdir, file), "r") as nwk_file:
                    newick = nwk_file.read()
                # TODO: warn user if dataset contains false newicks, check with is_newick from ptan
        if image and newick:
            dataset.append((image, newick))
        else:
            print(f"[{time}] Warning in load_dataset(): Subdir {subdir} did not contain image/nwk or image/txt pair.") 
    return dataset

def get_augmented_img_from_path(transform, image_path):
    try:
        image = cv2.imread(image_path)
    except:
        exit(f"[{time}] Error in get_augmented_img_from_path: Path could not be resolved.")
    augmented = transform(image=image)['image']
    augmented_rgb = cv2.cvtColor(augmented, cv2.COLOR_BGR2RGB)
    return augmented_rgb

def test_transforms(transform, image_path):
    augmented_rgb = get_augmented_img_from_path(transform, image_path)
    plt.figure(figsize=[image_size, image_size], dpi=100)
    plt.imshow(augmented_rgb)
    plt.axis('off')
    plt.show()
    
def perform_inference(image_path):
    image = cv2.imread(image_path)
    image_preprocessed = transform_1_ar(image=image)["image"]
    pixel_values = image_processor(image_preprocessed, return_tensors="pt").pixel_values

    # autoregressively generate caption (uses greedy decoding by default)
    generated_ids = model.generate(pixel_values)
    generated_text = tokenizer.batch_decode(generated_ids, skip_special_tokens=True)[0]
    print(generated_text)

def main():
    # argument parser for getting dataset
    parser = argparse.ArgumentParser(description="Module for finetuning a Vision Encoder Decoder Model. Dataset must be\
        created using module phylogenetic_tree_and_newicks.py")
    # arguments
    parser.add_argument("-d", "--dataset", required=False, type=str, help="Path to the dataset made with the\
        phylogenetic_trees_and_newicks.py module.")
    # get user arguments
    args = parser.parse_args()
    path_dataset = args.dataset 
    if path_dataset:
        dataset = load_dataset(path_dataset)
    else:
        print(f"[{time}] Error: Dataset could not be loaded.")
    ########## Testing ##########
    # ete3_20_taxa = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\data_collection\dataset\ete3_20_taxa\data_2025-06-25-00-25-56-658225\ete3_tree_2025-06-25-00-25-56-658225.jpg"
    # phylo_20_taxa = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\data_collection\dataset\phylo_20_taxa\data_2025-06-25-00-25-33-042017\phylo_tree_2025-06-25-00-25-33-042017.jpg"
    # phylo_rand_dist = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\data_collection\dataset\phylo_rand_distances\data_2025-06-25-00-24-55-590431\phylo_tree_2025-06-25-00-24-55-590431.jpg"
    # ete3_rand_dist = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\data_collection\dataset\ete3_rand_distances\data_2025-06-25-00-22-24-082030\ete3_tree_2025-06-25-00-22-24-082030.jpg"
    # # example images for transformation testing
    # test transformations
    # test_transforms(transform_2, ete3_20_taxa)
# execute main
main()



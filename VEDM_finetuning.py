# model and model configs
from transformers import VisionEncoderDecoderModel, AutoTokenizer, ViTImageProcessor, TrainingArguments, Trainer
from datasets import Dataset
import albumentations as A # for image augmentation
import matplotlib.pyplot as plt
import cv2 # 
import os # for getting dataset from filesystem
# import requests
from PIL import Image # 
# import torch
import datetime
import re
import numpy as np # jumping between np arrays (albumentations) and pytorch tensors (model)

# time for logging
time = datetime.datetime.now().strftime("%Y-%b-%d %H:%M:%S")

# image size in px
image_size = 2048

# maximum amount of tokens in one newick
max_token_length = 512

########## IMAGE AUGMENTATION ##########
# Basic transform for training the baseline. Uniform res of <image_size>x<image_size>. Breaks aspect ratio!
augment_1 = A.Resize(image_size,image_size)

# Basic transform for training the baseline. Uniform res of <image_size>x<image_size>. Keeps aspect ratio!
augment_1_ar = A.Compose(
    [
        A.LongestMaxSize(max_size=image_size, interpolation=1), # resize while keeping aspect ratio
        A.PadIfNeeded(min_height=image_size, min_width=image_size, border_mode=0), # padding for square res
        A.Normalize(),
        A.ToTensorV2() # output tensor with shape (Channels, Height, Width)
    ]
)

# Resize images to uniform resolution and augment them to get more samples and make model more robust 
# towards different conditions of the image. Only used after baseline is established i.e. second training round
augment_2 = A.Compose(
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
            A.GridDropout(ratio=0.1, unit_size_range=(500, 1000),p=1.0),
        ], p=0.1), # apply one type of dropout every 10th image TODO: every image instead in second training round
        A.Normalize(), # TODO: default values or not?
        A.ToTensorV2() # output tensor with shape (Channels, Height, Width)
        
    ]
)

def transform_1_ar(dataset):
    """
    transform wrapper for dataset object's .set_transform() 

    Args:
        dataset (huggingface dataset object): dataset created by the preprocess_dataset function

    Returns:
        _type_: _description_
    """
    # apply transform to all images in the dataset
    for i in range(len(dataset["pixel_values"])):
        # albumentations expects np array of shape (H,W,C) not (C,H,W) 
        pixel_values = np.transpose(dataset["pixel_values"][i], (1,2,0))
        dataset["pixel_values"][i] = augment_1_ar(image=pixel_values) # augment outputs tensors of shape (H,W,C) again
    return dataset

########## LOADING DATASET ##########
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
        exit(f"[{time}] Error in load_dataset: Given path is not a directory.")
    # iterate over each subdirectory in the dataset directory
    for subdir, dirs, files in os.walk(dataset_path):
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
            print(f"[{time}] Warning: Skipping subdir {subdir} because it does not have an ID in its name.")
            continue
        for file in files:
            # support other file formats as well in case ptan module uses other formats in the future
            if file.endswith(".jpg") or file.endswith(".png") or file.endswith(".jpeg") or file.endswith(".pdf"):
                image = os.path.join(subdir, file)
            elif file.endswith(".nwk") or file.endswith(".txt"):
                with open(os.path.join(subdir, file), "r") as nwk_file:
                    newick = nwk_file.read()
                # TODO: warn user if dataset contains false newicks?
        if id and image and newick:
            dataset["id"].append(id)
            dataset["image"].append(image)
            dataset["newick"].append(newick)
            id = ""
            image = ""
            newick = ""
        else:
            print(f"[{time}] Warning in load_dataset(): Subdir {subdir} is missing image path or Newick.")
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
        pixel_values = image_processor(pil_image, return_tensors="pt").pixel_values[0]
        label = tokenizer(
            dataset["labels"][i], 
            padding="max_length", # pad up to max length if necessary
            max_length=max_token_length, # set maximum of tokens in each newick string
            return_tensors="pt" # tensor is a pytorch tensor
        ).input_ids[0]
        # change tensor to np array for albumentations
        dataset["pixel_values"][i] = pixel_values.numpy()
        dataset["labels"][i] = label
    return Dataset.from_dict(dataset)

########## COMPUTE METRICS ##########
def compute_metrics(pred):
    # TODO: add metrics used in module comparison_of_newicks.py
    preds, label_ids = pred
    print("###############")
    print("preds:" + preds)
    print("label_ids:" + label_ids)
    print("###############")
      
########## FOR TESTING ##########
def get_augmented_img_from_path(image_path, augmentation=augment_1_ar):
    """
    Given the path to an image, opens it using OpenCV and applies a given albumentations transform.

    Args:
        image_path (str): path to image
        transform (albumentations transform, optional): Transform applied to given image. Defaults to augment_1_ar.

    Returns:
        pytorch tensor: numpy tensor containing the transformed image 
    """
    try:
        # loads image as numpy array in BGR format instead of RGB
        image = cv2.imread(image_path)
    except:
        exit(f"[{time}] Error in get_augmented_img_from_path: Path {image_path} could not be resolved.")
    # turn in into RGB format
    image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    # apply transform 
    augmented_tensor = augmentation(image=image_rgb)['image'] # tensor with shape (C, H, W)
    return augmented_tensor

def test_augmentation(augment, image_path):
    augmented_tensor = get_augmented_img_from_path(augmentation=augment, image_path=image_path)
    plt.figure(figsize=[image_size, image_size], dpi=100)
    # for matplotlib: turn the tensor into numpy array and transpose from (Channels, Height, Width) to (H, W, C)  
    plt.imshow(np.transpose(augmented_tensor.cpu().numpy(),(1,2,0))) 
    plt.axis('off')
    plt.show()
    
def get_number_tokens(newick, tokenizer):
    return len(tokenizer(newick, return_tensors="pt").input_ids[0])

########## PERFORMING INFERENCE ##########
def perform_inference_img_processor(image_path, model, image_processor, tokenizer):
    """
    Perform inference on an image using an image_processor e.g. the pretrained image processor used in finetuning

    Args:
        image_path (str): path to image
        image_processor (image processor object, optional): image processor e.g. HF ViTImageProcessor. Defaults to pretrained_image_processor.
    """
    image = cv2.imread(image_path)
    augmented_tensor = augment_1_ar(image=image)["image"] # tensor with shape (C, H, W)
    # turn the tensor into numpy array, transpose from (C, H, W) to (H, W, C) 
    augmented_np_arr = np.transpose(np.array(augmented_tensor),(1,2,0)) # TODO: test if it works this way
    pixel_values = image_processor(augmented_np_arr, return_tensors="pt").pixel_values
    generated_ids = model.generate(pixel_values)
    generated_text = tokenizer.batch_decode(generated_ids, skip_special_tokens=True)[0]
    print(generated_text)
    
def main():
    if __name__ == "__main__":
        # Finetuning
        # load a fine-tuned image captioning model, corresponding tokenizer, image processor and trainer
        model = VisionEncoderDecoderModel.from_pretrained("nlpconnect/vit-gpt2-image-captioning")
        pretrained_image_processor = ViTImageProcessor.from_pretrained("nlpconnect/vit-gpt2-image-captioning")
        tokenizer = AutoTokenizer.from_pretrained("nlpconnect/vit-gpt2-image-captioning") # TODO: gpt2 from huggingface instead?
        special_tokens = ["(", ")", ";", ",", ":"] # TODO: add special tokens or not?
        tokenizer.add_tokens(special_tokens, special_tokens=True)
        model.decoder.resize_token_embeddings(len(tokenizer))
        dataset_path = r".\data_collection\dataset_test\ete3_20_taxa"
        # r"""
        dataset = preprocess_dataset(
            dataset=load_dataset(dataset_path=dataset_path), 
            tokenizer=tokenizer, 
            image_processor=pretrained_image_processor
        )
        print(f"[{time}] main(): Dataset was loaded and preprocessed.")
        # use set_transform to apply image augmentations on the fly during training instead of map() to save disk space
        # train in 2 stages: baseline and then robustness => augment_1_ar then augment_2 (with dropout etc.)
        # TODO: album expects a np array not a tensor, turn tensor back into np array
        dataset.set_transform(transform=transform_1_ar)
        print(f"[{time}] main(): Dataset transform set to transform_1_ar.")
        # split dataset into 80% training and 20% evaluation
        split_dataset = dataset.train_test_split(test_size=0.2)
        train_dataset = split_dataset["train"]
        eval_dataset = split_dataset["test"]
        # set the training arguments 
        # taken from basic training example by huggingface: https://deepwiki.com/huggingface/transformers/3.1-trainer-class
        training_args = TrainingArguments(
            output_dir="./results",
            learning_rate=5e-5,
            per_device_train_batch_size=16,
            per_device_eval_batch_size=64,
            num_train_epochs=3,
            weight_decay=0.01,
            # TODO: load_best_model_at_end = True?
            # TODO: evaluation_strategy = "steps"?
        )

        trainer = Trainer(
            model=model,
            args=training_args,
            train_dataset=train_dataset,
            # TODO: data collator?
            eval_dataset=eval_dataset,
            tokenizer=tokenizer,
            # compute_metrics=compute_metrics, TODO: add compute_metrics function
        )
        print(f"[{time}] main(): Trainer was instantiated.")

        # Train the model
        train_results = trainer.train()
        print(f"[{time}] main(): trainer.train() finished.")
        trainer.save_model()
        trainer.log_metrics("train", train_results.metrics)
        trainer.save_metrics("train", train_results.metrics)
        trainer.save_state()
        # """
        ########## Testing ##########
        ete3_20_taxa = r".\data_collection\dataset_test\ete3_20_taxa\data_2025-06-25-00-25-53-945410\ete3_tree_2025-06-25-00-25-53-945410.jpg"
        phylo_20_taxa = r".\data_collection\dataset_test\phylo_20_taxa\data_2025-06-25-00-25-27-260802\phylo_tree_2025-06-25-00-25-27-260802.jpg"
        # # example images for transformation testing
        # long_newick = "(((((Sinorhizobium_fredii:1,Clostridium_baratii:1):1,((((Schaalia_odontolytica:1,Halanaerobium:1):1,Pylaiella:1):1,Trianthema:1):1,Arenaria:1):1):1,Bertholletia:1):1,Pastinaca:1,(((Trochodendron_aralioides:1,Alternaria_alternata:1):1,(Malacosoma_americanum:1,Xenopus_borealis:1):1):1,Coluber_constrictor:1):1):1,(((Afropavo_congensis:1,Paridae:1):1,Alouatta_seniculus:1):1,((Afipia_carboxidovorans:1,Desulfomicrobium_norvegicum:1):1,Centruroides_tecomanus:1):1):1);"
        # print(get_number_tokens(long_newick)) # 236 tokens without special tokens, 261 tokens with special tokens added 
        # test_augmentation(augment=augment_1_ar, image_path=phylo_20_taxa)
        perform_inference_img_processor(image_path=ete3_20_taxa, model=model, image_processor=pretrained_image_processor, tokenizer=tokenizer)
        
# execute main
main()



# model and model configs
from transformers import VisionEncoderDecoderModel, AutoTokenizer, ViTImageProcessor, TrainingArguments, Trainer
import cv2 # 
import numpy as np # jumping between np arrays (albumentations) and pytorch tensors (model)
# import a data collator
from transformers import default_data_collator
from extracting_phylogenies.image_augmentation import image_augmentation as augm


# image size in px
image_size = 2048

def transform_pad_resize(dataset):
    pad_resize = augm.pad_resize_augment_wrapper(to_tensor=True)
    # apply transform to all images in the dataset
    for i in range(len(dataset["pixel_values"])):
        # albumentations expects np array of shape (H,W,C) not (C,H,W) 
        pixel_values = np.transpose(dataset["pixel_values"][i], (1,2,0))
        dataset["pixel_values"][i] = pad_resize(image=pixel_values)["image"] # augment outputs tensors of shape (H,W,C) again
    # Debugging
    # print(100 * "#")
    # print("pixel_values AFTER AUGMENTATION")
    # print(dataset["pixel_values"][0])
    # print("TYPE OF pixel_values AFTER AUGMENTATION")
    # print(type(dataset["pixel_values"][0]))
    # print(100 * "#")
    return dataset

def transform_dropout(dataset):
    dropout = augm.dropout_augment_wrapper(to_tensor=True)
    # apply transform to all images in the dataset
    for i in range(len(dataset["pixel_values"])):
        # albumentations expects np array of shape (H,W,C) not (C,H,W) 
        pixel_values = np.transpose(dataset["pixel_values"][i], (1,2,0))
        dataset["pixel_values"][i] = dropout(image=pixel_values)["image"] # augment outputs tensors of shape (H,W,C) again
    return dataset


########## COMPUTE METRICS ##########
def compute_metrics(pred):
    # TODO: add metrics used in module newick_comparison.py
    preds, label_ids = pred
    print("###############")
    print("preds:" + preds)
    print("label_ids:" + label_ids)
    print("###############")

    
def get_number_tokens(newick, tokenizer):
    return len(tokenizer(newick, return_tensors="pt").input_ids[0])

########## PERFORMING INFERENCE ##########
def perform_inference_img_processor(augment, image_path, model, image_processor, tokenizer):
    image = cv2.imread(image_path)
    augmented_tensor = augment(image=image)["image"] # tensor with shape (C, H, W)
    # turn the tensor into numpy array, transpose from (C, H, W) to (H, W, C) 
    augmented_np_arr = np.transpose(np.array(augmented_tensor),(1,2,0)) # TODO: test if it works this way
    pixel_values = image_processor(augmented_np_arr, return_tensors="pt").pixel_values
    generated_ids = model.generate(pixel_values)
    generated_text = tokenizer.batch_decode(generated_ids, skip_special_tokens=True)[0]
    print(generated_text)
    
########## CUSTOM DATA COLLATOR ##########
def custom_collator(batch):
    """
    Returns the output of the default data collator but with interpolate_pos_encoding set to
    true in order to be able to use image resolutions that are not 224x224 (pretrained ViT default).
    See: https://discuss.huggingface.co/t/fine-tuning-vit-with-more-patches-higher-resolution/18731/3
    """
    default_collator_output = default_data_collator(batch)
    default_collator_output["interpolate_pos_encoding"] = True
    return default_collator_output

def main():
    # Finetuning
    # load a fine-tuned image captioning model, corresponding tokenizer, image processor and trainer
    model = VisionEncoderDecoderModel.from_pretrained("nlpconnect/vit-gpt2-image-captioning")
    pretrained_image_processor = ViTImageProcessor.from_pretrained("nlpconnect/vit-gpt2-image-captioning")
    tokenizer = AutoTokenizer.from_pretrained("nlpconnect/vit-gpt2-image-captioning") # TODO: gpt2 from huggingface instead?
    special_tokens = ["(", ")", ";", ",", ":"] # TODO: add special tokens or not?
    tokenizer.add_tokens(special_tokens, special_tokens=True)
    model.decoder.resize_token_embeddings(len(tokenizer))
    dataset_path = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\datasets\openai_test"
    r"""
    dataset = preprocess_dataset(
        dataset=load_dataset(dataset_path=dataset_path), 
        tokenizer=tokenizer, 
        image_processor=pretrained_image_processor
    )
    print(f"main(): Dataset was loaded and preprocessed.")
    # use set_transform to apply image augmentations on the fly during training instead of map() to save disk space
    # train in 2 stages: baseline and then robustness => augment_1_ar then augment_2 (with dropout etc.)
    dataset.set_transform(transform=transform_pad_resize)
    print(f"main(): Dataset transform set to transform_1_ar.")
    # split dataset into 80% training and 20% evaluation
    split_dataset = dataset.train_test_split(test_size=0.2)
    train_dataset = split_dataset["train"]
    eval_dataset = split_dataset["test"]
    print(f"Dataset split into train and split dataset.")
    # Debugging
    # print(100 * "#")
    # print("dataset: " + str(dataset))
    # print("Features:")
    # print(str(dataset.features))
    # print(100 * "#")
    # print("type of dataset: " + str(type(dataset)))
    # print(100 * "#")
    # print("First dataset entry:")
    # print("pixel_values: " + str(dataset[0]["pixel_values"]))
    # print("pixel_values type: " + str(type(dataset[0]["pixel_values"])))
    # print("labels: " + str(dataset[0]["labels"]))
    # print("labels type: " + str(type(dataset[0]["labels"])))
    # print(100 * "#")
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
    print(f"main(): Training arguments were set.")
    trainer = Trainer(
        model=model,
        args=training_args,
        train_dataset=train_dataset,
        data_collator=custom_collator,
        eval_dataset=eval_dataset,
        tokenizer=tokenizer,
        # compute_metrics=compute_metrics, TODO: add compute_metrics function
    )
    print(f"main(): Trainer was instantiated.")

    # Train the model
    train_results = trainer.train()
    print(f"main(): trainer.train() finished.")
    trainer.save_model()
    trainer.log_metrics("train", train_results.metrics)
    trainer.save_metrics("train", train_results.metrics)
    trainer.save_state()
    """
    
    
    test = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\datasets\openai_test\data190336566882\tree190336566882.jpg"
    augment = augm.pad_resize_augment_wrapper(to_tensor=True)
    perform_inference_img_processor(augment=augment,image_path=test, model=model, image_processor=pretrained_image_processor, tokenizer=tokenizer)
    # response: "a series of images showing a person standing in a room"
# execute main
main()



from transformers import BlipProcessor, BlipForQuestionAnswering, AutoConfig
from PIL import Image
import json, os
import torch
from extracting_phylogenies.hf_finetuning import finetuning_util as util 

def main():
    
    # test data
    test_data_dir = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\datasets\openai_test"

    # load dataset from dir
    dataset = util.load_dataset(test_data_dir)
    # get first 5 images for testing
    image_paths = dataset["image"][:5]
    
    # initialize HF blip
    config = AutoConfig.from_pretrained('Salesforce/blip-vqa-base', ignore_mismatched_sizes=True)
    model = BlipForQuestionAnswering._from_config(config)
    # model = model.cuda()
    processor  = BlipProcessor.from_pretrained('Salesforce/blip-vqa-base')
    
    # list for results
    results = []

    # question 
    question = "What does this image show?"
    question2 = "Give me all text labels in this image"
    # =>  ##gfrza isisile isis republished crashing mccormick reykjavik prairie ke isis republished crashing antwerp terraces pauline modifying hiscytes
    question3 = "How many lines are in this image?"
    for path in image_paths:
        
        image = Image.open(path).convert("RGB")

        # prepare inputs
        # encoding = processor(image, question, return_tensors="pt").to("cuda:0", torch.float16)
        encoding = processor(image, question3, return_tensors="pt")

        output = model.generate(**encoding)
        generated_text = processor.decode(output[0], skip_special_tokens=True)


        results.append(generated_text)

    print(results)
    
main()
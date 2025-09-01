from transformers import Qwen2VLForConditionalGeneration, Qwen2VLProcessor
from qwen_vl_utils import process_vision_info
import torch
from extracting_phylogenies.utilities import finetuning_util as util
from extracting_phylogenies.newick_extraction_openai import instructions as instr
from datasets import Image

# base model id 
model_id = "Qwen/Qwen2-VL-7B-Instruct"

# get base models
model = Qwen2VLForConditionalGeneration.from_pretrained(
    model_id,
    device_map=1,
    torch_dtype=torch.bfloat16,
)

processor = Qwen2VLProcessor.from_pretrained(model_id)


def generate_text_from_sample(model, processor, sample, max_new_tokens=1024):
    # Prepare the text input by applying the chat template
    text_input = processor.apply_chat_template(sample, tokenize=False, add_generation_prompt=True)
    image_inputs, _ = process_vision_info(sample)
    model_inputs = processor(
        text=[text_input],
        images=image_inputs,
        return_tensors="pt",
    ).to(model.device)
    generated_ids = model.generate(**model_inputs, min_new_tokens=64, max_new_tokens=max_new_tokens,)
    # Decode the output text
    output_text = processor.batch_decode(
        generated_ids, skip_special_tokens=True, clean_up_tokenization_spaces=False
    )
    return output_text[0]  

import re

def get_qwen2_output(sample, model, processor=processor):
    instructions = instr.instr_nwk_regular
    prompt = instr.prompt
    id = sample["id"]
    original_newick = sample["label"]
    sample_formatted = util.format_data_inference(sample["image"], instructions, prompt)
    generated_newick = re.search(
        r"(?<=assistant\n).+", 
        generate_text_from_sample(model=model, processor=processor, sample=sample_formatted)
    ).group()
    return id, original_newick, generated_newick


# load dataset
dataset_path = r"/home/henkelm/dataset_branchlengths/rand_dataset_branch_labels"
dataset = util.load_dataset(dataset_path, torch_dataset=True)
dataset = dataset.cast_column("image", Image(decode=True))

import os
import subprocess
from extracting_phylogenies.utilities import newick_util 

dir_path = "/home/henkelm/qwen2_base_branchlengths/"
tsv_path = "/home/henkelm/qwen2_branchlengths.tsv"


# create dirs for newick pairs 
for sample in dataset:
    id, original_newick, generated_newick = get_qwen2_output(sample=sample, model=model, processor=processor)
    print(id)
    
    # get the newick from the output
    if re.search(r"\(.+;", generated_newick):
        generated_newick = re.search(r"\(.+;", generated_newick).group()
    else:
        print("newick invalid")
        continue
    
    # postprocess 
    if not newick_util.is_balanced(generated_newick):
        print("Balancing out parentheses.")
        generated_newick = newick_util.balance_parentheses(generated_newick)
    # truncate
    index = generated_newick.find(";")
    generated_newick = generated_newick[:index+1]
    
    os.makedirs(os.path.join(dir_path, f"{id}"), exist_ok=True)
    original_newick_path = os.path.join(dir_path, f"{id}", f"original{id}.nwk")
    generated_newick_path = os.path.join(dir_path, f"{id}", f"base{id}.nwk")
    
    
    
    # skip invalid newicks
    if newick_util.is_newick(generated_newick):
        # remove duplicates
        generated_newick = newick_util.remove_duplicate_leaves(generated_newick)
        
    else:
        print("newick invalid")
        continue
    
    
    if not newick_util.is_newick(generated_newick):
        print("newick invalid")
        continue
    
    with open(original_newick_path, "w") as original:
        original.write(original_newick)
    with open(generated_newick_path, "w") as generated:
        generated.write(generated_newick)
        
    # get the params file
    params_path = os.path.join(dataset_path, f"data{id}", f"params{id}.tsv")
    
    # compare the newicks
    subprocess.run(
        [
            "python", 
            "-m", 
            "extracting_phylogenies.newick_comparison.newick_comparison", 
            "--original_newick", f"{original_newick_path}",
            "--generated_newick", f"{generated_newick_path}",
            "--outfile", f"{tsv_path}",
            "--params", f"{params_path}"
        ],
        check=True
    )

    

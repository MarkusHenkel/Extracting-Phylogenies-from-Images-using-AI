from transformers import Qwen2VLProcessor, Qwen2VLForConditionalGeneration
from PIL import Image
import torch
from extracting_phylogenies.newick_extraction_openai import instructions as instr
from extracting_phylogenies.utilities import finetuning_util as util

model_id = "Qwen/Qwen2-VL-2B-Instruct"

model = Qwen2VLForConditionalGeneration.from_pretrained(
    model_id,
    device_map=0,
    torch_dtype=torch.bfloat16,
)

processor = Qwen2VLProcessor.from_pretrained(model_id)

instructions = instr.instr_nwk_regular
prompt= instr.prompt

# load dataset
dataset_path = r"/home/henkelm/dataset_branchlengths"
dataset = util.load_dataset(dataset_path, torch_dataset=True)

image = dataset[0]["image"]

messages = [
    {
        "role": "system",
        "content": [
            {"type": "text", "text": instructions}
        ],
    },
    {
        "role": "user",
        "content": [
            {"type": "image", "image": image},
            {"type": "text", "text": prompt}
        ]
    },
]
inputs = processor.apply_chat_template(
	messages,
	add_generation_prompt=True,
	tokenize=True,
	return_dict=True,
	return_tensors="pt",
).to(model.device)

# TODO max_new_tokens higher?
outputs = model.generate(**inputs, max_new_tokens=512)
print(processor.decode(outputs[0][inputs["input_ids"].shape[-1]:]))


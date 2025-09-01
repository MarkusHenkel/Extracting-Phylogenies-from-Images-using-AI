import torch
from transformers import pipeline
from PIL import Image

pipeline = pipeline(
    task="document-question-answering",
    model="naver-clova-ix/donut-base-finetuned-docvqa",
    device=0,
    torch_dtype=torch.float16
)
image = Image.open(r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\datasets\openai_test\data190336566882\tree190336566882.jpg")

question1 = "List all text labels in this image seperated by commas." # returns one taxon
question2 = "What type of image is this?" # returns one taxon 
question3 = "Give me the Newick string corresponding to this phylogenetic tree." # returns one taxon
question4 = "Give me the topmost text label in this image" # returns one taxon 

pipeline_dict = pipeline(image=image, question=question4)
print(pipeline_dict)
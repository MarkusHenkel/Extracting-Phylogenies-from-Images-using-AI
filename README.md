# Extracting-Phylogenies-from-Images-using-AI
Repository for my Bachelors thesis on the extraction of phylogenies from images using AI namely the OpenAI API.
## General idea
### Pipeline
#### Data Generation part
1. Generate a random number x for the amount of taxa 
2. Generate x random numbers in the right range for the actual IDs of the taxa in the taxonomy
3. Generate the Newick tree using the ETE Toolkit
4. Save the Newick tree into a file so that Biopython can make an image out of the Newick tree using the draw function
5. Generate a random number for each edge in the image for randomizing the distances (Work in progress)
6. Use Regex? to put these random numbers at the right spots inside the newick tree inside the file
7. Draw the Newick tree using the draw function from Biopython
8. Save the drawn image using matplotlib
9. Put both image and text file containing newick tree into a zip file

#### API part
1. Take zip file as input 
2. take image out of zip file and give it to the API
3. let API generate a newick tree
4. put generated newick tree into a text file 
5. add text file to the zip file

#### comparison part
1. take zip folder as input 
2. take both text files and compare the strings 
3. analyse chatgpts mistakes my marking zip files with diverging newick trees

### TODOs

- [x] Taxa in the image are the IDs not the actual taxa they correspond to
- [x] parameter for preferred amount of taxa is not required and defaults to 10. How to default to 10 when amount_taxa is not specified??? 
- [x] How to save images generated? Save whatever image is generated in the opened window?
- [x] How should the image and the newick tree be saved? Put into a folder and then zipped? 
- [x] the save button of the window which shows the generated phylogenetic tree is stupid: the saved png will show exactly what is seen in the window. The window is not opened in fullscreen and will therefore be really cramped  just like the resulting png.
- [ ] catch error when specified amount of taxa is 0
- [x] make it so that if the specified maximum amount of taxa is 1 it is actually 1 and not 
- [ ] make randomizing the amount of taxa generated optional via a parameter
- [ ] create the folder, put the text file with newick and the image in it 
- [ ] create a bash file that executes the python script and takes argument like how many examples are supposed to be created, which type of example
### Problems
- [ ] phylo trees are not strictly binary
- [ ] taxa amount of 0 might cause issues
- [ ] issue with the taxa names in the tree, dont seem to be correct sometimes, maybe taxa names are truncated at spaces?
- [ ]  

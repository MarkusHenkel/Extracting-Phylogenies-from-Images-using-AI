from Bio import Phylo
from ete3 import NCBITaxa
from io import StringIO 
# ncbi taxonomy
ncbi = NCBITaxa()
# Random numbers for entries in the ncbi taxonomy
test = ncbi.get_topology([9222, 8123, 9023, 4002, 4728, 8473, 1923, 2130, 1023, 3493, 2139])
# write into newick format
test_newick = test.write()
print(test_newick)
# make the the newick tree (string) into a file so that it can be drawn by bio.phylo (not elegant?)
test_newick_file = StringIO(test_newick)
# draw the tree using bio.phylo and matplotlib 
test_newick_tree = Phylo.read(test_newick_file, "newick")
test_newick_tree.rooted = True
Phylo.draw(test_newick_tree)

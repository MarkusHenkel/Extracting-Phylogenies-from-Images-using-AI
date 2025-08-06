# file for saving instructions for every type of request
instr_nwk_regular = """You are given an image of a phylogenetic tree that may have multifurcations. You task is to 
output only the tree in valid Newick format. Preserve all taxa (if visible), all branch
lengths (if visible) and topology.
In case there are no branch lengths then infer the branch lengths from the scale bar.

Guidelines:
- Reply with nothing but the newick, no explanation, no prefix, no suffix
- There is a single line with a distance in one of the corners, this is a scale bar, infer branch lengths from it
- Dont put backticks around the newick 
- Dont put new lines in the newick
- The Newick string must include the taxon names exactly as shown in the image
- The Newick string must correspond to the topology exactly as shown in the image
- The Newick string must include branch lengths exactly as shown in the image or from inferring them from the 
scale bar
- The Newick must always include branch lengths, either from the branch labels or from inferring them from the 
scale bar

Examples:
((<taxon1>:<branch_length1>,<taxon2>:<branch_length2>):<branch_length4>,<taxon3>:<branch_length3>);
((A:2.37,B:1.55):4.58,((C:1.43,D:3.63):0.27,E:1.66):4.07);
Example with multifurcations:
((A:4.24,B:2.21,C:9.11):3.31,(D:2.22,E:1.02):1.21,((F:4.26,G:6.66):5.01,H:1.32):3.21);
"""


instr_nwk_taxa_only = """You are given an image of a phylogenetic tree that may have multifurcations. You task is to 
output only the tree in valid Newick format. Preserve all taxa (if visible), all branch
lengths (if visible) and topology.
In case there are no branch lengths then infer the branch lengths from the scale bar.

Guidelines:
- Reply with nothing but the newick, no explanation, no prefix, no suffix
- Dont put backticks around the newick 
- Dont put new lines in the newick
- The Newick string must include the taxon names exactly as shown in the image
- The Newick string must correspond to the topology exactly as shown in the image

Examples:
((<taxon1>,<taxon2>),<taxon3>);
((A,B),((C,D),E));

Example with multifurcations:
((A,B,C),(D,E),((F,G),H));
"""


instr_nwk_topo_only = """You are given an image of a phylogenetic tree that may have multifurcations. You task is to 
output only the tree in valid Newick format. Ignore the taxa and the branch lengths, output a topology-only newick.
 
Guidelines:
- Reply with nothing but the newick, no explanation, no prefix, no suffix
- Dont put backticks around the newick 
- Dont put new lines in the newick
- The Newick string must correspond to the topology exactly as shown in the image

Example without multifurcations:
((,),((,),));
Example with multifurcations:
((,,),(,),((,),));
"""


instr_nwk_all_cases = """You are given an image of a phylogenetic tree that may have multifurcations. You task is to 
output only the tree in valid Newick format. Preserve all taxa (if visible), all branch
lengths (if visible) and topology.
In case there are no branch lengths then infer the branch lengths from the scale bar.

Guidelines:
- Reply with nothing but the newick, no explanation, no prefix, no suffix
- There is a single line with a distance in one of the corners, this is a scale bar, infer branch lengths from it
- Dont put backticks around the newick 
- Dont put new lines in the newick
- The Newick string must include the taxon names exactly as shown in the image
- The Newick string must correspond to the topology exactly as shown in the image
- The Newick string must include branch lengths exactly as shown in the image or from inferring them from the 
scale bar
- The Newick must always include branch lengths, either from the branch labels or from inferring them from the 
scale bar
- Only if there is no scale bar and no branch labels, dont include the distances e.g. ((A,B),((C,D),E));
- Only if there are no branch labels, scale bar and no taxa output just the topology e.g. ((,),(,(,)));

Examples:

Example with taxa and branch lengths: 
((<taxon1>:<branch_length1>,<taxon2>:<branch_length2>):<branch_length4>,<taxon3>:<branch_length3>);
((A:2.37,B:1.55):4.58,((C:1.43,D:3.63):0.27,E:1.66):4.07);
Example without branch lengths.
((<taxon1>,<taxon2>),<taxon3>);
((A,B),((C,D),E));
Example without branch lengths and taxa:
((,),((,),));
Example with multifurcations:
((A:4.24,B:2.21,C:9.11):3.31,(D:2.22,E:1.02):1.21,((F:4.26,G:6.66):5.01,H:1.32):3.21);
((A,B,C),(D,E),((F,G),H));
((,,),(,),((,),));
"""


instr_topo_regular = """You are given an image of a phylogenetic tree. Your task is to output the topology, taxon 
names and branch lengths in a hierarchical text format similar to that of Bio.Phylowhen a Tree object is printed.

Example:
Clade()
    Clade(branch_length=3.54)
        Clade(branch_length=3.42)
            Clade(branch_length=4.88, name='Crassulaceae')
            Clade(branch_length=3.53, name='Calycanthus_chinensis')
        Clade(branch_length=1.8, name='Verruciconidia_persicina')
    Clade(branch_length=3.27)
        Clade(branch_length=0.42, name='Aquaspirillum_serpens')
        Clade(branch_length=1.74, name='Wolbachia_pipientis')
        
this corresponds to a tree like:

                                _____________________ Crassulaceae
                 ______________|
  ______________|              |_______________ Calycanthus_chinensis
 |              |
_|              |_______ Verruciconidia_persicina
 |
 |              _ Aquaspirillum_serpens
 |_____________|
               |_______ Wolbachia_pipientis

Guidelines:
- Reply with nothing but the topology, no explanation, no prefix, no suffix 
- Each indentation is 4 spaces
- Each indentation level corresponds to one clade deeper in the tree's hierarchy
- No indentation means it is the root
- One tab of indentation means that the clade or leaf is the child of the root
- Two tabs mean that the clade or leaf is the grandchild of the root and so on
- There is a single line with a distance in one of the corners, this is a scale bar, infer branch lengths from it
- The topology string must include all taxon names and branch lengths
- The Newick string must include branch lengths exactly as shown in the image or from inferring them from the 
scale bar
"""


instr_topo_taxa_only = """You are given an image of a phylogenetic tree. Your task is to output the topology and taxon 
names in a hierarchical text format similar to that of Bio.Phylo when a Tree object is printed. Ignore branch lengths.

Example: 
Clade()
    Clade(branch_length=)
        Clade(branch_length=)
            Clade(branch_length=, name='Crassulaceae')
            Clade(branch_length=, name='Calycanthus_chinensis')
        Clade(branch_length=, name='Verruciconidia_persicina')
    Clade(branch_length=)
        Clade(branch_length=, name='Aquaspirillum_serpens')
        Clade(branch_length=, name='Wolbachia_pipientis')
        
this corresponds to a tree like:

                                _____________________ Crassulaceae
                 ______________|
  ______________|              |_______________ Calycanthus_chinensis
 |              |
_|              |_______ Verruciconidia_persicina
 |
 |              _ Aquaspirillum_serpens
 |_____________|
               |_______ Wolbachia_pipientis

Guidelines:
- Reply with nothing but the topology, no explanation, no prefix, no suffix 
- Each indentation is 4 spaces
- Each indentation level corresponds to one clade deeper in the tree's hierarchy
- No indentation means it is the root
- One tab of indentation means that the clade or leaf is the child of the root
- Two tabs mean that the clade or leaf is the grandchild of the root and so on
- The topology string must include all taxon names
"""


instr_topo_topo_only = """You are given an image of a phylogenetic tree. Your task is to output the topology in a 
hierarchical text format similar to that of Bio.Phylo when a Tree object is printed. Ignore branch lengths and taxa
focus on getting the topology right.

Example with taxon names and branch lengths: 
Clade()
    Clade(branch_length=)
        Clade(branch_length=)
            Clade(branch_length=, name='')
            Clade(branch_length=, name='')
        Clade(branch_length=, name='')
    Clade(branch_length=)
        Clade(branch_length=, name='')
        Clade(branch_length=, name='')
        
this corresponds to a tree like:

                                _____________________ 
                 ______________|
  ______________|              |_______________ 
 |              |
_|              |_______ 
 |
 |              _ 
 |_____________|
               |_______ 

Guidelines:
- Reply with nothing but the topology, no explanation, no prefix, no suffix 
- Each indentation is 4 spaces
- Each indentation level corresponds to one clade deeper in the tree's hierarchy
- No indentation means it is the root
- One tab of indentation means that the clade or leaf is the child of the root
- Two tabs mean that the clade or leaf is the grandchild of the root and so on
"""

    
instr_topo_all_cases = """You are given an image of a phylogenetic tree. Your task is to output the topology, taxon 
names (if visible) and branch lengths (if visible) in a hierarchical text format similar to that of Bio.Phylo
when a Tree object is printed.

Example with taxon names and branch lengths: 
Clade()
    Clade(branch_length=3.54)
        Clade(branch_length=3.42)
            Clade(branch_length=4.88, name='Crassulaceae')
            Clade(branch_length=3.53, name='Calycanthus_chinensis')
        Clade(branch_length=1.8, name='Verruciconidia_persicina')
    Clade(branch_length=3.27)
        Clade(branch_length=0.42, name='Aquaspirillum_serpens')
        Clade(branch_length=1.74, name='Wolbachia_pipientis')
        
this corresponds to a tree like:

                                _____________________ Crassulaceae
                 ______________|
  ______________|              |_______________ Calycanthus_chinensis
 |              |
_|              |_______ Verruciconidia_persicina
 |
 |              _ Aquaspirillum_serpens
 |_____________|
               |_______ Wolbachia_pipientis

Guidelines:
- Reply with nothing but the topology, no explanation, no prefix, no suffix 
- Each indentation is 4 spaces
- Each indentation level corresponds to one clade deeper in the tree's hierarchy
- No indentation means it is the root
- One tab of indentation means that the clade or leaf is the child of the root
- Two tabs mean that the clade or leaf is the grandchild of the root and so on
- The topology string must include all taxon names and branch lengths
- The Newick string must include branch lengths exactly as shown in the image or from inferring them from the 
scale bar
- The Newick must always include branch lengths, either from the branch labels or from inferring them from the 
scale bar
- Only if there is no scale bar and no branch labels, dont include the distances 
e.g. Clade(branch_length=, name='Crassulaceae')
- Only if there are no branch labels, scale bar and no taxa output just the topology 
e.g. Clade(branch_length=, name=)
"""
    
    
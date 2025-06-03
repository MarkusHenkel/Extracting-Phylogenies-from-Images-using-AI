from data_collection import phylogenetic_trees_and_newicks
import argparse
import datetime

def main():
    argument_parser = argparse.ArgumentParser(
        description="""Comparison of an AI-generated Newick and the actual Newick.
        Modes:\n
        \tgenerate_image: Use the phylogenetic_trees_and_newicks script to generate """
    )
    argument_parser.add_argument("--mode", required=True, choices=["generate_image"], 
                                 help="Specify which functionality you want to use.")
    argument_parser.add_argument("-a", "--actual_newick", required=False, type=str,
                                 help="Path to the .nwk of the actual Newick.")
    argument_parser.add_argument("-g", "--generated_newick", required=False, type=str,
                                 help="Path to the .nwk of the AI-generated Newick.")
    argument_parser.add_argument("-o", "--outfile", required=False, type=str,
                                 help="Path where output is created at.")
    args = argument_parser.parse_args()
    path_actual_newick = args.actual_newick
    path_generated_newick = args.generated_newick
    path_outfile = args.outfile
    mode = args.mode
    if mode == "generate_image":
        if not path_actual_newick or not path_generated_newick or not path_outfile:
            exit(f"""[{datetime.datetime.now}] When using '--mode generate_image' please provide the path where the image 
                 is supposed to be saved at and either an actual newick or a generated newick.""")
        if path_generated_newick:
            with open(path_generated_newick, "r") as nwk_file:
                newick = nwk_file.read()
        # Generate the image of the phylogenetic tree that represents the AI generated Newick
        phylogenetic_trees_and_newicks.save_newick_image(newick=newick, outfile_path=path_outfile)

main()
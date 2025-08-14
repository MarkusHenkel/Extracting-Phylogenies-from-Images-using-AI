import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb

def main():
    path_branchlengths = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\results\rand_dataset_branch_labels_4.1_5.tsv"
    path_no_branchlenghts = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\results\rand_dataset_no_branch_labels_4.1_5.tsv"
    branchlength_dataframe = pd.read_csv(path_branchlengths, sep="\t", header=0)
    no_branchlength_dataframe = pd.read_csv(path_no_branchlenghts, sep="\t", header=0)
    merged_dataframe = pd.concat([branchlength_dataframe, no_branchlength_dataframe], ignore_index=True)
    
    ######### gpt5 and 4.1 rf dist ratio vs. taxa count
    merged_dataframe["newick2_grouped"] = merged_dataframe["newick2"].map(lambda x: "gpt-5" if "gpt-5" in x else "gpt-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=merged_dataframe, x="count_taxa1", y="rf_ratio", hue="newick2_grouped")
    plt.title("GPT-5 and GPT-4.1 RF Distance Distribution")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean RF Distance")
    # plt.show()
    
    
    ############ gpt 4.1 and 5, taxa count vs. correct edge ratio
    merged_dataframe["newick2_grouped"] = merged_dataframe["newick2"].map(lambda x: "gpt-5" if "gpt-5" in x else "gpt-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=merged_dataframe, x="count_taxa1", y="correct_edge_ratio", hue="newick2_grouped")
    plt.title("GPT-5 and GPT-4.1 Correct Edge Ratio Distribution")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean Correct Edge Ratio")
    # plt.show()
    
    ############ gpt 4.1 and 5, taxa count vs. correct taxa ratio
    merged_dataframe["newick2_grouped"] = merged_dataframe["newick2"].map(lambda x: "gpt-5" if "gpt-5" in x else "gpt-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=merged_dataframe, x="count_taxa1", y="correct_taxa_ratio", hue="newick2_grouped")
    plt.title("GPT-5 and GPT-4.1 Correct Taxa Ratio Distribution")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean Correct Taxa Ratio")
    # plt.show()
    
    ################## gpt 4.1 and 5, leaf parent branch lengths difference branch lengths displayed 
    branchlength_dataframe["newick2_grouped"] = branchlength_dataframe["newick2"].map(lambda x: "gpt-5" if "gpt-5" in x else "gpt-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=branchlength_dataframe, x="count_taxa1", y="mean_abs_diff_leaf_dists", hue="newick2_grouped")
    plt.title("GPT-4.1 and GPT-5 leaf-parent branch length difference on trees with branch labels")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean of absolute leaf-parent branch length difference")
    # plt.show()
    
    ################## gpt 4.1 and 5, leaf parent branch lengths difference no branch lengths displayed 
    no_branchlength_dataframe["newick2_grouped"] = no_branchlength_dataframe["newick2"].map(lambda x: "gpt-5" if "gpt-5" in x else "gpt-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=no_branchlength_dataframe, x="count_taxa1", y="mean_abs_diff_leaf_dists", hue="newick2_grouped")
    plt.title("GPT-4.1 and GPT-5 leaf-parent branch length difference on trees without branch labels")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean of absolute leaf-parent branch length difference")
    plt.show()
    
    ################## gpt 4.1 and 5, leaf leaf distance difference branch lengths displayed 
    branchlength_dataframe["newick2_grouped"] = branchlength_dataframe["newick2"].map(lambda x: "gpt-5" if "gpt-5" in x else "gpt-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=branchlength_dataframe, x="count_taxa1", y="mean_pairwise_dist_diff", hue="newick2_grouped")
    plt.title("GPT-4.1 and GPT-5 pairwise leaf-leaf distance on trees with branch labels")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean of mean pairwise leaf-leaf distance")
    plt.show()
    
    ################## gpt 4.1 and 5, leaf leaf distance difference no branch lengths displayed 
    no_branchlength_dataframe["newick2_grouped"] = no_branchlength_dataframe["newick2"].map(lambda x: "gpt-5" if "gpt-5" in x else "gpt-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=no_branchlength_dataframe, x="count_taxa1", y="mean_pairwise_dist_diff", hue="newick2_grouped")
    plt.title("GPT-4.1 and GPT-5 pairwise leaf-leaf distance difference on trees without branch labels")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean of mean pairwise leaf-leaf distance")
    plt.show()
main()
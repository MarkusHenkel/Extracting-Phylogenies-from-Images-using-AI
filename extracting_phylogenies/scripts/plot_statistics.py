import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb

def main():
    path_branchlengths = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\results_topo\finetuned_nwk_branchlengths.tsv"
    path_no_branchlenghts = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\results_topo\finetuned_nwk_no_branchlengths.tsv"
    
    # path_branchlengths = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\qwen2_branchlengths.tsv"
    # path_no_branchlenghts = r"C:\Users\marku\Desktop\StudiumVault\Semester6\Bachelorarbeit\Code\Extracting-Phylogenies-from-Images-using-AI\qwen2_no_branchlengths.tsv"
    
    branchlength_dataframe = pd.read_csv(path_branchlengths, sep="\t", header=0)
    no_branchlength_dataframe = pd.read_csv(path_no_branchlenghts, sep="\t", header=0)
    merged_dataframe = pd.concat([branchlength_dataframe, no_branchlength_dataframe], ignore_index=True)
    
    merged_dataframe["Model"] = merged_dataframe["newick2"].map(lambda x: "GPT-5" if "gpt-5" in x else "GPT-4.1-ft" if "gpt-4.1_finetuned" in x else "GPT-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=merged_dataframe, x="count_taxa1", y="rf_ratio", hue="Model", marker="o")
    plt.title("RF Distance Ratio Distribution of OpenAI models")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean RF Distance Ratio per tree")
    plt.savefig("RF Distance Ratio Distribution of OpenAI models")
    
    
    merged_dataframe["Model"] = merged_dataframe["newick2"].map(lambda x: "GPT-5" if "gpt-5" in x else "GPT-4.1-ft" if "gpt-4.1_finetuned" in x else "GPT-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=merged_dataframe, x="count_taxa1", y="correct_edge_ratio", hue="Model", marker="o")
    plt.title("Correct Edge Ratio Distribution of OpenAI models")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean Correct Edge Ratio per tree")
    plt.savefig("Correct Edge Ratio Distribution of OpenAI models")
    
    merged_dataframe["Model"] = merged_dataframe["newick2"].map(lambda x: "GPT-5" if "gpt-5" in x else "GPT-4.1-ft" if "gpt-4.1_finetuned" in x else "GPT-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=merged_dataframe, x="count_taxa1", y="correct_taxa_ratio", hue="Model", marker="o")
    plt.title("Correct Taxa Ratio Distribution of OpenAI models")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean Correct Taxa Ratio per tree")
    plt.savefig("Correct Taxa Ratio Distribution of OpenAI models")
    
    branchlength_dataframe["Model"] = branchlength_dataframe["newick2"].map(lambda x: "GPT-5" if "gpt-5" in x else "GPT-4.1-ft" if "gpt-4.1_finetuned" in x else "GPT-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=branchlength_dataframe, x="count_taxa1", y="mean_abs_diff_leaf_dists", hue="Model", marker="o")
    plt.title("Leaf-to-parent branch length difference on trees with branch labels of OpenAI models")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean Absolute Leaf-to-Parent Branch Length Difference per tree")
    plt.savefig("Leaf-to-parent branch length difference on trees with branch labels of OpenAI models")
    
    no_branchlength_dataframe["Model"] = no_branchlength_dataframe["newick2"].map(lambda x: "GPT-5" if "gpt-5" in x else "GPT-4.1-ft" if "gpt-4.1_finetuned" in x else "GPT-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=no_branchlength_dataframe, x="count_taxa1", y="mean_abs_diff_leaf_dists", hue="Model", marker="o")
    plt.title("Leaf-to-Parent Branch Length Difference on trees without branch labels of OpenAI models")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean Absolute Leaf-to-Parent Branch Length Difference per tree")
    plt.savefig("Leaf-to-Parent Branch Length Difference on trees without branch labels of OpenAI models")
    
    branchlength_dataframe["Model"] = branchlength_dataframe["newick2"].map(lambda x: "GPT-5" if "gpt-5" in x else "GPT-4.1-ft" if "gpt-4.1_finetuned" in x else "GPT-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=branchlength_dataframe, x="count_taxa1", y="mean_pairwise_dist_diff", hue="Model", marker="o")
    plt.title("Pairwise leaf-to-leaf distance on trees with branch labels of OpenAI models")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean of Mean Pairwise Leaf-to-Leaf Distance per tree")
    plt.savefig("Pairwise leaf-to-leaf distance on trees with branch labels of OpenAI models")
    
    no_branchlength_dataframe["Model"] = no_branchlength_dataframe["newick2"].map(lambda x: "GPT-5" if "gpt-5" in x else "GPT-4.1-ft" if "gpt-4.1_finetuned" in x else "GPT-4.1")
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=no_branchlength_dataframe, x="count_taxa1", y="mean_pairwise_dist_diff", hue="Model", marker="o")
    plt.title("Pairwise Leaf-to-Leaf Distance Difference on trees without branch labels of OpenAI models")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean of Mean Pairwise Leaf-to-Leaf Distance per tree")
    plt.savefig("Pairwise Leaf-to-Leaf Distance Difference on trees without branch labels of OpenAI models")
    
    # expected amount taxa - actual amount taxa
    merged_dataframe["Model"] = merged_dataframe["newick2"].map(lambda x: "GPT-5" if "gpt-5" in x else "GPT-4.1-ft" if "gpt-4.1_finetuned" in x else "GPT-4.1")
    merged_dataframe["difference_taxa"] = abs(merged_dataframe["count_taxa1"] - merged_dataframe["count_taxa2"] )
    plt.figure(figsize=(10, 6))
    sb.lineplot(data=merged_dataframe, x="count_taxa1", y="difference_taxa", hue="Model", marker="o")
    plt.title("Absolute Difference of expected and actual Amount of Taxa of OpenAI models")
    plt.xlabel("Count of Taxa")
    plt.ylabel("Mean Absolute Taxa Count Difference per tree")
    plt.savefig("Absolute Difference of expected and actual Amount of Taxa of OpenAI models")
    # plt.show()
    
    
    
    
    
    
    
    
    
    # merged_dataframe["Model"] = merged_dataframe["newick2"].map(lambda x: "Qwen2-VL-7B-ft" if "finetuned" in x else "Qwen-VL-7B")
    # plt.figure(figsize=(10, 6))
    # sb.lineplot(data=merged_dataframe, x="count_taxa1", y="rf_ratio", hue="Model", marker="o")
    # plt.title("RF Distance Ratio Distribution of base and finetuned Qwen2")
    # plt.xlabel("Count of Taxa")
    # plt.ylabel("Mean RF Distance Ratio per tree")
    # plt.savefig("RF Distance Ratio Distribution of base and finetuned Qwen2")
    # # plt.show()
    
    # merged_dataframe["Model"] = merged_dataframe["newick2"].map(lambda x: "Qwen2-VL-7B-ft" if "finetuned" in x else "Qwen-VL-7B")
    # plt.figure(figsize=(10, 6))
    # sb.lineplot(data=merged_dataframe, x="count_taxa1", y="correct_edge_ratio", hue="Model", marker="o")
    # plt.title("Correct Edge Ratio Distribution of base and finetuned Qwen2")
    # plt.xlabel("Count of Taxa")
    # plt.ylabel("Mean Correct Edge Ratio per tree")
    # plt.savefig("Correct Edge Ratio Distribution of base and finetuned Qwen2")
    # # plt.show()
    
    # ############ gpt 4.1 and 5, taxa count vs. correct taxa ratio
    # merged_dataframe["Model"] = merged_dataframe["newick2"].map(lambda x: "Qwen2-VL-7B-ft" if "finetuned" in x else "Qwen-VL-7B")
    # plt.figure(figsize=(10, 6))
    # sb.lineplot(data=merged_dataframe, x="count_taxa1", y="correct_taxa_ratio", hue="Model", marker="o")
    # plt.title("Correct Taxa Ratio Distribution of base and finetuned Qwen2")
    # plt.xlabel("Count of Taxa")
    # plt.ylabel("Mean Correct Taxa Ratio per tree")
    # plt.savefig("Correct Taxa Ratio Distribution of base and finetuned Qwen2")
    # # plt.show()
    
    # ################## gpt 4.1 and 5, leaf parent branch lengths difference branch lengths displayed 
    # branchlength_dataframe["Model"] = branchlength_dataframe["newick2"].map(lambda x: "Qwen2-VL-7B-ft" if "finetuned" in x else "Qwen-VL-7B")
    # plt.figure(figsize=(10, 6))
    # sb.lineplot(data=branchlength_dataframe, x="count_taxa1", y="mean_abs_diff_leaf_dists", hue="Model", marker="o")
    # plt.title("Leaf-to-parent branch length difference on trees with branch labels of base and finetuned Qwen2")
    # plt.xlabel("Count of Taxa")
    # plt.ylabel("Mean Absolute Leaf-to-Parent Branch Length Difference per tree")
    # plt.savefig("Leaf-to-parent branch length difference on trees with branch labels of base and finetuned Qwen2")
    # # plt.show()

    # ################## gpt 4.1 and 5, leaf parent branch lengths difference no branch lengths displayed 
    # no_branchlength_dataframe["Model"] = no_branchlength_dataframe["newick2"].map(lambda x: "Qwen2-VL-7B-ft" if "finetuned" in x else "Qwen-VL-7B")
    # plt.figure(figsize=(10, 6))
    # sb.lineplot(data=no_branchlength_dataframe, x="count_taxa1", y="mean_abs_diff_leaf_dists", hue="Model", marker="o")
    # plt.title("Leaf-to-Parent Branch Length Difference on trees without branch labels of base and finetuned Qwen2")
    # plt.xlabel("Count of Taxa")
    # plt.ylabel("Mean Absolute Leaf-to-Parent Branch Length Difference per tree")
    # plt.savefig("Leaf-to-Parent Branch Length Difference on trees without branch labels of base and finetuned Qwen2")
    # # plt.show()

    # ################## gpt 4.1 and 5, leaf leaf distance difference branch lengths displayed 
    # branchlength_dataframe["Model"] = branchlength_dataframe["newick2"].map(lambda x: "Qwen2-VL-7B-ft" if "finetuned" in x else "Qwen-VL-7B")
    # plt.figure(figsize=(10, 6))
    # sb.lineplot(data=branchlength_dataframe, x="count_taxa1", y="mean_pairwise_dist_diff", hue="Model", marker="o")
    # plt.title("Pairwise leaf-to-leaf distance on trees with branch labels of base and finetuned Qwen2")
    # plt.xlabel("Count of Taxa")
    # plt.ylabel("Mean of ean Pairwise Leaf-to-Leaf Distance per tree")
    # plt.savefig("Pairwise leaf-to-leaf distance on trees with branch labels of base and finetuned Qwen2")
    # # plt.show()

    # ################## gpt 4.1 and 5, leaf leaf distance difference no branch lengths displayed 
    # no_branchlength_dataframe["Model"] = no_branchlength_dataframe["newick2"].map(lambda x: "Qwen2-VL-7B-ft" if "finetuned" in x else "Qwen-VL-7B")
    # plt.figure(figsize=(10, 6))
    # sb.lineplot(data=no_branchlength_dataframe, x="count_taxa1", y="mean_pairwise_dist_diff", hue="Model", marker="o")
    # plt.title("Pairwise Leaf-to-Leaf Distance Difference on trees without branch labels of base and finetuned Qwen2")
    # plt.xlabel("Count of Taxa")
    # plt.ylabel("Mean of Mean Pairwise Leaf-to-Leaf Distance per tree")
    # plt.savefig("Pairwise Leaf-to-Leaf Distance Difference on trees without branch labels of base and finetuned Qwen2")
    # # plt.show()
    
    # # expected amount taxa - actual amount taxa
    # merged_dataframe["Model"] = merged_dataframe["newick2"].map(lambda x: "Qwen2-VL-7B-ft" if "finetuned" in x else "Qwen-VL-7B")
    # merged_dataframe["difference_taxa"] = abs(merged_dataframe["count_taxa1"] - merged_dataframe["count_taxa2"] )
    # plt.figure(figsize=(10, 6))
    # sb.lineplot(data=merged_dataframe, x="count_taxa1", y="difference_taxa", hue="Model", marker="o")
    # plt.title("Absolute Difference of expected and actual Amount of Taxa of base and finetuned Qwen2")
    # plt.xlabel("Count of Taxa")
    # plt.ylabel("Mean Absolute Taxa Count Difference per tree")
    # plt.savefig("Absolute Difference of expected and actual Amount of Taxa of base and finetuned Qwen2")
    # # plt.show()
    
main()
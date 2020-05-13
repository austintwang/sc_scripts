import numpy as np
import os
import pickle
import pandas as pd

def read_data(plasma_data, coloc_data, clusters, gene_name):
    # print(coloc_data) ####
    data = []
    if not "clusters" in coloc_data:
        return data
    for c in clusters:
        plasma_clust = plasma_data.get(c, None)
        coloc_clust = coloc_data["clusters"].get(c, None)
        if plasma_clust is None or coloc_clust is None:
            continue
        # print(plasma_clust.keys()) ####
        # print(coloc_clust.keys()) ####
        data_clust = [
            gene_name, 
            c, 
            np.mean(plasma_clust.get("causal_set_indep", np.nan)), 
            np.mean(plasma_clust.get("causal_set_eqtl", np.nan)),
            coloc_clust.get("h4_indep_eqtl"),
            coloc_clust.get("h4_ase_eqtl"),
            coloc_clust.get("h4_eqtl_eqtl")
        ]
        # print(data_clust) ####
        data.append(data_clust)
    return data

def load_clusters(cluster_map_path):
    with open(cluster_map_path, "rb") as cluster_map_file:
        cluster_map = pickle.load(cluster_map_file)
    return cluster_map.keys()

def make_violin(
        df,
        var, 
        model_flavors,
        model_names, 
        model_colors,
        title, 
        result_path,
    ):
    sns.set(style="whitegrid", font="Roboto")
    plt.figure(figsize=(4,2))

    palette = [model_colors[m] for m in model_flavors]
    names = [model_names[m] for m in model_flavors]
    chart = sns.violinplot(
        x=var, 
        y="Model", 
        data=df, 
        order=model_flavors, 
        palette=palette,
        cut=0,
        scale="width"
    )
    ax = plt.gca()
    for art in ax.get_children():
        if isinstance(art, matplotlib.collections.PolyCollection):
            art.set_edgecolor((0., 0., 0.))
    plt.xlim(0., 1.)
    chart.set_yticklabels([model_names[m] for m in model_flavors])
    plt.ylabel("")
    plt.title(title)
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()

def plot_sets(df, out_dir):
    clusters = {
        "_all": "All Cells",
        "Ex": "Excitatory Neuron",
        "Oligo": "Oligodendrocyte",
        "Astro": "Astroglia",
        "In": "Inhibitory Neuron",
        "Endo": "Endothelial",
        "OPC": "Oligodendrocyte Progenitor",
        "Per": "Per"
    }
    model_map = {
        "PP4Joint": "PLASMA/C-J",
        "PP4AS": "PLASMA/C-AS",
        "PP4QTL": "QTL-Only"
    }
    var_dists = "PP4 Score"
    model_flavors = model_map.keys()
    pal = sns.color_palette()
    model_colors = {
        "PP4Joint": pal[0],
        "PP4AS": pal[4],
        "PP4QTL": pal[7],
    }
    for cluster in clusters.keys():
        df_clust = df.loc[df["Cluster"] == cluster]
        df_dists = pd.melt(
            df.loc[df["Cluster"] == cluster], 
            id_vars=["Gene"], 
            value_vars=model_map.keys(),
            var_name="Model",
            value_name=var_dists
        )
        title = clusters[cluster]
        make_violin(
            df_dists,
            var_dists, 
            model_flavors,
            model_map, 
            model_colors,
            title, 
            os.path.join(out_dir, "pp4s_{0}.svg".format(cluster)),
        )

def interpret_genes(genes_dir, gwas_name, cluster_map_path, out_dir):
    clusters = load_clusters(cluster_map_path)
    genes = os.listdir(genes_dir)
    data_lst = []
    for g in genes:
        gene_dir = os.path.join(genes_dir, g)
        plasma_path = os.path.join(gene_dir, "plasma.pickle")
        if os.path.isdir(plasma_path):
            plasma_path = os.path.join(plasma_path, "output.pickle")
        coloc_path = os.path.join(gene_dir, "coloc_{0}.pickle".format(gwas_name))
        if os.path.isdir(coloc_path):
            coloc_path = os.path.join(coloc_path, "output.pickle")
        try:
            with open(plasma_path, "rb") as plasma_file:
                plasma_data = pickle.load(plasma_file)
            with open(coloc_path, "rb") as coloc_file:
                coloc_data = pickle.load(coloc_file)
        except FileNotFoundError:
            continue

        data = read_data(plasma_data, coloc_data, clusters, g)
        data_lst.extend(data)

    cols = [
        "Gene", 
        "Cluster", 
        "CredibleSetPropIndep",
        "CredibleSetPropGWAS",
        "PP4Joint",
        "PP4AS",
        "PP4QTL"
    ]
    data_df = pd.DataFrame(data_lst, columns=cols)
    data_df.sort_values(by=["PP4Joint"], ascending=False, inplace=True)
    out_dir_gwas = os.path.join(out_dir, gwas_name)
    os.makedirs(out_dir_gwas, exist_ok=True)
    data_df.to_csv(os.path.join(out_dir_gwas, "data.csv"), index=False)
    with open(os.path.join(out_dir_gwas, "data.txt"), "w") as txt_file:
        data_df.to_string(txt_file)
    plot_sets(data_df, out_dir_gwas)


if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    # cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map.pickle")
    # genes_dir_kellis = os.path.join(data_path_kellis, "genes")

    # out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis/alz.txt"

    # interpret_genes(genes_dir_kellis, "alz", cluster_map_path_kellis, out_path_kellis)

    cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")
    out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/colocalization"
    gwass = [
        "AlzheimersMaternal",
        "AlzheimersPaternal",
        "AlzheimersProxyMetaIGAP",
        "BD",
        "BDSCZ",
        "CD",
        "DepressedAffect",
        "Depression",
        "IBD",
        "Intelligence",
        "MDD",
        "Neuroticism",
        "ReactionTime",
        "SCZ",
        "SCZvsBD",
        "UC",
        "VerbalNumericReasoning",
        "Worry"
    ]
    for g in gwass:
        interpret_genes(genes_dir_kellis, g, cluster_map_path_kellis, out_path_kellis)

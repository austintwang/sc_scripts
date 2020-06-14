#!/usr/bin/env python3

import sys
import numpy as np
import scipy.stats
import os
import pickle
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['agg.path.chunksize'] = 10000
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def read_data(plasma_data, coloc_data, clusters, gene_id, gene_name):
    # print(coloc_data.keys()) ####
    data = []
    if not "clusters" in coloc_data:
        return data
    top_z = np.nanmax(np.abs(coloc_data["z_beta"]))
    top_nlp = np.nan_to_num(-np.log10(scipy.stats.norm.sf(abs(top_z))*2))
    num_informative = np.sum(coloc_data.get("informative_snps", np.nan))
    for c in clusters:
        plasma_clust = plasma_data.get(c, None)
        coloc_clust = coloc_data["clusters"].get(c, None)
        if plasma_clust is None or coloc_clust is None:
            continue
        # print(plasma_clust.keys()) ####
        # print(coloc_clust.keys()) ####
        data_clust = [
            gene_id,
            gene_name, 
            c, 
            np.mean(plasma_clust.get("causal_set_indep", np.nan)), 
            np.mean(coloc_data.get("causal_set_eqtl", np.nan)),
            top_nlp,
            coloc_clust.get("h4_indep_eqtl"),
            coloc_clust.get("h4_ase_eqtl"),
            coloc_clust.get("h4_eqtl_eqtl"),
            num_informative
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
    plt.xlim(0, 1)
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
        df_dists = pd.melt(
            df, 
            # df.loc[np.logical_and(df["Cluster"] == cluster, df["GWASSig"] >= -np.log10(1))], 
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

def calc_sumstats(df, out_dir, thresh):
    df_coloc = df.loc[
        df["PP4Joint"] >= thresh
    ]
    df_ncoloc = df.loc[
        df["PP4Joint"] < thresh
    ]
    clusters = {
        "_all": "All Cells",
        "Ex": "Excitatory Neuron",
        "In": "Inhibitory Neuron",
        "Oligo": "Oligodendrocyte",
        "OPC": "Oligodendrocyte Progenitor",
        "Astro": "Astroglia",
        # "Endo": "Endothelial",
        # "Per": "Per"
    }
    coloc_data = {}
    dfs_clust = {}
    for i in clusters.keys():
        df_clust = df_coloc.loc[df_coloc["Cluster"] == i]
        coloc_data[i] = df_clust.count()["PP4Joint"]
        dfs_clust[i] = df_clust

    df_nall = df_ncoloc.loc[df_ncoloc["Cluster"] == "_all"]
    diff_data = {}
    for i in clusters.keys():
        df_diff = pd.merge(
            dfs_clust[i], 
            df_nall, 
            on=["Gene"], 
            suffixes=["_clust", "_all"]
        )
        diff_data[i] = df_diff.count()["PP4Joint_clust"]

    outs = ["Cluster\tNumColoc\tNumColocDiff"]
    for i in clusters.keys():
        outs.append(f"{i}\t{coloc_data[i]}\t{diff_data[i]}\n")

    with open(os.path.join(out_dir, "sumstats.txt"), "w") as out_file:
        out_file.writelines(outs)

def interpret_genes(genes_dir, genes_map_dir, gwas_name, cluster_map_path, out_dir, status_path):
    with open(status_path, "w") as status_file:
        status_file.write("")

    with open(genes_map_dir, "rb") as genes_map_file:
        genes_map = pickle.load(genes_map_file)

    clusters = load_clusters(cluster_map_path)
    genes = os.listdir(genes_dir)
    # genes = genes[:500] ####
    data_lst = []
    for g in genes:
        gene_dir = os.path.join(genes_dir, g)
        plasma_path = os.path.join(gene_dir, "combined", "plasma_i0.pickle")
        coloc_path = os.path.join(gene_dir, "coloc", f"{gwas_name}.pickle")
        try:
            with open(plasma_path, "rb") as plasma_file:
                plasma_data = pickle.load(plasma_file)
            with open(coloc_path, "rb") as coloc_file:
                coloc_data = pickle.load(coloc_file)
        except (FileNotFoundError, pickle.UnpicklingError):
            continue

        data = read_data(plasma_data, coloc_data, clusters, g, genes_map.get(g.split(".")[0]))
        data_lst.extend(data)

    cols = [
        "Gene", 
        "GeneName",
        "Cluster", 
        "CredibleSetPropIndep",
        "CredibleSetPropGWAS",
        "GWASSig",
        "PP4Joint",
        "PP4AS",
        "PP4QTL",
        "UsableSNPCount",
    ]
    data_df = pd.DataFrame(data_lst, columns=cols)
    data_df.sort_values(by=["PP4Joint"], ascending=False, inplace=True)
    out_dir_gwas = os.path.join(out_dir, gwas_name)
    os.makedirs(out_dir_gwas, exist_ok=True)
    data_df.to_csv(os.path.join(out_dir_gwas, "data.csv"), index=False)
    with open(os.path.join(out_dir_gwas, "data.txt"), "w") as txt_file:
        data_df.to_string(txt_file)
    calc_sumstats(data_df, out_dir_gwas, 0.1)
    plot_sets(data_df, out_dir_gwas)

    data_sig_df = data_df.loc[data_df["GWASSig"] >= -np.log10(5e-8)]
    out_dir_sig_gwas = os.path.join(out_dir, f"{gwas_name}_sig")
    os.makedirs(out_dir_sig_gwas, exist_ok=True)
    data_sig_df.to_csv(os.path.join(out_dir_sig_gwas, "data.csv"), index=False)
    with open(os.path.join(out_dir_sig_gwas, "data.txt"), "w") as txt_file:
        data_sig_df.to_string(txt_file)
    calc_sumstats(data_sig_df, out_dir_sig_gwas, 0.1)
    plot_sets(data_sig_df, out_dir_sig_gwas)

    with open(status_path, "w") as status_file:
        status_file.write("Complete")

if __name__ == '__main__':
    interpret_genes(*sys.argv[1:])

    # data_path_kellis = "/agusevlab/awang/sc_kellis"
    # # cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map.pickle")
    # # genes_dir_kellis = os.path.join(data_path_kellis, "genes")

    # # out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis/alz.txt"

    # # interpret_genes(genes_dir_kellis, "alz", cluster_map_path_kellis, out_path_kellis)

    # cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    # genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")
    # out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/colocalization"
    # gwass = [
    #     "AlzheimersMaternal",
    #     "AlzheimersPaternal",
    #     "AlzheimersProxyMetaIGAP",
    #     "BD",
    #     "BDSCZ",
    #     "CD",
    #     "DepressedAffect",
    #     "Depression",
    #     "IBD",
    #     "Intelligence",
    #     "MDD",
    #     "Neuroticism",
    #     "ReactionTime",
    #     "SCZ",
    #     "SCZvsBD",
    #     "UC",
    #     "VerbalNumericReasoning",
    #     "Worry"
    # ]
    # for g in gwass:
    #     print(g) ####
    #     interpret_genes(genes_dir_kellis, g, cluster_map_path_kellis, out_path_kellis)

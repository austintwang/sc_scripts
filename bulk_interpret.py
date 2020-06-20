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

def read_data_bulk(plasma_data, bulk_data, clusters, gene_name):
    # print(bulk_data.get("run_error")) ####
    # print(bulk_data.get("traceback")) ####
    data = []
    if not "clusters" in bulk_data:
        return data
    # top_snp = np.nanargmax(plasma_clust["ppas_indep"])
    # top_z = np.nanmax(np.abs(bulk_data["z_beta"]))
    # top_nlp = np.nan_to_num(-np.log10(scipy.stats.norm.sf(abs(top_z))*2))
    num_informative_bulk = np.sum(bulk_data.get("informative_snps", np.nan))
    for c in clusters:
        plasma_clust = plasma_data.get(c, None)
        bulk_clust = bulk_data["clusters"].get(c, None)
        if plasma_clust is None or bulk_clust is None:
            continue
        if "ppas_indep" not in plasma_clust:
            # print(plasma_clust) ####
            continue
        # print(plasma_clust.keys()) ####
        # print(coloc_clust.keys()) ####
        num_informative_plasma = plasma_clust["num_snps_informative"]
        try:
            top_snp = np.nanargmax(plasma_clust["ppas_indep"])
        except ValueError:
            continue
        top_z_beta = plasma_clust["z_beta"][top_snp]
        top_nlp_beta = -scipy.stats.norm.logsf(np.abs(top_z_beta)) / np.log(10) - np.log10(2)
        top_z_phi = plasma_clust["z_phi"][top_snp]
        top_nlp_phi = -scipy.stats.norm.logsf(np.abs(top_z_phi)) / np.log(10) - np.log10(2)
        top_z_bulk = bulk_data["z_beta"][top_snp]
        top_nlp_bulk = -scipy.stats.norm.logsf(np.abs(top_z_bulk)) / np.log(10) - np.log10(2)
        data_clust = [
            gene_name, 
            c, 
            np.mean(plasma_clust.get("causal_set_indep", np.nan)), 
            np.mean(bulk_data.get("causal_set_eqtl", np.nan)),
            top_z_bulk,
            top_nlp_bulk,
            top_z_phi,
            top_nlp_phi,
            top_z_beta,
            top_nlp_beta,
            bulk_clust.get("h0_indep_eqtl"),
            bulk_clust.get("h0_ase_eqtl"),
            bulk_clust.get("h0_eqtl_eqtl"),
            bulk_clust.get("h0_fmb_fmb"),
            bulk_clust.get("h1_indep_eqtl"),
            bulk_clust.get("h1_ase_eqtl"),
            bulk_clust.get("h1_eqtl_eqtl"),
            bulk_clust.get("h1_fmb_fmb"),
            bulk_clust.get("h2_indep_eqtl"),
            bulk_clust.get("h2_ase_eqtl"),
            bulk_clust.get("h2_eqtl_eqtl"),
            bulk_clust.get("h2_fmb_fmb"),
            bulk_clust.get("h3_indep_eqtl"),
            bulk_clust.get("h3_ase_eqtl"),
            bulk_clust.get("h3_eqtl_eqtl"),
            bulk_clust.get("h3_fmb_fmb"),
            bulk_clust.get("h4_indep_eqtl"),
            bulk_clust.get("h4_ase_eqtl"),
            bulk_clust.get("h4_eqtl_eqtl"),
            bulk_clust.get("h4_fmb_fmb"),
            num_informative_bulk,
            num_informative_plasma,
        ]
        # print(data_clust) ####
        data.append(data_clust)
    return data

def load_clusters(cluster_map_path):
    with open(cluster_map_path, "rb") as cluster_map_file:
        cluster_map = pickle.load(cluster_map_file)
    return cluster_map.keys()

def make_df_bulk(run_name, bulk_name, genes_dir, cluster_map_path):
    clusters = load_clusters(cluster_map_path)
    genes = os.listdir(genes_dir)
    genes = genes[:500] ####
    data_lst = []
    for g in genes:
        gene_dir = os.path.join(genes_dir, g)
        plasma_path = os.path.join(gene_dir, run_name, "plasma_i0.pickle")
        coloc_path = os.path.join(gene_dir, "bulk_qtl", f"{bulk_name}_out.pickle")
        try:
            with open(plasma_path, "rb") as plasma_file:
                plasma_data = pickle.load(plasma_file)
            with open(coloc_path, "rb") as coloc_file:
                coloc_data = pickle.load(coloc_file)
        except (FileNotFoundError, pickle.UnpicklingError):
            continue

        data = read_data_bulk(plasma_data, coloc_data, clusters, g)
        # print(data) ####
        data_lst.extend(data)
    # print(data_lst) ####

    cols = [
        "Gene", 
        "Cluster", 
        "CredibleSetPropIndep",
        "CredibleSetPropGWAS",
        "TopSNPZBulk",
        "TopSNPNLPBulk",
        "TopSNPZPhi",
        "TopSNPNLPPhi",
        "TopSNPZBeta",
        "TopSNPNLPBeta",
        "PP0Joint",
        "PP0AS",
        "PP0QTL",
        "PP0FINEMAP",
        "PP1Joint",
        "PP1AS",
        "PP1QTL",
        "PP1FINEMAP",
        "PP2Joint",
        "PP2AS",
        "PP2QTL",
        "PP2FINEMAP",
        "PP3Joint",
        "PP3AS",
        "PP3QTL",
        "PP3FINEMAP",
        "PP4Joint",
        "PP4AS",
        "PP4QTL",
        "PP4FINEMAP",
        "NumInformativeBulk",
        "NumInformativePlasma"
    ]
    data_df = pd.DataFrame(data_lst, columns=cols)
    data_df.sort_values(by=["PP4Joint"], ascending=False, inplace=True)
    return data_df

def make_heatmap(arr, order, title, result_path):
    heat_data = pd.DataFrame(data=arr, index=order, columns=["Bulk"])
    print(heat_data) ####
    sns.set(style="whitegrid", font="Roboto")
    f, ax = plt.subplots(figsize=(5, 5))
    sns.heatmap(heat_data, annot=True, fmt=".2g", square=True, cbar=False, annot_kws={"size": 10})
    plt.title(title)
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()

def plot_xcells(df, out_dir, stat_name):
    sn1 = stat_name
    df_tr_sig = df.loc[
        df[f"TopSNPNLP{sn1}"] >= -np.log10(0.05/df["NumInformativePlasma"])
    ]
    # print(df_tr_sig) ####
    # print(df) ####

    clusters = {
        "_all": "All Cells",
        "Ex": "Excitatory Neuron",
        "In": "Inhibitory Neuron",
        "Oligo": "Oligodendrocyte",
        "OPC": "Oligodendrocyte Progenitor",
        "Astro": "Astroglia",
        "Microglia": "Microglia",
        # "Endo": "Endothelial",
        # "Per": "Per"
    }
    cluster_order = list(clusters.keys())
    storey_pis = np.zeros((len(cluster_order), 1,),)
    for ind_i, i in enumerate(cluster_order):
        df_merged = df_tr_sig.loc[df_tr_sig["Cluster"] == i]
        # print(df_merged["TopSNPNLPBulk"].dtype) ####
        num_sig_train = np.count_nonzero(~np.isnan(pd.to_numeric(df_merged["TopSNPNLPBulk"])))
        num_sig_test = np.sum(df_merged["TopSNPNLPBulk"] >= -np.log10(0.05))
        print(num_sig_train, num_sig_test)
        storey_pis[ind_i, 0] = num_sig_test / num_sig_train


    title = "Bulk Replication Storey Pi"
    make_heatmap(storey_pis, cluster_order, title, os.path.join(out_dir, f"storey_pi_{sn1}.svg"))

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

def plot_sets(df, out_dir, hyp):
    h = hyp
    clusters = {
        "_all": "All Cells",
        "Ex": "Excitatory Neuron",
        "Oligo": "Oligodendrocyte",
        "Astro": "Astroglia",
        "In": "Inhibitory Neuron",
        "Endo": "Endothelial",
        "Microglia": "Microglia",
        "OPC": "Oligodendrocyte Progenitor",
        "Per": "Per"
    }
    model_map = {
        f"PP{h}Joint": "PLASMA/C-J",
        f"PP{h}AS": "PLASMA/C-AS",
        f"PP{h}QTL": "QTL-Only",
        f"PP{h}FINEMAP": "FINEMAP"
    }
    var_dists = "PP4 Score"
    model_flavors = model_map.keys()
    pal = sns.color_palette()
    model_colors = {
        f"PP{h}Joint": pal[0],
        f"PP{h}AS": pal[4],
        f"PP{h}QTL": pal[7],
        f"PP{h}FINEMAP": pal[3],
    }
    for cluster in clusters.keys():
        df_dists = pd.melt(
            df.loc[df["Cluster"] == cluster], 
            # df.loc[np.logical_and(df["Cluster"] == cluster, df["GWASSig"] >= -np.log10(1))], 
            id_vars=["Gene"], 
            value_vars=model_map.keys(),
            var_name="Model",
            value_name=var_dists
        )
        # df_dists[f"PP{h}Joint"] = np.nan_to_num(df_dists[f"PP{h}Joint"])
        # df_dists[f"PP{h}AS"] = np.nan_to_num(df_dists[f"PP{h}AS"])
        # df_dists[f"PP{h}QTL"] = np.nan_to_num(df_dists[f"PP{h}QTL"])
        # df_dists[f"PP{h}FINEMAP"] = np.nan_to_num(df_dists[f"PP{h}FINEMAP"])
        title = clusters[cluster]
        print(df_dists) ####
        make_violin(
            df_dists,
            var_dists, 
            model_flavors,
            model_map, 
            model_colors,
            title, 
            os.path.join(out_dir, f"pp{h}s_{cluster}.svg"),
        )

def get_info_xval(run_name, bulk_name, genes_dir, cluster_map_path, out_dir_base):
    out_dir = os.path.join(out_dir_base, bulk_name)
    df = make_df_bulk(run_name, bulk_name, genes_dir, cluster_map_path)
    plot_xcells(df, out_dir, "Phi")
    plot_xcells(df, out_dir, "Beta")
    plot_sets(df, out_dir, "0")
    plot_sets(df, out_dir, "1")
    plot_sets(df, out_dir, "2")
    plot_sets(df, out_dir, "3")
    plot_sets(df, out_dir, "4")
    csv_path = os.path.join(out_dir, "cluster_info.csv")
    df.to_csv(csv_path, index=False, na_rep="None")
    txt_path = os.path.join(out_dir, "cluster_info.txt")
    with open(txt_path, "w") as txt_file:
        df.to_string(txt_file)

if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"

    cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")

    out_dir_base_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429"

    get_info_xval("combined", "rosmap", genes_dir_kellis, cluster_map_path_kellis, out_dir_base_kellis)



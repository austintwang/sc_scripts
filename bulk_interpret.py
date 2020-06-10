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
    # print(coloc_data.keys()) ####
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
        top_snp = np.nanargmax(plasma_clust["ppas_indep"])
        top_z_beta = plasma_clust["z_beta"][top_snp]
        top_nlp_beta = np.nan_to_num(-np.log10(scipy.stats.norm.sf(abs(top_z_beta))*2))
        top_z_phi = plasma_clust["z_phi"][top_snp]
        top_nlp_phi = np.nan_to_num(-np.log10(scipy.stats.norm.sf(abs(top_z_phi))*2))
        top_z_bulk = bulk_data["z_beta"][top_snp]
        top_nlp_bulk = np.nan_to_num(-np.log10(scipy.stats.norm.sf(abs(top_z_bulk))*2))
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
            bulk_clust.get("h4_indep_eqtl"),
            bulk_clust.get("h4_ase_eqtl"),
            bulk_clust.get("h4_eqtl_eqtl"),
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
        data_lst.extend(data)

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
        "PP4Joint",
        "PP4AS",
        "PP4QTL",
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
    print(df_tr_sig) ####

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
    cluster_order = list(clusters.keys())
    storey_pis = np.zeros((len(cluster_order), 1,),)
    for ind_i, i in enumerate(cluster_order):
        df_merged = df_tr_sig.loc[df_tr_sig["Cluster"] == i]
        num_sig_train = df_merged.shape[0]
        num_sig_test = np.sum(df_merged["TopSNPNLPBulk"] >= -np.log10(0.05))
        storey_pis[ind_i, 0] = num_sig_test / num_sig_train


    title = "Bulk Replication Storey Pi"
    make_heatmap(storey_pis, cluster_order, title, os.path.join(out_dir, f"storey_pi_{sn1}.svg"))

def get_info_xval(run_name, bulk_name, genes_dir, cluster_map_path, out_dir_base):
    out_dir = os.path.join(out_dir_base, bulk_name)
    df = make_df_bulk(run_name, bulk_name, genes_dir, cluster_map_path)
    plot_xcells(df, out_dir, "Phi")
    plot_xcells(df, out_dir, "Beta")
    csv_path = os.path.join(out_dir, "cluster_info.csv")
    df.to_csv(csv_path, sep="\t", index=False, na_rep="None")
    txt_path = os.path.join(out_dir, "cluster_info.txt")
    with open(txt_path, "w") as txt_file:
        df.to_string(txt_file)

if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"

    cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")

    out_dir_base_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429"

    get_info_xval("combined", "rosmap", genes_dir_kellis, cluster_map_path_kellis, out_dir_base_kellis)



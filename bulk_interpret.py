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

def read_data(plasma_data, clusters, gene_name, top_snps=None):
    # print(coloc_data) ####
    # print(gene_name) ####
    data = []
    for c in clusters:
        plasma_clust = plasma_data.get(c, None)
        if plasma_clust is None:
            # print(plasma_data) ####
            continue
        # if top_snps:
        #     print(plasma_clust) ####
        if "causal_set_indep" in plasma_clust:
            # print(plasma_clust["ppas_indep"]) ####
            ppa = True
            if top_snps is None:
                try:
                    top_snp = np.nanargmax(plasma_clust["ppas_indep"])
                    plasma_clust["snp_ids"][top_snp]
                except (ValueError, KeyError):
                    ppa = False
            else:
                print(gene_name) ####
                try:
                    # top_snp = plasma_clust["snp_ids"].index(top_snps[c])
                    top_snp = np.argwhere(plasma_clust["snp_ids"] == top_snps[c])[0][0]
                except (ValueError, KeyError, IndexError):
                    print("False") ####
                    ppa = False
            data_clust = [
                gene_name, 
                c, 
                plasma_clust["avg_counts_total"],
                plasma_clust["avg_counts_total_scaled"],
                plasma_clust["avg_counts_mapped"],
                np.mean(plasma_clust["overdispersion"]),
                plasma_clust["avg_num_cells"],
                plasma_clust["effective_sample_size"],
                plasma_clust["sample_size"],
                plasma_clust["num_snps_informative"],
                plasma_clust["num_snps_total"],
                np.sum(plasma_clust.get("causal_set_indep", np.nan)), 
                np.sum(plasma_clust.get("causal_set_ase", np.nan)),
                np.sum(plasma_clust.get("causal_set_eqtl", np.nan)),
                np.sum(plasma_clust.get("causal_set_indep", np.nan)) / plasma_clust["num_snps_total"], 
                np.sum(plasma_clust.get("causal_set_ase", np.nan)) / plasma_clust["num_snps_total"],
                np.sum(plasma_clust.get("causal_set_eqtl", np.nan)) / plasma_clust["num_snps_total"],
                plasma_clust["ppas_indep"][top_snp] if ppa else np.nan,
                plasma_clust["z_phi"][top_snp] if ppa else np.nan,
                plasma_clust["phi"][top_snp] if ppa else np.nan,
                np.nan_to_num(-np.log10(scipy.stats.norm.sf(abs(plasma_clust["z_phi"][top_snp]))*2) if ppa else np.nan),
                # np.count_nonzero(-np.log10(scipy.stats.norm.sf(abs(plasma_clust["z_phi"]))*2/plasma_clust["num_snps_informative"]) >= 1.301) if ppa else 0,
                plasma_clust["z_beta"][top_snp] if ppa else np.nan,
                plasma_clust["beta"][top_snp] if ppa else np.nan,
                np.nan_to_num(-np.log10(scipy.stats.norm.sf(abs(plasma_clust["z_beta"][top_snp]))*2) if ppa else np.nan),
                # np.count_nonzero(-np.log10(scipy.stats.norm.sf(abs(plasma_clust["z_beta"]))*2/plasma_clust["num_snps_informative"]) >= 1.301) if ppa else 0,
                plasma_clust["snp_ids"][top_snp] if ppa else None,
                plasma_clust.get("split", np.nan),
            ]
            # print(plasma_clust["snp_ids"][top_snp] if ppa else None) ####
            # print(data_clust) ####
            data.append(data_clust)
    return data

def make_df(run_name, split, genes_dir, cluster_map_path, top_snps_dict):
    clusters = load_clusters(cluster_map_path)
    genes = os.listdir(genes_dir)
    # genes = genes[:100] ####
    data_lst = []
    for g in genes:
        if (top_snps_dict is not None) and( g not in top_snps_dict):
            continue
        gene_dir = os.path.join(genes_dir, g)
        plasma_path = os.path.join(gene_dir, run_name, "plasma_{0}.pickle")
        try:
            with open(plasma_path.format(split), "rb") as plasma_file:
                plasma_data = pickle.load(plasma_file)
        except FileNotFoundError:
            continue

        data = read_data(plasma_data, clusters, g, top_snps=(top_snps_dict[g] if top_snps_dict is not None else None))
        data_lst.extend(data)

    cols = [
        "Gene", 
        "Cluster", 
        "MeanTotalCoverage",
        "MeanTotalCoverageScaled",
        "MeanMappedCoverage",
        "MeanOverdispersion",
        "MeanCellCount",
        "UsableSampleSize",
        "TotalSampleSize",
        "UsableSNPCount",
        "TotalSNPCount",
        "CredibleSetSizeJoint", 
        "CredibleSetSizeAS",
        "CredibleSetSizeQTL",
        "CredibleSetPropJoint", 
        "CredibleSetPropAS",
        "CredibleSetPropQTL",
        "TopSNPPosterior",
        "TopSNPZPhi",
        "TopSNPPhi",
        "TopSNPNLPPhi",
        # "NumSigPhi",
        "TopSNPZBeta",
        "TopSNPBeta",
        "TopSNPNLPBeta",
        # "NumSigBeta",
        "TopSNPID",
        "Split",
    ]

    data_df = pd.DataFrame(data_lst, columns=cols)
    print(np.count_nonzero(data_df["TopSNPNLPPhi"] >= -np.log10(0.05 / data_df["UsableSNPCount"])))
    print(np.count_nonzero(data_df["TopSNPNLPBeta"] >= -np.log10(0.05 / data_df["UsableSNPCount"])))
    return data_df

def make_scatter(
        df,
        var_x,
        var_y,
        var_h,
        x_label,
        y_label, 
        h_label,
        x_lim,
        y_lim,
        title, 
        result_path,
    ):
    df_rn = df.rename(columns={var_x: x_label, var_y: y_label, var_h: h_label})
    # print(df_rn) ####
    sns.set(style="whitegrid", font="Roboto")
    f, ax = plt.subplots(figsize=(5, 5))
    # ax.set(xscale="log", yscale="log")

    sns.scatterplot(
        x=x_label, 
        y=y_label, 
        hue=h_label,
        hue_norm=(0, 10),
        data=df_rn, 
    )
    if x_lim is not None:
        plt.xlim(-x_lim, x_lim)
    if y_lim is not None:
        plt.ylim(-y_lim, y_lim)
    plt.title(title)
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()

def make_heatmap(arr, order, title, result_path):
    heat_data = pd.DataFrame(data=arr, index=order, columns=order)
    sns.set(style="whitegrid", font="Roboto")
    f, ax = plt.subplots(figsize=(5, 5))
    sns.heatmap(heat_data, annot=True, fmt=".2g", square=True, cbar=False, annot_kws={"size": 10})
    plt.title(title)
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()

def plot_xcells(df_train, df_test, out_dir, stat_name_train, stat_name_test, cutoff_train, cutoff_test):
    sn1 = stat_name_train
    sn2 = stat_name_test
    df_tr_sig = df_train.loc[
        df_train[f"TopSNPNLP{sn1}"] >= -np.log10(0.05/df_train["UsableSNPCount"])
    ]
    # df_ts_sig = df_test.loc[
    #     np.logical_and(
    #         df_train["TopSNPNLPPhi"] >= -np.log10(0.05/df_train["UsableSNPCount"]),
    #         abs(df_train["TopSNPPhi"]) <= 5
    #     )
    # ]
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
    storey_pis = np.zeros((len(cluster_order), len(cluster_order),),)
    for ind_i, i in enumerate(cluster_order):
        for ind_j, j in enumerate(cluster_order):
            df_merged = pd.merge(
                df_tr_sig.loc[df_tr_sig["Cluster"] == i], 
                df_test.loc[df_test["Cluster"] == j], 
                on=["Gene"], 
                suffixes=["_train", "_test"]
            )

            num_sig_train = df_merged.shape[0]
            num_sig_test = np.sum(df_merged[f"TopSNPNLP{sn2}_test"] >= -np.log10(0.05))
            storey_pis[ind_i, ind_j] = num_sig_test / num_sig_train

            # make_scatter(
            #     df_merged,
            #     f"TopSNP{sn1}_train",
            #     f"TopSNP{sn2}_test",
            #     f"TopSNPNLP{sn2}_test",
            #     "Train Effect Size",
            #     "Test Effect Size", 
            #     "Test -log10 P",
            #     cutoff_train,
            #     cutoff_test,
            #     "{0} to {1}".format(clusters[i], clusters[j]), 
            #     os.path.join(out_dir, "xcell_{0}_{1}.svg".format(i, j)),
            # )
            # print(i, j) ####
            # print(xw) ####
            # print(yw) ####

    title = "Cross-Cell Cross-Validation Storey Pi"
    make_heatmap(storey_pis, cluster_order, title, os.path.join(out_dir, "xcell_stats_pi.svg"))

def get_info_xval(run_name, num_splits, genes_dir, cluster_map_path, out_dir):
    df_train = make_df(run_name, "x0", genes_dir, cluster_map_path, None)
    top_snps_train = {}
    for index, row in df_train.iterrows():
        top_snps_train.setdefault(row["Gene"], {})[row["Cluster"]] = row["TopSNPID"]
    # print(top_snps_train) ####
    df_test = make_df(run_name, "x1", genes_dir, cluster_map_path, top_snps_train)
    df_comb = pd.merge(df_train, df_test, on=["Gene", "Cluster"], suffixes=["_train", "_test"])
    # print(df_train) ####
    # print(df_test) ####
    # print(df_comb) ####
    plot_xval(df_comb, os.path.join(out_dir, "xvals"))
    plot_xcells(df_train, df_test, os.path.join(out_dir, "xcells_phi"), "Phi", "Phi", 5, 5)
    plot_xcells(df_train, df_test, os.path.join(out_dir, "xcells_beta"), "Beta", "Beta", 30, 30)
    plot_xcells(df_train, df_test, os.path.join(out_dir, "xcells_phi_beta"), "Phi", "Beta", 5, 30)
    plot_xcells(df_train, df_test, os.path.join(out_dir, "xcells_beta_phi"), "Beta", "Phi", 30, 5)
    df_train.to_csv(os.path.join(out_dir, "train.csv"), sep="\t", index=False, na_rep="None")
    df_test.to_csv(os.path.join(out_dir, "test.csv"), sep="\t", index=False, na_rep="None")   

if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    # cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map.pickle")
    # genes_dir_kellis = os.path.join(data_path_kellis, "genes")

    # out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis/cluster_info.txt"

    # get_info(genes_dir_kellis, cluster_map_path_kellis, out_path_kellis)

    cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")

    out_dir_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429"

    # get_info("combined", genes_dir_kellis, cluster_map_path_kellis, out_dir_kellis)

    get_info_xval("split", 2, genes_dir_kellis, cluster_map_path_kellis, out_dir_kellis)

    # get_info_xval_nfold("split5", 5, genes_dir_kellis, cluster_map_path_kellis, out_dir_kellis)


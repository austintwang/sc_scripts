import numpy as np
import os
import pickle
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['agg.path.chunksize'] = 10000
import seaborn as sns
import pandas as pd

def read_data(plasma_data, clusters, gene_name):
    # print(coloc_data) ####
    data = []
    for c in clusters:
        plasma_clust = plasma_data.get(c, None)
        if plasma_clust is None:
            continue
        if "causal_set_indep" in plasma_clust:
            # print(plasma_clust["ppas_indep"]) ####
            ppa = True
            try:
                top_snp = np.nanargmax(plasma_clust["ppas_indep"])
            except ValueError:
                ppa = False
            data_clust = [
                gene_name, 
                c, 
                plasma_clust["avg_counts_total"],
                plasma_clust["avg_counts_total_scaled"],
                plasma_clust["avg_counts_mapped"],
                np.mean(plasma_clust["overdispersion"]),
                np.percentile(plasma_clust["overdispersion"], 25),
                np.percentile(plasma_clust["overdispersion"], 50),
                np.percentile(plasma_clust["overdispersion"], 75),
                plasma_clust["avg_num_cells"],
                plasma_clust["effective_sample_size"],
                plasma_clust["sample_size"],
                plasma_clust["num_snps_informative"],
                plasma_clust["num_snps_total"],
                np.sum(plasma_clust["causal_set_indep"]), 
                np.sum(plasma_clust["causal_set_ase"]),
                np.sum(plasma_clust["causal_set_eqtl"]),
                np.sum(plasma_clust["causal_set_indep"]) / plasma_clust["num_snps_total"], 
                np.sum(plasma_clust["causal_set_ase"]) / plasma_clust["num_snps_total"],
                np.sum(plasma_clust["causal_set_eqtl"]) / plasma_clust["num_snps_total"],
                plasma_clust["ppas_indep"][top_snp] if ppa else np.nan,
                plasma_clust["z_phi"][top_snp] if ppa else np.nan,
                plasma_clust["phi"][top_snp] if ppa else np.nan,
                plasma_clust["z_beta"][top_snp] if ppa else np.nan,
                plasma_clust["phi"][top_snp] if ppa else np.nan,
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

def make_thresh_barplot(
        df,
        var, 
        model_flavors,
        model_names, 
        threshs,
        thresh_data,
        thresh_data_models,
        title, 
        result_path,
    ):
    sns.set(style="whitegrid", font="Roboto")
    plt.figure(figsize=(4,2))
    palette = sns.cubehelix_palette(len(threshs), rot=-.25, light=.7)

    for i, t in enumerate(reversed(threshs)):
        estimator = lambda x: np.mean((x <= t).astype(int))
        # print(df[var]) ####
        # print(df[var].dtype) ####
        chart = sns.barplot(
            x=var, 
            y="Model", 
            data=df, 
            label=t, 
            order=model_flavors, 
            color=palette[i], 
            estimator=estimator,
            ci=None
        )
        plt.xlabel("Proportion of Loci")
        chart.set_yticklabels([model_names[m] for m in model_flavors])

    last_marker = [None for _ in range(len(thresh_data_models))]
    for i, t in enumerate(thresh_data[:-1]):
        for j, x in enumerate(t):
            if thresh_data_models[j] in model_flavors:
                xval = float(x)
                if (last_marker[j] is None and xval >= 0.04) or (last_marker[j] and (xval - last_marker[j]) >= 0.08):
                    plt.text(
                        xval,
                        model_flavors.index(thresh_data_models[j]),
                        threshs[i],
                        size="xx-small",
                        weight="medium",
                        ha="center",
                        va="center",
                        bbox={"boxstyle":"round", "pad":.25, "fc":"white", "ec":"white"}
                    )
                    last_marker[j] = xval

    plt.ylabel("")
    plt.title(title)
    # plt.legend()
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()

def plot_sets(df, out_dir):
    clusters = df["Cluster"].unique()
    clusters = {
        "_All": "All Cells",
        "Ex": "Excitatory Neuron",
        "Oligo": "Oligodendrocyte"
        "Astro": "Astroglia",
        "In": "Inhibitory Neuron",
        "Endo": "Endothelial",
        "OPC": "Oligodendrocyte Progenitor",
        "Per": "Per"
    }
    model_map = {
        "CredibleSetPropJoint": "PLASMA-J",
        "CredibleSetPropAS": "PLASMA-AS",
        "CredibleSetPropQTL": "QTL-Only"
    }
    var = "Credible Set Proportion"
    model_flavors = model_map.keys()
    model_names = model_map
    pal = sns.color_palette()
    model_colors = {
        "CredibleSetPropJoint": pal[0],
        "CredibleSetPropAS": pal[4],
        "CredibleSetPropQTL": pal[7],
    }
    for cluster in clusters.keys():
        df_clust = pd.melt(
            df.loc[df["Cluster"] == cluster], 
            id_vars=["Gene"], 
            value_vars=model_map.keys(),
            var_name="Model",
            val_name=var
        )
        title = clusters[cluster]
        make_violin(
            df_clust,
            var, 
            model_flavors,
            model_names, 
            model_colors,
            title, 
            os.path.join(out_dir, "sets_{0}.svg".format(cluster)),
        )
        make_thresh_barplot(
            df_clust,
            var, 
            model_flavors,
            model_names, 
            threshs,
            thresh_data,
            thresh_data_models,
            title, 
            os.path.join(out_dir, "thresh_{0}.svg".format(cluster)),
        )


def get_info(genes_dir, cluster_map_path, out_dir):
    clusters = load_clusters(cluster_map_path)
    genes = os.listdir(genes_dir)
    data_lst = []
    for g in genes:
        gene_dir = os.path.join(genes_dir, g)
        plasma_path = os.path.join(gene_dir, "plasma.pickle")
        if os.path.isdir(plasma_path):
            plasma_path = os.path.join(plasma_path, "output.pickle")
        try:
            with open(plasma_path, "rb") as plasma_file:
                plasma_data = pickle.load(plasma_file)
        except FileNotFoundError:
            continue

        data = read_data(plasma_data, clusters, g)
        data_lst.extend(data)

    cols = [
        "Gene", 
        "Cluster", 
        "MeanTotalCoverage",
        "MeanTotalCoverageScaled",
        "MeanMappedCoverage",
        "MeanOverdispersion",
        "25PercOverdispersion",
        "50PercOverdispersion",
        "75PercOverdispersion",
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
        "TopSNPZBeta",
        "TopSNPBeta",
    ]
    data_df = pd.DataFrame(data_lst, columns=cols)
    data_df.sort_values(by=["TopSNPPosterior"], ascending=False, inplace=True)
    csv_path = os.path.join(out_dir, "cluster_info.csv")
    data_df.to_csv(csv_path, sep="\t", index=False, na_rep="None")
    txt_path = os.path.join(out_dir, "cluster_info.txt")
    data_df.to_string(txt_path)
    plot_sets(df, out_dir)


if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    # cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map.pickle")
    # genes_dir_kellis = os.path.join(data_path_kellis, "genes")

    # out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis/cluster_info.txt"

    # get_info(genes_dir_kellis, cluster_map_path_kellis, out_path_kellis)

    cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")

    out_dir_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429"

    get_info(genes_dir_kellis, cluster_map_path_kellis, out_dir_kellis)

import numpy as np
import os
import pickle
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
                plasma_clust["ppas_indep"][top_snp] if ppa else np.nan,
                plasma_clust["z_phi"][top_snp] if ppa else np.nan,
                plasma_clust["z_beta"][top_snp] if ppa else np.nan,
            ]
            # print(data_clust) ####
            data.append(data_clust)
    return data

def load_clusters(cluster_map_path):
    with open(cluster_map_path, "rb") as cluster_map_file:
        cluster_map = pickle.load(cluster_map_file)
    return cluster_map.keys()

def get_info(genes_dir, cluster_map_path, out_path):
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
        "CredibleSetSize", 
        "TopSNPPosterior",
        "TopSNPZPhi",
        "TopSNPZBeta",
    ]
    data_df = pd.DataFrame(data_lst, columns=cols)
    data_df.sort_values(by=["TopSNPPosterior"], ascending=False, inplace=True)
    data_df.to_csv(out_path, sep="\t", index=False, na_rep="None")


if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    # cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map.pickle")
    # genes_dir_kellis = os.path.join(data_path_kellis, "genes")

    # out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis/cluster_info.txt"

    # get_info(genes_dir_kellis, cluster_map_path_kellis, out_path_kellis)

    cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")

    out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/cluster_info.txt"

    get_info(genes_dir_kellis, cluster_map_path_kellis, out_path_kellis)

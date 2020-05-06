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
        if "causal_set_indep" in plasma_clust and "h4_indep_eqtl" in coloc_clust:
            data_clust = [gene_name, c, np.mean(plasma_clust["causal_set_indep"]), coloc_clust["h4_indep_eqtl"]]
            # print(data_clust) ####
            data.append(data_clust)
    return data

def load_clusters(cluster_map_path):
    with open(cluster_map_path, "rb") as cluster_map_file:
        cluster_map = pickle.load(cluster_map_file)
    return cluster_map.keys()

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

    pp4_name = "{0} PP4".format(gwas_name)
    cols = ["Gene", "Cluster", "Credible Set Prop", pp4_name]
    data_df = pd.DataFrame(data_lst, columns=cols)
    data_df.sort_values(by=[pp4_name], ascending=False, inplace=True)
    data_df.to_csv(os.path.join(out_dir, gwas_name + ".csv"), index=False)
    with open(os.path.join(out_dir, gwas_name + ".txt"), "w") as txt_file:
        data_df.to_string(txt_file)


if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    # cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map.pickle")
    # genes_dir_kellis = os.path.join(data_path_kellis, "genes")

    # out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis/alz.txt"

    # interpret_genes(genes_dir_kellis, "alz", cluster_map_path_kellis, out_path_kellis)

    cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")

    out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429"

    interpret_genes(genes_dir_kellis, "alz", cluster_map_path_kellis, out_path_kellis)

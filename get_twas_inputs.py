import os
import pickle
import numpy as np

def write_gene(gene_name, gene_path_base, out_path_base):
    gene_path = os.path.join(gene_path_base, gene_name)
    plasma_path = os.path.join(gene_path, "combined", "plasma_i0.pickle") # Update after plasma rerun
    with open(plasma_path, "rb") as plasma_file:
        plasma_data = pickle.load(plasma_file)

    out_gene_dir = os.path.join(out_path_base, gene_name)
    os.makedirs(out_gene_dir, exist_ok=True)
    print(gene_name) ####

    for cluster, result in plasma_data.items():
        # print(cluster) ####
        # print(result.keys()) ####
        cluster_dir = os.path.join(out_gene_dir, cluster)
        # print(result) ####
        os.makedirs(cluster_dir, exist_ok=True)
        try:
            np.savetxt(os.path.join(cluster_dir, "hapA.txt"), result["hap_A"], fmt='%i')
            np.savetxt(os.path.join(cluster_dir, "hapB.txt"), result["hap_B"], fmt='%i')
            np.savetxt(os.path.join(cluster_dir, "countsA.txt"),  result["counts_A"])
            np.savetxt(os.path.join(cluster_dir, "countsB.txt"), result["counts_B"])
            np.savetxt(os.path.join(cluster_dir, "countsTotal.txt"), result["total_exp"])
            np.savetxt(os.path.join(cluster_dir, "sampleErr.txt"), np.sqrt(result["imbalance_errors"]))
            np.savetxt(os.path.join(cluster_dir, "snpIds.txt"), result["snp_ids"], fmt="%s")
        except KeyError:
            continue

def get_twas_inputs(gene_path_base, out_path_base):
    for gene in os.listdir(gene_path_base):
        try:
            write_gene(gene, gene_path_base, out_path_base)
        except FileNotFoundError:
          continue

if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")
    out_dir_kellis = os.path.join(data_path_kellis, "twas_429")
    get_twas_inputs(genes_dir_kellis, out_dir_kellis)
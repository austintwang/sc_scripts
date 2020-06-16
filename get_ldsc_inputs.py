#!/usr/bin/env python3

import os
import pickle
import sys
import numpy as np
import pandas as pd

def load_gene(data, cluster, gene, genes_dir):
    plasma_path = os.path.join(genes_dir, gene, "combined", "plasma_i0.pickle")
    try:
        with open(plasma_path, "rb") as plasma_file:
            plasma_data = pickle.load(plasma_file)
    except (FileNotFoundError, pickle.UnpicklingError) as e:
        # print(e) ####
        return 

    plasma_clust = plasma_data.get(cluster)
    # print(plasma_clust) ####
    if cluster is None:
        return

    informative_snps = plasma_clust["informative_snps"]
    cred = plasma_clust.get("causal_set_indep")
    if cred is None:
        return

    # print(informative_snps) ####
    for ind in informative_snps:
        causal = cred[ind]
        if causal == 0:
            continue
        contig = plasma_data["_gen"]["snp_pos"][ind][0]
        # print(plasma_data["_gen"]["snp_pos"][ind]) ####
        pos = int(plasma_data["_gen"]["snp_pos"][ind][1]) + 1
        rsid = plasma_data["_gen"]["snp_id"][ind]
        ppa = plasma_clust["ppa_indep"][ind]
        data.append([contig, pos, pos + 1, rsid, ppa, gene])

def write_bed(data, out_path):
    with open(out_path, "w") as out_file:
        out_file.writelines(f"{' '.join(map(str, i))}\n" for i in data)

def load_cluster(cluster, clusters_dir, genes_dir, out_dir):
    cluster_path = os.path.join(clusters_dir, f"{cluster}.csv")
    df = pd.read_csv(cluster_path, sep="\t")
    data = []
    # print(df.columns) ####
    for gene in df["Gene"]:
        load_gene(data, cluster, gene, genes_dir)

    out_path = os.path.join(out_dir, f"{cluster}.bed")
    write_bed(data, out_path)

def get_ldsc_inputs(clusters_dir, genes_dir, out_dir):
    for i in os.listdir(clusters_dir):
        cluster = i.split(".")[0]
        load_cluster(cluster, clusters_dir, genes_dir, out_dir)

if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    genes_dir = os.path.join(data_path_kellis, "genes_429")
    out_dir = os.path.join(data_path_kellis, "ldsc_429")
    clusters_dir = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/cell_type_spec"
    get_ldsc_inputs(clusters_dir, genes_dir, out_dir)


#!/usr/bin/env python3

import os
import pickle
import sys
import numpy as np

def load_cell_counts(cell_counts, barcodes_map):
    if isinstance(next(iter(barcodes_map.keys())), str):
        well_only = True
    else:
        well_only = False

    count_data = []
    for k, v in cell_counts.items():
        if well_only:
            bar_key = k[1]
        else:
            bar_key = k

        sample = barcodes_map[bar_key]
        counts_A, counts_B, counts_all = v
        data = (k, counts_A, counts_B, counts_all, sample)
        count_data.append(data)

    return np.stack([np.array(i, dtype="object") for i in count_data])


def write_gene(gene_name, gene_path_base, barcodes_map_path, out_path_base):
    with open(barcodes_map_path, "rb") as barcodes_map_file:
        barcodes_map = pickle.load(barcodes_map_file)
    gene_path = os.path.join(gene_path_base, gene_name)
    plasma_path = os.path.join(gene_path, "combined", "plasma_i0.pickle") # Update after plasma rerun
    with open(plasma_path, "rb") as plasma_file:
        plasma_data = pickle.load(plasma_file)
    gene_data_path = os.path.join(gene_path, "gene_data.pickle")
    with open(gene_data_path, "rb") as gene_data_file:
        gene_data = pickle.load(gene_data_file)

    out_gene_dir = os.path.join(out_path_base, gene_name)
    os.makedirs(out_gene_dir, exist_ok=True)
    # print(gene_name) ####
    # print(plasma_data.keys()) ####
    cell_counts = load_cell_counts(gene_data["cell_counts"], barcodes_map)
    np.savetxt(os.path.join(out_gene_dir, "cell_counts.txt"), cell_counts, fmt="%s")

    np.savetxt(os.path.join(out_gene_dir, "hapA.txt"), plasma_data["_gen"]["hap_A"], fmt='%i')
    np.savetxt(os.path.join(out_gene_dir, "hapB.txt"), plasma_data["_gen"]["hap_B"], fmt='%i')
    pos = plasma_data["_gen"]["snp_pos"]
    # print(pos) ####
    sid = plasma_data["_gen"]["snp_ids"]
    sal = plasma_data["_gen"]["snp_alleles"]
    snp_data = np.stack(
        [np.array([i[0][0], int(i[0][1]) + 1, i[1], i[2][0], i[2][1]], dtype='object') for i in zip(pos, sid, sal)],
    )
    np.savetxt(os.path.join(out_gene_dir, "snps.txt"), snp_data, fmt="%s")

    for cluster, result in plasma_data.items():
        cluster_dir = os.path.join(out_gene_dir, cluster)
        if cluster == "_gen":
            continue
        # print(cluster) ####
        # print(result.keys()) ####
        # print(result) ####
        try:
            counts_A = result["counts_A"]
            counts_B = result["counts_B"]
            total_exp = result["total_exp"]
            errs = np.sqrt(result["imbalance_errors"])
        except KeyError:
            continue

        os.makedirs(cluster_dir, exist_ok=True)
        # for i in os.listdir(cluster_dir):
        #     os.remove(os.path.join(cluster_dir, i))
        np.savetxt(os.path.join(cluster_dir, "countsA.txt"),  counts_A)
        np.savetxt(os.path.join(cluster_dir, "countsB.txt"), counts_B)
        np.savetxt(os.path.join(cluster_dir, "countsTotal.txt"), total_exp)
        np.savetxt(os.path.join(cluster_dir, "sampleErr.txt"), errs)

        try:
            z_phi = result["z_phi"]
            phi = result["phi"]
        except KeyError:
            continue

        np.savetxt(os.path.join(cluster_dir, "zPhi.txt"), z_phi)
        np.savetxt(os.path.join(cluster_dir, "phi.txt"), phi)

def get_twas_inputs(gene, gene_path_base, out_path_base, barcodes_map_path, status_path):
    with open(status_path, "w") as status_file:
        status_file.write("")

    try:
        write_gene(gene, gene_path_base, barcodes_map_path, out_path_base)
    except FileNotFoundError as e:
        print(e) ####
        pass

    with open(status_path, "w") as status_file:
        status_file.write("Complete")

if __name__ == '__main__':
    get_twas_inputs(*sys.argv[1:])


# def get_twas_inputs(gene_path_base, out_path_base):
#     for gene in os.listdir(gene_path_base):
#         try:
#             write_gene(gene, gene_path_base, out_path_base)
#         except FileNotFoundError:
#           continue

# if __name__ == '__main__':
#     data_path_kellis = "/agusevlab/awang/sc_kellis"
#     genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")
#     out_dir_kellis = os.path.join(data_path_kellis, "twas_429")
#     get_twas_inputs(genes_dir_kellis, out_dir_kellis)
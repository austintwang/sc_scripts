#!/usr/bin/env python3

import os
import sys
import pickle
import gzip
import glob
import numpy as np
import scipy.stats

def cluster_norm(arr):
    return arr / np.mean(arr, axis=1, keepdims=True)

def rank_norm(arr):
    return np.argsort(-arr, axis=0) / arr.shape[0]

def genes_center(arr):
    return arr - np.mean(arr, axis=0, keepdims=True)

def logtrans(arr):
    return np.log2(counts_arr + 1)

def regress_pca(arr, num_pc):
    u, s, vh = np.linalg.svd(arr)
    pcs = np.hstack([np.ones((u.shape[0], 1),), u[:,:num_pc]])
    regs, *rest = np.linalg.lstsq(pcs, arr)
    res = arr - pcs @ regs
    return res

def process(arr, flags_list):
    flag_map = {
        "c": cluster_norm,
        "r": rank_norm,
        "m": genes_center,
        "l": logtrans,
        "f": lambda x: regress_pca(x, 5),
        "t": lambda x: regress_pca(x, 10),
    }
    processed = {}
    for flags in flags_list:
        print(flags) ####
        arr_p = arr
        for f in flags:
            arr_p = flag_map[f](arr_p)
        processed[flags] = arr_p

    return processed

# def process(counts_arr):
#     # counts_norm = counts_arr / np.mean(counts_arr, axis=1, keepdims=True)
#     # logtrans = np.log2(counts_norm + 1)
#     logtrans = np.log2(counts_arr + 1)
#     # logtrans = np.argsort(-counts_arr, axis=0) / counts_arr.shape[0] ####
#     # logtrans = np.argsort(-counts_norm, axis=0) / counts_arr.shape[0] ####
#     logtrans = logtrans - np.mean(logtrans, axis=0, keepdims=True)
#     u, s, vh = np.linalg.svd(logtrans)
#     print(s[:10] / np.sum(s)) ####
#     pcs = np.hstack([np.ones((u.shape[0], 1),), u[:,:10]])
#     regs, *rest = np.linalg.lstsq(pcs, logtrans)
#     res = logtrans - pcs @ regs
#     ss_res = np.sum(res**2, axis=0) ####
#     ss_tot = np.sum((logtrans - np.mean(logtrans, axis=0, keepdims=True))**2, axis=0) ####
#     r2 = 1 - ss_res / ss_tot ####
#     print(np.nanquantile(r2, [0.25, 0.5, 0.75])) ####
#     # print(ss_res, ss_tot) ####
#     return res

def load_data(counts_paths, col_paths, row_names):
    counts_agg_dict = {}
    counts_dict = {}
    for counts_path, col_path in zip(counts_paths, col_paths):
        with gzip.open(col_path, "r") as col_file:
            col_names = col_file.read().decode('utf-8').strip().split("\n")
        # print(col_names) ####
        counts_agg_arr = np.zeros(len(col_names))
        counts_arr = np.zeros((len(col_names), len(row_names)),)

        with gzip.open(counts_path, "r") as counts_file:
            for i, cl, gl in zip(range(len(row_names)), counts_file, row_names):
                counts_gene = np.fromiter(map(float, cl.decode('utf-8').strip().split(" ")), float)
                counts_arr[:, i] = counts_gene
                counts_agg_arr += counts_gene

        for sample, counts in zip(col_names, counts_agg_arr):
            counts_agg_dict.setdefault(sample, 0)
            counts_agg_dict[sample] += counts

        for sample, counts in zip(col_names, counts_arr):
            counts_dict.setdefault(sample, 0)
            counts_dict[sample] += counts

    samples = list(counts_dict.keys())
    counts_arr = np.stack([counts_dict[i] for i in samples])
    counts_agg_arr = np.stack([counts_agg_dict[i] for i in samples])

    return samples, counts_arr, counts_agg_arr

def parse(counts_paths, col_paths, row_names, out_dir, agg_out_dir, name, flags):
    col_names_all, counts_all, counts_agg_all = load_data(counts_paths, col_paths, row_names)
    print(counts_all.shape) ####
    print(counts_all) ####

    processed = process(counts_all, flags)
    # counts_agg_out = counts_out.sum(axis=1)
    # print(counts_out) ####
    for flags, counts_out in processed.items():
        for i, gl in enumerate(row_names):
            counts_dct = dict(zip(col_names_all, counts_out[:,i]))
            # counts_dct_raw = dict(zip(col_names_all, counts_all[:,i]))
            gene = gl.strip()
            out_pattern = os.path.join(out_dir, gene + ".*")
            out_match = glob.glob(out_pattern)
            # print(out_match) ####
            if len(out_match) == 0:
                continue
            gene_counts_dir = os.path.join(out_match[0], "processed_counts")
            os.makedirs(gene_counts_dir, exist_ok=True)
            with open(os.path.join(gene_counts_dir, name), "wb") as out_file:
                pickle.dump(counts_dct, out_file)
            # with open(os.path.join(gene_counts_dir, name + '_raw'), "wb") as out_file:
            #     pickle.dump(counts_dct_raw, out_file)

    # counts_agg_dct = dict(zip(col_names_all, counts_agg_out))
    # counts_agg_dct_raw = dict(zip(col_names_all, counts_agg_all))
    # with open(os.path.join(agg_out_dir, name), "wb") as agg_out_file:
    #     pickle.dump(counts_agg_dct, agg_out_file)
    # with open(os.path.join(agg_out_dir, name + '_raw'), "wb") as agg_out_file:
    #     pickle.dump(counts_agg_dct_raw, agg_out_file)

def load_counts(base_path, rows_path, genes_dir, agg_out_dir):
    with gzip.open(rows_path, "rb") as row_file:
        row_names = row_file.read().decode('utf-8').strip().split("\n")
    counts_paths = glob.glob(os.path.join(base_path, "*.s1.gz"))
    col_paths = [i.replace(".s1.gz", ".cols.gz") for i in counts_paths]
    # print(counts_paths) ####
    # print(col_paths) ####
    parse(counts_paths, col_paths, row_names, genes_dir, agg_out_dir, "_all")

if __name__ == '__main__':
    base_dir = "/agusevlab/awang/sc_kellis"
    genes_dir = os.path.join(base_dir, "genes_429")
    agg_out_dir = os.path.join(base_dir, "agg_counts_processed")

    data_dir = os.path.join(base_dir, "snRNAseq_PFC_eQTL")
    counts_dir = os.path.join(data_dir, "result", "aggregate", "merged", "broad")
    rows_path = os.path.join(data_dir, "auxiliary", "all_features.gz")

    flags = []
    for c1 in ["", "c"]:
        for gn in ["", "r", "l"]:
            for pc in ["", "f", "t"]:
                for c2 in ["", "c"]:
                    flags.append(f"{c1}{gn}m{pc}{c2}")

    load_counts(counts_dir, rows_path, genes_dir, agg_out_dir, flags)


#!/usr/bin/env python3

import os
import sys
import pickle
import gzip
import glob
import numpy as np

def cluster_norm(arr):
    return arr / np.mean(arr, axis=1, keepdims=True)

def null_inv(arr):
    return np.ones(arr.shape) / np.mean(arr, axis=1, keepdims=True)

def rank_norm(arr):
    return np.argsort(-arr, axis=0) / arr.shape[0]

def genes_center(arr):
    return arr - np.mean(arr, axis=0, keepdims=True)

def logtrans(arr):
    return np.log2(arr + 1)

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
        "n": null_inv
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
#     logtrans = logtrans - np.mean(logtrans, axis=0, keepdims=True)
#     u, s, vh = np.linalg.svd(logtrans)
#     # print(s[:10]) ####
#     pcs = np.hstack([np.ones((u.shape[0], 1),), u[:,:10]])
#     regs, *rest = np.linalg.lstsq(pcs, logtrans)
#     res = logtrans - pcs @ regs
#     # ss_res = np.sum(res**2, axis=0) ####
#     ss_tot = np.sum((logtrans - np.mean(logtrans, axis=0, keepdims=True))**2, axis=0) ####
#     # print(1 - ss_res / ss_tot) ####
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

def parse(counts_paths, col_paths, row_names, out_dir, agg_out_dir, name, flags_list):
    # counts_agg_arrs = []
    # counts_arrs = []
    # col_names_all = []
    # for counts_path, col_path in zip(counts_paths, col_paths):
    #     with gzip.open(col_path, "r") as col_file:
    #         col_names = col_file.read().decode('utf-8').strip().split("\n")
    #     # print(col_names) ####
    #     counts_agg_arr = np.zeros(len(col_names))
    #     counts_arr = np.zeros((len(col_names), len(row_names)),)

    #     with gzip.open(counts_path, "r") as counts_file:
    #         for i, cl, gl in zip(range(len(row_names)), counts_file, row_names):
    #             counts_gene = np.fromiter(map(float, cl.decode('utf-8').strip().split(" ")), float)
    #             counts_arr[:, i] = counts_gene
    #             counts_agg_arr += counts_gene

    #     counts_agg_arrs.append(counts_agg_arr)
    #     counts_arrs.append(counts_arr)
    #     col_names_all.extend(col_names)

    # counts_agg_all = np.concatenate(counts_agg_arrs)
    # counts_all = np.concatenate(counts_arrs , axis=0)
    # # print(counts_all.shape) ####
    # # print(counts_all) ####

    col_names_all, counts_all, counts_agg_all = load_data(counts_paths, col_paths, row_names)

    processed = process(counts_all, flags_list)
    # counts_agg_out = counts_out.sum(axis=1)
    # print(counts_out) ####
    
    for i, gl in enumerate(row_names):
        out_data = {}
        for flags, counts_out in processed.items():
            counts_dct = dict(zip(col_names_all, counts_out[:,i]))
            out_data[flags] = counts_dct
        # counts_dct_raw = dict(zip(col_names_all, counts_all[:,i]))
        gene = gl.strip()
        out_pattern = os.path.join(out_dir, gene + ".*")
        out_match = glob.glob(out_pattern)
        if len(out_match) == 0:
            continue
        gene_counts_dir = os.path.join(out_match[0], "processed_counts")
        os.makedirs(gene_counts_dir, exist_ok=True)
        with open(os.path.join(gene_counts_dir, f"{name}.pickle"), "wb") as out_file:
            pickle.dump(out_data, out_file)
        # with open(os.path.join(gene_counts_dir, name + '_raw'), "wb") as out_file:
        #     pickle.dump(counts_dct, out_file)

    # counts_agg_dct = dict(zip(col_names_all, counts_agg_out))
    # counts_agg_dct_raw = dict(zip(col_names_all, counts_agg_all))
    # with open(os.path.join(agg_out_dir, name), "wb") as agg_out_file:
    #     pickle.dump(counts_agg_dct, agg_out_file)
    # with open(os.path.join(agg_out_dir, name + '_raw'), "wb") as agg_out_file:
    #     pickle.dump(counts_agg_dct_raw, agg_out_file)

def load_counts(name, patterns, base_path, rows_path, genes_dir, agg_out_dir, *args):
    with gzip.open(rows_path, "rb") as row_file:
        row_names = row_file.read().decode('utf-8').strip().split("\n")
    col_paths = []
    for p in patterns.split(","):
        counts_paths = glob.glob(os.path.join(base_path, p + ".s1.gz"))
        col_paths.extend(i.replace(".s1.gz", ".cols.gz") for i in counts_paths)
    # print(counts_paths) ####
    # print(col_paths) ####
    parse(counts_paths, col_paths, row_names, genes_dir, agg_out_dir, name, args)

if __name__ == '__main__':
    load_counts(*sys.argv[1:])


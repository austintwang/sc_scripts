#!/usr/bin/env python3

import os
import sys
import pickle
import gzip
import glob
import numpy as np

def process(counts_arr):
    logtrans = np.log2(counts_arr + 1)
    logtrans = logtrans - np.mean(logtrans, axis=0, keepdims=True)
    u, s, vh = np.linalg.svd(logtrans)
    pcs = np.hstack([np.ones((u.shape[0], 1),), u[:,:10]])
    regs, *rest = np.linalg.lstsq(pcs, logtrans)
    res = logtrans - pcs @ regs
    return res

def parse(counts_paths, col_paths, row_names, out_dir, agg_out_dir, name):
    counts_agg_arrs = []
    counts_arrs = []
    col_names_allt = []
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

        counts_agg_arrs.append(counts_agg_arr)
        counts_arrs.append(counts_arr)
        col_names_all.extend(col_names)

    counts_agg_all = np.concatenate(counts_agg_arrs)
    counts_all = np.concatenate(counts_arrs , axis=0)
    print(counts_all.shape) ####
    print(counts_all) ####

    counts_out = process(counts_all)
    counts_agg_out = counts_out.sum(axis=1)
    print(counts_out) ####
    for i, gl in enumerate(row_names):
        counts_dct = dict(zip(col_names_all, counts_out[:,i]))
        counts_dct_raw = dict(zip(col_names_all, counts_all[:,i]))
        gene = gl.strip()
        out_pattern = os.path.join(out_dir, gene + ".*")
        out_match = glob.glob(out_pattern)
        if len(out_match) == 0:
            continue
        gene_counts_dir = os.path.join(out_match[0], "processed_counts")
        os.makedirs(gene_counts_dir, exist_ok=True)
        with open(os.path.join(gene_counts_dir, name), "wb") as out_file:
            pickle.dump(counts_dct, out_file)
        with open(os.path.join(gene_counts_dir, name + '_raw'), "wb") as out_file:
            pickle.dump(counts_dct, out_file)

    counts_agg_dct = dict(zip(col_names_all, counts_agg_out))
    counts_agg_dct_raw = dict(zip(col_names_all, counts_agg_all))
    with open(os.path.join(agg_out_dir, name), "wb") as agg_out_file:
        pickle.dump(counts_agg_dct, agg_out_file)
    with open(os.path.join(agg_out_dir, name + '_raw'), "wb") as agg_out_file:
        pickle.dump(counts_agg_dct_raw, agg_out_file)

def load_counts(name, pattern, base_path, rows_path, genes_dir, agg_out_dir):
    with gzip.open(rows_path, "rb") as row_file:
        row_names = row_file.read().decode('utf-8').strip().split("\n")
    counts_paths = glob.glob(os.path.join(base_path, pattern + ".s1.gz"))
    col_paths = [i.replace(".s1.gz", ".cols.gz") for i in counts_paths]
    # print(counts_paths) ####
    # print(col_paths) ####
    parse(counts_paths, col_paths, row_names, genes_dir, agg_out_dir, name)

if __name__ == '__main__':
    load_counts(*sys.argv[1:])


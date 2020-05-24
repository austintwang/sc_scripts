#!/usr/bin/env python3

import os
import sys
import pickle
import gzip
import glob
import numpy as np

def process(counts_arr):
    logtrans = np.log2(counts_arr + 1)
    u, s, vh = np.linalg.svd(logtrans)
    pcs = np.hstack(np.array([[1]]), u[:,:10])
    regs = np.linalg.lstsq(pcs, logtrans)
    res = logtrans - pcs @ regs
    return res

def parse(counts_path, col_path, row_names, out_dir, agg_out_dir, file_name):
    with gzip.open(col_path, "r") as col_file:
        col_names = col_file.read().decode('utf-8').strip().split("\n")
    # print(col_names) ####
    counts_agg_arr = np.zeros(len(col_names))
    counts_arr = np.zeros(len(col_names), len(row_names))

    with gzip.open(counts_path, "r") as counts_file:
        for i, cl, gl in zip(range(len(row_names)), counts_file, row_names):
            counts_gene = np.fromiter(map(float, cl.decode('utf-8').strip().split(" ")), float)
            counts_arr[:, i] = counts_gene
            counts_agg_arr += counts_gene

    counts_out = process(counts_arr)
    for i, gl in enumerate(row_names):
        counts_dct = dict(zip(col_names, counts_out[:,i]))
        gene = gl.strip()
        out_pattern = os.path.join(out_dir, gene + ".*")
        out_match = glob.glob(out_pattern)
        if len(out_match) == 0:
            continue
        gene_counts_dir = os.path.join(out_match[0], "total_counts")
        os.makedirs(gene_counts_dir, exist_ok=True)
        with open(os.path.join(gene_counts_dir, file_name), "wb") as out_file:
            pickle.dump(counts_dct, out_file)

    counts_agg_dct = dict(zip(col_names, counts_agg_arr))
    with open(os.path.join(agg_out_dir, file_name), "wb") as agg_out_file:
        pickle.dump(counts_agg_dct, agg_out_file)

def load_counts(file_name, base_path, rows_path, genes_dir, agg_out_dir):
    with gzip.open(rows_path, "rb") as row_file:
        row_names = row_file.read().decode('utf-8').strip().split("\n")
    counts_path = os.path.join(base_path, file_name + ".s1.gz")
    col_path = os.path.join(base_path, file_name + ".cols.gz")
    parse(counts_path, col_path, row_names, genes_dir, agg_out_dir, file_name)

if __name__ == '__main__':
    load_counts(*sys.argv[1:])


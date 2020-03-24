import os
import sys
import pickle
import gzip
import glob
import numpy as np

def parse(counts_path, col_path, row_names, out_dir):
    with gzip.open(col_path, "r") as col_file:
        col_names = col_file.read().decode('utf-8').strip().split("\n")
    counts_agg_arr = np.zeros(len(col_names))
    with gzip.open(counts_path, "r") as counts_file:
        for cl, gl in zip(counts_file, row_names):
            counts_lst = map(float, cl.decode('utf-8').strip().split(" "))
            counts_dct = {ind: val for ind, val in enumerate(counts_lst) if val != 0}
            print(counts_lst) ####
            counts_agg_arr += np.array(counts_lst)
            gene = gl.strip()
            out_pattern = os.path.join(out_dir, gene + ".*")
            out_match = glob.glob(out_pattern)
            if len(out_match) == 0:
                continue
            with open(out_match[0]) as out_file:
                pickle.dump(counts_dct, out_file)
    counts_agg_dct = dict(zip(col_names, counts_agg_arr))
    return counts_agg_dct

def load_counts(counts_dir, rows_path, genes_dir, agg_out_path):
    with gzip.open(rows_path, "rb") as row_file:
        row_names = row_file.read().decode('utf-8').strip().split("\n")
    counts_paths = glob.glob(os.path.join(counts_dir, "*.s1.gz"))
    counts_agg = {}
    for counts_path in counts_paths:
        col_path = os.path.splitext(os.path.splitext(counts_path)[0])[0] + ".cols.gz"
        agg = parse(counts_path, col_path, row_names, genes_dir)
        counts_agg.update(agg)
    with open(agg_out_path, "wb") as agg_out_file:
        pickle.dump(counts_agg, agg_out_file)

if __name__ == '__main__':
    base_dir = "/agusevlab/awang/sc_kellis"
    genes_dir = os.path.join(base_dir, "genes")
    agg_out_path = os.path.join(base_dir, "agg_counts.pickle")

    data_dir = os.path.join(base_dir, "snRNAseq_PFC_eQTL")
    counts_dir = os.path.join(data_dir, "result", "aggregate", "merged", "broad")
    rows_path = os.path.join(data_dir, "auxiliary", "all_features.gz")

    load_counts(counts_dir, rows_path, genes_dir, agg_out_path)


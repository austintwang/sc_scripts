#!/usr/bin/env python3

import os
import pickle
import sys
import numpy as np
import pandas as pd

def load_gene(data, clusters, gene, genes_dir):
    gene_dir = os.path.join(genes_dir, gene)
    gene_path = os.path.join(gene_dir, "gene_data.pickle")

    try:
        with open(gene_path, "rb") as gene_file:
            gene_data = pickle.load(gene_file)
    except (FileNotFoundError, pickle.UnpicklingError) as e:
        # print(e) ####
        return 

    try:
        total_counts = gene_data["total_counts"]
        contig, tss_pos = gene_data["tss"]
    except Exception as e:
        print(e)
        return

    expression_dct = {}
    for c in clusters:
        try:
            # print(total_counts[c].keys()) ####
            counts = total_counts[c]["cm"]
        except (KeyError, TypeError) as e:
            print(e)
            continue
        for sample, phen in counts.items():
            expression_dct.setdefault(sample, {})[c] = phen

    expression_lst = []
    for gene, cluster_data in expression_dct.items():
        expression_lst.append([cluster_data.get(c, np.nan) for c in clusters])
    if len(expression_lst) == 0:
        return

    # print(expression_lst) ####
    expression_arr = np.array(expression_lst)
    exp_sum = np.nansum(expression_arr, axis=1, keepdims=True)
    exp_valid = ~np.isnan(expression_arr)
    num_valid = np.sum(exp_valid, axis=1, keepdims=True)
    exp_rest = (exp_sum - expression_arr) / (num_valid - 1)
    exp_diff = expression_arr - exp_rest
    diff_mean = np.nanmean(exp_diff, axis=0)
    diff_std = np.nanstd(exp_diff, axis=0)
    num_samps = np.sum(exp_valid, axis=0)
    t_scores = diff_mean / (diff_std / np.sqrt(num_samps))

    for i, c in enumerate(clusters):
        data[c].append([contig, tss_pos - 100000, tss_pos + 100000, gene, np.abs(t_scores[i])])

def write_bed(data, out_path):
    with open(out_path, "w") as out_file:
        out_file.writelines(f"{' '.join(map(str, i))}\n" for i in sorted(data))

def get_ldsc_inputs_total(clusters, genes_dir, out_dir, gene_list_path):
    data = {c: [] for c in clusters}
    with open(gene_list_path, "rb") as gene_list_file:
        gene_list = pickle.load(gene_list_file)
    gene_list = gene_list[:500] ####
    
    for gene in gene_list:
        load_gene(data, clusters, gene, genes_dir)

    for c in clusters:
        out_path = os.path.join(out_dir, f"t_{cluster}.bed")
        write_bed(data[c], out_path)


if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    genes_dir = os.path.join(data_path_kellis, "genes_429")
    out_dir = os.path.join(data_path_kellis, "ldsc_429_exp")
    gene_list_path = os.path.join(data_path_kellis, "list_429_all.pickle")
    clusters_dir = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/cell_type_spec"
    clusters = ["_all", "Ex", "Oligo", "Astro", "In", "Endo", "Microglia", "OPC", "Per"]

    get_ldsc_inputs_total(clusters, genes_dir, out_dir, gene_list_path)


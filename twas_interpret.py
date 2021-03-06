#!/usr/bin/env python3

import sys
import numpy as np
import scipy.stats
import os
import pickle
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['agg.path.chunksize'] = 10000
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def read_data(data_path, namemap):
    data_lst = []
    with open(data_path) as data_file:
        for line in data_file:
            cols = line.strip().split()
            gene = namemap.get(cols[0].split('.')[0], cols[0])
            test = cols[1]
            z = np.nan if cols[2] == "NA" else float(cols[2])
            data_lst.append([gene, test, z])

    cols = ["Gene", "Test", "Z"]
    data_df = pd.DataFrame(data_lst, columns=cols)
    return data_df

def plot_heatmap(df, result_path):
    df_plot = df.pivot(index="Gene", columns="Test", values="Z")
    # print(np.logical_not(np.isnan(df_plot).all(1))) ####
    # df_plot = df_plot[np.logical_not(np.isnan(df_plot).all(1))]
    # print(df_plot.to_numpy()) ####
    mask = np.isnan(df_plot)
    df_filled = np.abs(df_plot.fillna(df_plot.mean()))
    # print(df_filled) ####

    sns.set(style="whitegrid", font="Roboto")
    g = sns.clustermap(df_filled, mask=mask, annot=True, annot_kws={"size": 10, "weight": "medium"})
    g.savefig(result_path, bbox_inches='tight')

def twas_interpret(in_dir, in_files, namemap_path, out_dir):
    with open(namemap_path, 'rb') as namemap_file:
        namemap = pickle.load(namemap_file)
    for i in in_files:
        data_path = os.path.join(in_dir, i)
        gwas_name = i.split(".")[1]
        os.makedirs(os.path.join(out_dir, i), exist_ok=True)
        result_path = os.path.join(out_dir, i, "tops.svg")
        df = read_data(data_path, namemap)
        plot_heatmap(df, result_path)

if __name__ == '__main__':
    in_dir = "/agusevlab/awang/sc_kellis/twas_res"
    in_files = ["PLOT.AlzheimersProxyMetaIGAP_Marioni2018"]
    out_dir = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/twas"
    namemap_path = "/agusevlab/awang/ensembl/id_to_name.pickle"
    twas_interpret(in_dir, in_files, namemap_path, out_dir)

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

def plot_heatmap(df, title, result_path):
    df_plot = df.pivot(index="Study", columns="Cluster", values="Enrichment")

    sns.set(style="whitegrid", font="Roboto")
    g = sns.clustermap(df_plot, annot=True, annot_kws={"size": 10, "weight": "medium"})
    g.title(title)
    g.savefig(result_path, bbox_inches='tight')

def ldsc_interpret(in_dir, name, out_dir):
    in_path = os.path.join(in_dir, f"{name}.csv")
    df = pd.read_csv(in_path)
    for thresh, window in params:
        df_sub = df.loc[np.logical_and(df["Threshold"] == thresh, df["Window"] == window)]
        title = f"Top {thresh}, {window} kb window"
        result_path = os.path.join(out_dir, "name", f"heatmap_t_{thresh}_w_{window}.svg")
        plot_heatmap(df_sub, title, result_path)

if __name__ == '__main__':
    in_dir = "/agusevlab/awang/sc_kellis/ldsc_res/agg/"
    name = "results_chisq"
    out_dir = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/ldsc"
    ldsc_interpret(in_dir, name, out_dir)

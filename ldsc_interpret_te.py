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

STUDY_NAMES = {
    "AlzheimersMaternal_Marioni2018": "Alzheimers Maternal",
    "AlzheimersPaternal_Marioni2018": "Alzheimers Paternal",
    "AlzheimersProxyMetaIGAP_Marioni2018": "Alzheimers Proxy",
    "BD_Ruderfer2018": "Bipolar (Ruderfer)",
    "BDSCZ_Ruderfer2018": "Bipolar + Schizophrenia",
    "CD_deLange2017": "Crohn's",
    "DepressedAffect_Nagel2018": "Depressed Affect",
    "Depression_Nagel2018": "Depression (Nagel)",
    "IBD_deLange2017": "Inflammatory Bowel",
    "Intelligence_SavageJansen2018": "Intelligence",
    "MDD_Wray2018": "Depression (Wray)",
    "Neuroticism_Nagel2018": "Neuroticism",
    "PASS_Alzheimers_Jansen2019": "Alzheimers",
    "PASS_BIP_Stahl2019": "Bipolar (Stahl)",
    "PASS_Schizophrenia_Pardinas2018": "Schizophrenia (Pardinas)",
    "ReactionTime_Davies2018": "Reaction Time",
    "SCZ_Ruderfer2018": "Schizophrenia (Ruderfer)",
    "SCZvsBD_Ruderfer2018": "Bipolar vs. Schizophrenia",
    "UC_deLange2017": "Ulcerative Colitis",
    "VerbalNumericReasoning_Davies2018": "Verbal & Numeric Reasoning",
    "Worry_Nagel2018": "Worry",
}

STUDY_ORDER = [
    "PASS_Alzheimers_Jansen2019",
    "AlzheimersProxyMetaIGAP_Marioni2018",
    "AlzheimersMaternal_Marioni2018",
    "AlzheimersPaternal_Marioni2018",
    "BD_Ruderfer2018",
    "PASS_BIP_Stahl2019",
    "BDSCZ_Ruderfer2018",
    "SCZvsBD_Ruderfer2018",
    "PASS_Schizophrenia_Pardinas2018",
    "SCZ_Ruderfer2018",
    "MDD_Wray2018",
    "Depression_Nagel2018",
    "DepressedAffect_Nagel2018",
    "Neuroticism_Nagel2018",
    "Worry_Nagel2018",
    "Intelligence_SavageJansen2018",
    "VerbalNumericReasoning_Davies2018",
    "ReactionTime_Davies2018",
    "IBD_deLange2017",
    "UC_deLange2017",
    "CD_deLange2017",
]

def plot_heatmap(df, title, result_path, var, fmt=None):
    df_plot = df.pivot(index="Study", columns="Cluster", values=var).sort_index()
    # print(df_plot) ####
    df_plot = df_plot.reindex(STUDY_ORDER)
    df_plot.rename(index=STUDY_NAMES, inplace=True)
    sns.set(style="whitegrid", font="Roboto")
    # print(df_plot) ####
    kws_extras = {} if fmt is None else {"fmt": fmt}
    g = sns.heatmap(df_plot, annot=True, cmap="vlag", center=1, vmin=-1, vmax=5, annot_kws={"size": 10, "weight": "medium"}, **kws_extras)
    # g.fig.suptitle(title)
    plt.title(title)
    # g.savefig(result_path, bbox_inches='tight')
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()
    plt.close()

def ldsc_interpret(in_dir, name, params, out_dir):
    in_path = os.path.join(in_dir, f"{name}.csv")
    df = pd.read_csv(in_path)
    os.makedirs(os.path.join(out_dir, name), exist_ok=True)
    df["EnrichmentZ"] = df["Enrichment"] / df["EnrichmentStdError"]
    df["EnrichmentSig"] = df["Enrichment"].mask(df["EnrichmentP"] <= 0.05)
    # print(in_path) ####
    # print(df.to_string()) ####
    for cutoff, window in params:
        # print(thresh, window) ####
        df_sub = df.loc[np.logical_and(df["Cutoff"] == cutoff, df["Window"] == window)]
        title = f"Top {cutoff}, {window} kb window"
        result_path = os.path.join(out_dir, name, f"heatmap_c_{cutoff}_w_{window}.svg")
        plot_heatmap(df_sub, title, result_path, "Enrichment")
        result_path = os.path.join(out_dir, name, f"heatmap_p_c_{cutoff}_w_{window}.svg")
        plot_heatmap(df_sub, title, result_path, "EnrichmentP", fmt=".2e")
        result_path = os.path.join(out_dir, name, f"heatmap_z_c_{cutoff}_w_{window}.svg")
        plot_heatmap(df_sub, title, result_path, "EnrichmentZ")
        result_path = os.path.join(out_dir, name, f"heatmap_se_c_{cutoff}_w_{window}.svg")
        plot_heatmap(df_sub, title, result_path, "EnrichmentStdError")
        result_path = os.path.join(out_dir, name, f"heatmap_sg_c_{cutoff}_w_{window}.svg")
        plot_heatmap(df_sub, title, result_path, "EnrichmentSig")


if __name__ == '__main__':
    in_dir = "/agusevlab/awang/sc_kellis/ldsc_res/agg/"
    # name = "results_chisq"
    out_dir = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/ldsc"

    name = "results_total_expression"
    params = [(i, j) for i in [200, 1000] for j in [100]]
    ldsc_interpret(in_dir, name, params, out_dir)


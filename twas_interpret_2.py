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

def get_genename(eid, namemap):
    return namemap.get(eid.split('.')[0], eid)

def plot_heatmap(df, title, result_path, var):
    df_plot = df.pivot(index="Gene Name", columns="Cluster", values=var).sort_index()
    # print(df_plot) ####
    # df_plot = df_plot.reindex(STUDY_ORDER)
    # df_plot.rename(index=STUDY_NAMES, inplace=True)
    sns.set(style="whitegrid", font="Roboto")
    # print(df_plot) ####
    g = sns.heatmap(df_plot, center=0, annot=True, annot_kws={"size": 10, "weight": "medium"})
    # g.fig.suptitle(title)
    plt.title(title)
    # g.savefig(result_path, bbox_inches='tight')
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()
    plt.close()

def ldsc_interpret(in_dir, name, namemap_path, out_dir):
    with open(namemap_path, 'rb') as namemap_file:
        namemap = pickle.load(namemap_file)

    in_path = os.path.join(in_dir, f"{name}.csv")
    df = pd.read_csv(in_path)
    df["Gene Name"] = df["Gene"].map(lambda x: get_genename(x, namemap))
    # print(in_path) ####
    # print(df.to_string()) ####
    for name, group in df.groupby(["Study", "Regression"]):
        study, model = name
        study_name = STUDY_NAMES[study]
        title = f"{study_name}, {model} TWAS"

        result_path = os.path.join(out_dir, name, f"heatmap_{study}_m_{model}.svg")
        plot_heatmap(df_sub, title, result_path, "Z-Score")
        result_path = os.path.join(out_dir, name, f"heatmap_p_{study}_m_{model}.svg")
        plot_heatmap(df_sub, title, result_path, "p-Value")


if __name__ == '__main__':
    in_dir = "/agusevlab/awang/sc_kellis/ldsc_res/twas_res_2/"
    # name = "results_chisq"
    name = "SCTWAS_RESULTS"
    out_dir = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/twas/heatmaps"
    namemap_path = "/agusevlab/awang/ensembl/id_to_name.pickle"
    ldsc_interpret(in_dir, name, namemap_path, out_dir)


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

CLUSTERS =  {
    "_all": "All Cells",
    "_neur": "All Neurons",
    "_glia": "All Glia",
    "Ex": "Excitatory Neuron",
    "Oligo": "Oligodendrocyte",
    "Astro": "Astroglia",
    "In": "Inhibitory Neuron",
    "Endo": "Endothelial",
    "Microglia": "Microglia",
    "OPC": "Oligodendrocyte Progenitor",
    "Per": "Per"
}

def plot_clusters(clusters_df, rsid, gene_name, out_path):
    # print(clusters_df) ####
    title = f"{gene_name}, {rsid}"
    sns.set(style="whitegrid", font="Roboto")
    sns.violinplot(x="Cluster", y="Allelic Fraction", data=clusters_df)
    plt.title(title)
    plt.savefig(out_path, bbox_inches='tight')
    plt.clf()

def plot_fractions(props, tcounts, rsid, allele, gene_name, cluster_name, out_path):
    title = f"{gene_name}, {rsid}, {cluster_name}"
    sns.set(style="whitegrid", font="Roboto")

    x = np.arange(props.size)
    sns.scatterplot(x=x, y=props, hue=props, size=tcounts, legend=False, palette="vlag", hue_norm=(0.,1.))

    plt.axvline((props.size - 1) / 2)
    plt.axhline(0.5)
    plt.ylim(0, 1)
    plt.ylabel(f"Allelic Fraction For {allele}")
    plt.title(title)
    plt.savefig(out_path, bbox_inches='tight')
    plt.clf()

def calc_fractions(gene_id, rsid, gene_data, finemap_data, gene_map, out_dir):
    gene_name = gene_map.get(gene_id.split(".")[0], gene_id)

    snp_ids = finemap_data["_gen"]["snp_ids"]
    snp_alleles = finemap_data["_gen"]["snp_alleles"]
    hap_A = finemap_data["_gen"]["hap_A"]
    hap_B = finemap_data["_gen"]["hap_B"]
    snp_idx = np.nonzero(snp_ids == rsid)[0][0]
    # print(snp_idx) ####
    alleles = snp_alleles[snp_idx]
    phases = np.squeeze(hap_A[:,snp_idx] - hap_B[:,snp_idx])
    hets = (phases != 0)
    # print(phases) ####
    # print(hets) ####
    # print(hap_B.shape) ####

    cluster_data = []

    for cluster, fm_res in finemap_data.items():
        print(cluster)
        if cluster == "_gen":
            continue
        try:
            counts_A = fm_res["counts_A"]
            counts_B = fm_res["counts_B"]
            z_scr = fm_res["z_phi"][snp_idx]
        except KeyError as e:
            print(e)
            continue

        tcounts = counts_A + counts_B
        prop_A = (counts_A + 0.001) / (counts_A + counts_B + 0.001)
        # print(np.nansum(prop_A)) ####
        prop_alt = (prop_A * phases) % 1
        select = (hets & ~np.isnan(prop_alt) & (counts_A >= 1) & (counts_B >= 1))
        prop_hets = prop_alt[select]
        tcounts_hets = tcounts[select]
        
        if np.isnan(z_scr):
            print("z nan")
            continue
        # print(fm_res["z_phi"]) ####
        # print(snp_idx) ####
        direction = int(np.sign(z_scr))
        prop_eff = (prop_hets * direction) % 1
        ranks = np.argsort(prop_eff)[::-1]
        prop_sorted = prop_eff[ranks]
        # print(prop_sorted) ####
        tcounts_sorted = tcounts_hets[ranks]
        # print(direction) ####
        allele_eff = alleles[(1-direction)//2]

        cluster_name = CLUSTERS[cluster]
        out_path = os.path.join(out_dir, f"{gene_id}_{rsid}_{cluster}.svg")
        plot_fractions(prop_sorted, tcounts_sorted, rsid, allele_eff, gene_name, cluster_name, out_path)
        for i in prop_sorted:
            cluster_data.append([cluster, i])

    clusters_df = pd.DataFrame.from_records(cluster_data, columns=["Cluster", "Allelic Fraction"])
    out_path = os.path.join(out_dir, f"{gene_id}_{rsid}_violins.svg")
    plot_clusters(clusters_df, rsid, gene_name, out_path)


def gene_interpret(genes, data_dir, genes_map_path, run_name_plasma, out_dir):
    with open(genes_map_path, "rb") as genes_map_file:
        genes_map = pickle.load(genes_map_file)

    for gene in genes:
        gene_dir = os.path.join(data_dir, gene)
        gene_path = os.path.join(gene_dir, "gene_data.pickle")
        finemap_path = os.path.join(gene_dir, run_name_plasma, "plasma_i0.pickle")
        with open(gene_path, "rb") as gene_file:
            gene_data = pickle.load(gene_file)
        with open(finemap_path, "rb") as finemap_file:
            finemap_data = pickle.load(finemap_file)

        top_hits = np.argpartition(-finemap_data["_all"]["z_phi"]**2, 5)
        # print(top_hits) ####
        out_dir_frac = os.path.join(out_dir, "fractions")
        os.makedirs(out_dir_frac, exist_ok=True)
        for hit in top_hits[:5]:
            top_rsid = finemap_data["_gen"]["snp_ids"][hit]
            calc_fractions(gene, top_rsid, gene_data, finemap_data, genes_map, out_dir_frac)


if __name__ == '__main__':
    genes = [
        "ENSG00000120885.21_3",
        "ENSG00000087586.17_2",
        "ENSG00000100599.15_3",
        "ENSG00000167384.10_2",
        "ENSG00000168077.13_2",
        "ENSG00000104859.14_3",
        "ENSG00000168958.19_3",
        "ENSG00000047315.15_2",
        "ENSG00000090581.9_3",
        "ENSG00000171858.17_2",
        "ENSG00000100302.6_2",
        "ENSG00000005486.16_2",
        "ENSG00000144834.13_3",
    ]
    data_dir = "/agusevlab/awang/sc_kellis/genes_429"
    gene_map_path = "/agusevlab/awang/ensembl/id_to_name.pickle"
    run_name_plasma = "combined"
    out_dir = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/genes"

    gene_interpret(genes, data_dir, gene_map_path, run_name_plasma, out_dir)

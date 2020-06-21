import os
import pickle
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker 
import pandas as pd 
import scipy.stats

# def region_plotter(regions, bounds, color):
#     def region_plot(*args, **kwargs):
#         for p, q in regions:
#             if p < bounds[0]:
#                 start = bounds[0]
#             else:
#                 start = p
#             if q > bounds[1]:
#                 end = bounds[1]
#             else:
#                 end = q
#             plt.axvspan(start, end, facecolor=color, linewidth=0, alpha=0.2)

#     return region_plot

def plot_manhattan(pp_df, gene_name, out_dir):
    sns.set(style="ticks", font="Roboto")

    pal = sns.xkcd_palette(["dark slate blue", "blood red"])

    g = sns.FacetGrid(
        pp_df, 
        row="Cluster", 
        column="Source",
        # hue="Causal",
        # hue_kws={"marker":["o", "o", "D"]},
        palette=pal,
        margin_titles=True, 
        height=5, 
        aspect=1
    )

    # for k, v in regions.items():
    #     if k in annot_colormap:
    #         g.map(region_plotter(v, bounds, annot_colormap[k]))

    g.map(
        sns.scatterplot, 
        "Position", 
        "-log10 p-Value",
        # size="Causal", 
        legend=False,
        # color=".3", 
        linewidth=0,
        # hue_order=[2, 1, 0],
        # sizes={0:9, 1:12},
        s=9
    )

    x_formatter = matplotlib.ticker.ScalarFormatter()
    for i, ax in enumerate(g.fig.axes): 
        ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
        ax.xaxis.set_major_formatter(x_formatter)
    
    # plt.subplots_adjust(top=0.9, bottom = 0.13, right = 0.96)
    g.fig.suptitle(gene_name)
    plt.savefig(os.path.join(out_dir, "{0}.svg".format(gene_name)))
    plt.clf()
    plt.close()

def analyze_locus(gene_id, plasma_data, bulk_data, gene_map, out_dir):
    clusters = {
        "_all": "All Cells",
        "Ex": "Excitatory Neuron",
        "Oligo": "Oligodendrocyte",
        "Astro": "Astroglia",
        "In": "Inhibitory Neuron",
        "Endo": "Endothelial",
        "Microglia": "Microglia",
        "OPC": "Oligodendrocyte Progenitor",
        "Per": "Per"
    }
    # print(bulk_data.keys()) ####
    for clust, clust_name in clusters.items():
        plasma_clust = plasma_data.get(clust)
        if plasma_clust is None:
            continue
        # print(plasma_clust.get("run_error")) ####
        # print(plasma_clust.keys()) ####
        pp_lst = []
        try:
            for spos, z_beta, z_phi, z_bulk in zip(plasma_data["_gen"]["snp_pos"], plasma_clust["z_beta"], plasma_clust["z_phi"], bulk_data["z_beta"]):
                pos = int(spos[1]) + 1
                pp_data = [
                    [pos, clust_name, -scipy.stats.norm.logsf(np.abs(z_beta)) / np.log(10) - np.log10(2), "Single-Cell Total"],
                    [pos, clust_name, -scipy.stats.norm.logsf(np.abs(z_phi)) / np.log(10) - np.log10(2), "Single-Cell AS"],
                    [pos, clust_name, -scipy.stats.norm.logsf(np.abs(z_bulk)) / np.log(10) - np.log10(2), "Bulk"],
                ]
                pp_lst.extend(pp_data)
        except KeyError as e:
            print(e)
            continue

    pp_cols = [
        "Position", 
        "Cluster"
        "-log10 p-Value", 
        "Source"
    ]

    pp_df = pd.DataFrame(pp_lst, columns=pp_cols)

    gene_name = genes_map.get(gene_id.split(".")[0], gene_id)
    plot_manhattan(pp_df, gene_name, out_dir)

def analyze_list(gene_ids, genes_dir, gene_map_path, bulk_name, out_dir):
    with open(gene_map_path, "rb") as gene_map_file:
        gene_map = pickle.load(gene_map_file)

    data_lst = []
    for g in gene_ids:
        gene_dir = os.path.join(genes_dir, g)
        plasma_path = os.path.join(gene_dir, "combined", "plasma_i0.pickle")
        bulk_path = os.path.join(gene_dir, "bulk_qtl", f"{bulk_name}_out.pickle")
        try:
            with open(plasma_path, "rb") as plasma_file:
                plasma_data = pickle.load(plasma_file)
            with open(bulk_path, "rb") as bulk_file:
                bulk_data = pickle.load(bulk_file)
        except (FileNotFoundError, pickle.UnpicklingError):
            continue

        analyze_locus(g, plasma_data, bulk_data, gene_map, out_dir)


if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"

    genes_dir = os.path.join(data_path_kellis, "genes_429")
    gene_map_path = "/agusevlab/awang/ensembl/id_to_name.pickle"
    out_dir = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/rosmap/manhattan"

    gene_ids = [
        # "ENSG00000162601.9_2",
        "ENSG00000152642.10_3",
        "ENSG00000130749.9_2",
        "ENSG00000204120.14_3",
        "ENSG00000114388.12_2",
        "ENSG00000113070.7_2"
    ]
    bulk_name = "rosmap"

    analyze_list(gene_ids, genes_dir, gene_map_path, bulk_name, out_dir)




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

def plot_manhattan(pp_df, gene_name, gene_id, out_dir):
    # print(pp_df) ####
    sns.set(style="ticks", font="Roboto")

    pal = sns.xkcd_palette(["dark slate blue", "blood red"])

    g = sns.FacetGrid(
        pp_df, 
        row="Cluster", 
        col="Source",
        # hue="Causal",
        # hue_kws={"marker":["o", "o", "D"]},
        palette=pal,
        margin_titles=True, 
        height=1.7, 
        aspect=3
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
    plt.savefig(os.path.join(out_dir, f"{gene_id}.svg"))
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
    pp_lst = []
    for clust, clust_name in clusters.items():
        plasma_clust = plasma_data.get(clust)
        if plasma_clust is None:
            continue
        # print(plasma_clust.get("run_error")) ####
        # print(plasma_clust.keys()) ####
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
            print(clust)
            # print(plasma_clust.keys())
            continue

    pp_cols = [
        "Position", 
        "Cluster",
        "-log10 p-Value", 
        "Source"
    ]
    # print(pp_lst) ####

    pp_df = pd.DataFrame(pp_lst, columns=pp_cols)

    gene_name = gene_map.get(gene_id.split(".")[0], gene_id)
    plot_manhattan(pp_df, gene_name, gene_id, out_dir)

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
        "ENSG00000005059.15_3",
        "ENSG00000005486.16_2",
        "ENSG00000011485.14_2",
        "ENSG00000029534.19_3",
        "ENSG00000034510.5_2",
        "ENSG00000043093.13_3",
        "ENSG00000047365.11_3",
        "ENSG00000049656.13_2",
        "ENSG00000050748.17_3",
        "ENSG00000054282.15_2",
        "ENSG00000055332.16_3",
        "ENSG00000058866.14_3",
        "ENSG00000059145.18_3",
        "ENSG00000064886.13_3",
        "ENSG00000065308.4_2",
        "ENSG00000067334.13_2",
        "ENSG00000069424.14_2",
        "ENSG00000076864.19_3",
        "ENSG00000079134.11_3",
        "ENSG00000079805.16_3",
        "ENSG00000088205.12_2",
        "ENSG00000089048.14_2",
        "ENSG00000089327.14_3",
        "ENSG00000091844.7_2",
        "ENSG00000101361.16_3",
        "ENSG00000101363.12_2",
        "ENSG00000101489.19_3",
        "ENSG00000101577.9_3",
        "ENSG00000103275.19_3",
        "ENSG00000103351.12_3",
        "ENSG00000103494.12_3",
        "ENSG00000104213.12_2",
        "ENSG00000105278.10_2",
        "ENSG00000105698.15_2",
        "ENSG00000105993.14_3",
        "ENSG00000106009.15_2",
        "ENSG00000106049.8_2",
        "ENSG00000106483.11_2",
        "ENSG00000109654.14_3",
        "ENSG00000111907.20_3",
        "ENSG00000112697.15_2",
        "ENSG00000112992.16_3",
        "ENSG00000113360.16_3",
        "ENSG00000113649.11_2",
        "ENSG00000114120.11_3",
        "ENSG00000115306.15_3",
        "ENSG00000115355.15_2",
        "ENSG00000115806.12_2",
        "ENSG00000115827.13_2",
        "ENSG00000115840.13_2",
        "ENSG00000116661.9_2",
        "ENSG00000116731.22_3",
        "ENSG00000116962.14_2",
        "ENSG00000117114.19_3",
        "ENSG00000117505.12_2",
        "ENSG00000118197.13_3",
        "ENSG00000118217.5_3",
        "ENSG00000122257.18_3",
        "ENSG00000122644.12_3",
        "ENSG00000123600.19_3",
        "ENSG00000124222.22_3",
        "ENSG00000124357.12_2",
        "ENSG00000125651.13_3",
        "ENSG00000125844.15_3",
        "ENSG00000125868.15_2",
        "ENSG00000125898.12_2",
        "ENSG00000125901.5_2",
        "ENSG00000129911.8_2",
        "ENSG00000130522.5_3",
        "ENSG00000131781.12_3",
        "ENSG00000132639.12_3",
        "ENSG00000133275.15_1",
        "ENSG00000134986.13_3",
        "ENSG00000135924.15_2",
        "ENSG00000136205.16_2",
        "ENSG00000137558.7_2",
        "ENSG00000142327.12_3",
        "ENSG00000143878.9_2",
        "ENSG00000144357.16_3",
        "ENSG00000144455.13_3",
        "ENSG00000144567.10_3",
        "ENSG00000144834.13_3",
        "ENSG00000145191.12_3",
        "ENSG00000146350.13_3",
        "ENSG00000146963.17_3",
        "ENSG00000149527.17_3",
        "ENSG00000150938.9_2",
        "ENSG00000152583.12_3",
        "ENSG00000155545.19_2",
        "ENSG00000157212.18_2",
        "ENSG00000157617.16_3",
        "ENSG00000158792.15_2",
        "ENSG00000159228.12_2",
        "ENSG00000159905.14_2",
        "ENSG00000160439.15_3",
        "ENSG00000162728.4_2",
        "ENSG00000163075.12_2",
        "ENSG00000163510.13_2",
        "ENSG00000163635.17_3",
        "ENSG00000163637.12_3",
        "ENSG00000163714.17_3",
        "ENSG00000164066.12_2",
        "ENSG00000164128.6_3",
        "ENSG00000164209.16_2",
        "ENSG00000164402.13_3",
        "ENSG00000164733.20_3",
        "ENSG00000164741.14_3",
        "ENSG00000166501.12_3",
        "ENSG00000167378.8_3",
        "ENSG00000168090.9_2",
        "ENSG00000168314.17_3",
        "ENSG00000168411.13_3",
        "ENSG00000168491.9_2",
        "ENSG00000169398.19_3",
        "ENSG00000169676.5_3",
        "ENSG00000169715.14_2",
        "ENSG00000169851.15_3",
        "ENSG00000170631.14_3",
        "ENSG00000170873.18_3",
        "ENSG00000171492.14_3",
        "ENSG00000171621.13_2",
        "ENSG00000171720.9_2",
        "ENSG00000172954.13_2",
        "ENSG00000173175.14_3",
        "ENSG00000173852.14_3",
        "ENSG00000175221.14_3",
        "ENSG00000176349.11_2",
        "ENSG00000176490.4_3",
        "ENSG00000177335.10_2",
        "ENSG00000177873.12_2",
        "ENSG00000178623.11_3",
        "ENSG00000178719.16_3",
        "ENSG00000180228.12_2",
        "ENSG00000180537.12_3",
        "ENSG00000181450.17_3",
        "ENSG00000181652.19_2",
        "ENSG00000182108.9_3",
        "ENSG00000184076.13_3",
        "ENSG00000184203.7_3",
        "ENSG00000184313.19_3",
        "ENSG00000185658.13_3",
        "ENSG00000185760.15_3",
        "ENSG00000186204.14_3",
        "ENSG00000186952.14_3",
        "ENSG00000187189.10_3",
        "ENSG00000187193.8_3",
        "ENSG00000187905.10_2",
        "ENSG00000188976.10_2",
        "ENSG00000196569.11_2",
        "ENSG00000196683.10_2",
        "ENSG00000196776.14_2",
        "ENSG00000197302.10_3",
        "ENSG00000197928.10_3",
        "ENSG00000198440.9_3",
        "ENSG00000198892.6_2",
        "ENSG00000204120.14_3",
        "ENSG00000204257.14_3",
        "ENSG00000204469.12_2",
        "ENSG00000204520.12_3",
        "ENSG00000205560.12_3",
        "ENSG00000213020.9_3",
        "ENSG00000221823.10_3",
        "ENSG00000223802.7_3",
        "ENSG00000227500.9_1",
        "ENSG00000249087.6_2",
        "ENSG00000259030.6_3",
    ]
    bulk_name = "rosmap"

    analyze_list(gene_ids, genes_dir, gene_map_path, bulk_name, out_dir)




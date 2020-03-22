import os
import sys
import pickle
import numpy as np
import scipy.optimize
from scipy.special import gammaln

def nlog_likelihood(x, n, overdispersion):
    # print(overdispersion) ####
    r = (1 / overdispersion) - 1
    ll = np.sum(
        - (gammaln(n+1) + gammaln(x+r/2) + gammaln(n-x+r/2) + gammaln(r)) 
        + (gammaln(x+1) + gammaln(n-x+1) + gammaln(r/2) + gammaln(r/2) + gammaln(n+r))
    )
    return ll

def fit_overdispersion(samples):
    print(samples) ####
    x = samples[:,0]
    n = np.sum(samples, axis=1)
    res = scipy.optimize.minimize_scalar(
        lambda ovd: nlog_likelihood(x, n, ovd), 
        method="bounded", 
        bounds=(0., 1.)
    )
    if not res.success:
        print(res.message)
    return res.x

def make_cell_map(cluster_map):
    cell_map = {}
    for cluster, cells in cluster_map.items():
        for cell in cells:
            cell_map.setdefault(cell, []).append(cluster)
    return cell_map

def load_gene(clusters, gene, cell_map, barcodes_map, gene_dir, min_counts=0):
    data_path = os.path.join(gene_dir, "gene_data.pickle")
    try:
        with open(data_path, "rb") as data_file:
            data = pickle.load(data_file)
    except FileNotFoundError:
        return
    for cell, counts in data["cell_counts"].items():
        allele_counts = counts[:2]
        if np.sum(allele_counts) < min_counts:
            continue
        # ind_map = {val: ind for ind, val in enumerate(data["samples"])}
        sample = barcodes_map[cell]
        for cluster in cell_map[cell]:
            clust_data = clusters.setdefault(cluster, {})
            clust_data_ind = clust_data.setdefault(sample, {})
            gene_data = clust_data_ind.setdefault(gene, np.array([0, 0]))
            gene_data += allele_counts
            # clust_data_ind.append(counts[:2])
            # for idx, ind_counts in enumerate(counts):
            #     clust_data_ind = clust_data.setdefault(data["samples"][idx], [])
            #     clust_data_ind.append(ind_counts[:2])

def calc_overdispersions(clusters):
    overdispersions = {}
    for cluster, counts in clusters.items():
        for sample_name, genes_data in counts.items():
            samples = np.stack(genes_data.values(), axis=0)
            overdispersion = fit_overdispersion(samples)
            print(overdispersion) ####
            overdispersions.setdefault(cluster, {}).setdefault(sample_name)
            overdispersions[cluster][sample_name] = overdispersion
    return overdispersions

def get_overdispersions(data_dir, cluster_map_path, barcodes_map_path, out_path):
    with open(cluster_map_path, "rb") as cluster_map_file:
        cluster_map = pickle.load(cluster_map_file)
    with open(barcodes_map_path, "rb") as barcodes_map_file:
        barcodes_map = pickle.load(barcodes_map_file)
    cell_map = make_cell_map(cluster_map)

    genes = os.listdir(data_dir)
    clusters = {}
    for g in genes:
        load_gene(clusters, g, cell_map, barcodes_map, os.path.join(data_dir, g))

    overdispersions = calc_overdispersions(clusters)

    with open(out_path, "wb") as out_file:
        pickle.dump(overdispersions, out_file)

if __name__ == '__main__':
    base_dir = "/agusevlab/awang/sc_kellis"
    data_dir = os.path.join(base_dir, "genes")
    cluster_map_path = os.path.join(base_dir, "cluster_map.pickle")
    barcodes_map_path = os.path.join(base_dir, "metadata.pickle")
    out_path = os.path.join(base_dir, "overdispersions.pickle")
    get_overdispersions(data_dir, cluster_map_path, barcodes_map_path, out_path)
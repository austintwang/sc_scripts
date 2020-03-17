import os
import sys
import pickle
import numpy as np
import scipy.optimize
from scipy.special import gammaln

def log_likelihood(samples, overdispersion):
    r = overdispersion - 1
    x = samples[:,0]
    n = np.sum(samples, axis=1)
    ll = np.sum(
        - (gammaln(n+1) + gammaln(x+r/2) + gammaln(n-x+r/2) + gammaln(r)) 
        + (gammaln(x+1) + gammaln(n-x+1) + gammaln(r/2) + gammaln(r/2) + gammaln(n+r))
    )
    return ll

def fit_overdispersion(samples):
    res = scipy.optimize.minimize_scalar(
        lambda x: log_likelihood(samples, x), 
        method="bounded", 
        bounds=(0., 1.)
    )
    return res.x

def make_cell_map(cluster_map):
    cell_map = {}
    for cluster, cells in cluster_map.items():
        for cell in cells:
            cell_map.setdefault(cell, []).append(cluster)
    return cell_map

def load_gene(clusters, cell_map, gene_dir):
    data_path = os.path.join(gene_dir, "gene_data.pickle")
    try:
        with open(data_path, "rb") as data_file:
            data = pickle.load(data_file)
    except FileNotFoundError:
        return
    for cell, counts in data["cell_counts"].items():
        ind_map = {val: ind for ind, val in data["samples"]}
        for cluster in cell_map[cell]:
            clust_data = clusters.setdefault(cluster, {})
            for idx, ind_counts in enumerate(counts):
                clust_data_ind = clust_data.setdefault(ind_map[idx], [])
                clust_data.append(ind_counts[:2])

def calc_overdispersions(clusters):
    overdispersions = {}
    for cluster, counts in clusters.items():
        for sample_name, ind_counts in counts.items():
            samples = np.stack(ind_counts, axis=0)
            overdispersion = fit_overdispersion(samples)
            overdispersions.setdefault(cluster, {}).setdefault(sample_name)
            overdispersions[cluster][sample_name] = overdispersion
    return overdispersions

def get_overdispersions(data_dir, cluster_map_path, out_path):
    with open(cluster_map_path, "rb") as cluster_map_file:
        cluster_map = pickle.load(cluster_map_file)
    cell_map = make_cell_map(cluster_map)

    genes = os.listdir(data_dir)
    clusters = {}
    for g in genes:
        load_gene(clusters, cell_map, os.path.join(data_dir, g))

    overdispersions = calc_overdispersions(clusters)

    with open(out_path, "wb") as out_file:
        pickle.dump(overdispersions, out_file)

if __name__ == '__main__':
    base_dir = "/agusevlab/awang/sc_kellis"
    data_dir = os.path.join(base_dir, "genes")
    cluster_map_path = os.path.join(base_dir, "cluster_map.pickle")
    out_path = os.path.join(base_dir, "overdispersions.pickle")
    get_overdispersions(data_dir, cluster_map_path, out_path)
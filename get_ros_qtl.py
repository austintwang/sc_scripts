import os
import pickle
import gzip
import glob
import numpy as np
import scipy.stats

def get_ros_data(bulk_path, names_path, out_dir):
    with open(names_path, "rb") as names_file:
        gene_to_id = pickle.load(names_file)

    genes = {}
    not_found = set()
    with open(bulk_path) as bulk_file:
        colnames = next(bulk_file).strip().split()
        snpid = colnames.index("SNPid")
        feature = colnames.index("featureName")
        spearman = colnames.index("SpearmanRho")
        pval = colnames.index("pValue")
        for line in bulk_file:
            data = line.split()
            marker = data[snpid]
            gene = gene_to_id.get(data[feature])
            if gene is None:
                not_found.add(data[feature])
                continue
            zscr = scipy.stats.norm.ppf(float(data[pval]) / 2) * np.sign(float(data[spearman]))
            genes.setdefault(gene, {})[marker] = zscr

    print(not_found)

    for gene, markers in genes.items():
        matches = glob.glob(os.path.join(out_dir, gene + "*"))
        if len(matches) == 0:
            continue
        target_dir = os.path.join(out_dir, matches[0], "bulk_qtl")
        os.makedirs(target_dir, exist_ok=True)
        with open(os.path.join(target_dir, "rosmap_in.pickle"), "wb") as out_file:
            pickle.dump(markers, out_file)
        print(target_dir) ####

if __name__ == '__main__':

    bulk_path = "/agusevlab/awang/xqtl/eQTLs_all.txt"
    out_dir = "/agusevlab/awang/sc_kellis/genes_429"
    names_path = "/agusevlab/awang/ensembl/GRCh37.p13/name_to_id.pickle"

    get_ros_data(bulk_path, names_path, out_dir)
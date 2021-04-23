import os
import pickle
import gzip
import glob
import re
import numpy as np
import pandas as pd
# import scipy.stats

def add_hits(res_path, res_name, hits):
    study = res_name.split(".")[0]
    with open(res_path) as res_file:
        for line in res_file:
            cols = line.strip().split("\t")
            line_id = cols[1] 
            cluster = line_id.split(".")[-3]
            gene = cols[2]
            model = cols[15]

            hits.add((cluster, gene, model),)

def add_data(res_path, res_name, data_lst, hits):
    # print(res_path) ####
    # print(res_name) ####
    study = res_name.split(".")[0]
    with open(res_path) as res_file:
        header = next(res_file)
        cols = {val: ind for ind, val in enumerate(header.strip().split())}
        for line in res_file:
            vals = line.strip().split()
            line_id = vals[cols["FILE"]]
            cluster = line_id.split(".")[-3]
            gene = vals[cols["ID"]]
            model = vals[cols["MODEL"]]
            zstr = vals[cols["TWAS.Z"]]
            z = np.nan if zstr == "NA" else float(zstr)
            pstr = vals[cols["TWAS.P"]]
            p = np.nan if pstr == "NA" else float(pstr)

            sig = int((cluster, gene, model) in hits)
            data_lst.append([study, model, cluster, gene, z, p, sig])

def load_ldsc_out(name, res_dir_base, out_dir):
    res_dir = os.path.join(res_dir_base, name)
    hits = set()
    for i in glob.glob(os.path.join(res_dir, "*.ase.bonf_top")):
        basename = os.path.basename(i)
        add_hits(i, basename, hits)
    data_lst = []
    for i in glob.glob(os.path.join(res_dir, "*.twas")):
        basename = os.path.basename(i)
        add_data(i, basename, data_lst, hits)
    cols = [
        "Study", 
        "Regression", 
        "Cluster", 
        "Gene", 
        "Z-Score", 
        "p-Value", 
        "Hits",
    ]
    df = pd.DataFrame(data_lst, columns=cols)
    print(df) ####

    csv_path = os.path.join(out_dir, f"{name}.csv")
    df.to_csv(csv_path, index=False, na_rep="nan")

if __name__ == '__main__':
    res_dir_base = "/agusevlab/awang/sc_kellis/twas_res_2/"
    out_dir = os.path.join(res_dir_base, "agg")

    load_ldsc_out("SCTWAS_RESULTS", res_dir_base, out_dir)
    
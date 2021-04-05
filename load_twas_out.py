import os
import pickle
import gzip
import glob
import re
import numpy as np
import pandas as pd
# import scipy.stats

def add_data(res_path, res_name, data_lst):
    study = res_name.split(".")[0]
    with open(res_path) as res_file:
        for line in res_file:
            cols = line.strip().split("\t")
            line_id = cols[0] 
            cluster = line_id.split(".")[-3]
            gene = cols[3]
            model = cols[14]
            z = int(cols[17])
            p = int(cols[18])

            data_lst.append(study, model, cluster, gene, z, p)

def load_ldsc_out(name, res_dir_base, out_dir):
    res_dir = os.path.join(res_dir_base, name)
    data_lst = []
    for i in glob.glob(os.path.join(res_dir, "*.ase.bonf_top")):
        res_path = os.path.join(res_dir, i)
        add_data(res_path, i, data_lst)
    cols = [
        "Study", 
        "Regression", 
        "Cluster", 
        "Gene", 
        "Z-Score", 
        "p-Value", 
    ]
    df = pd.DataFrame(data_lst, columns=cols)

    csv_path = os.path.join(out_dir, f"{name}.csv")
    df.to_csv(csv_path, index=False, na_rep="None")

if __name__ == '__main__':
    res_dir_base = "/agusevlab/awang/sc_kellis/twas_res_2/"
    out_dir = os.path.join(res_dir_base, "agg")

    load_ldsc_out("SCTWAS_RESULTS", res_dir_base, out_dir)
    
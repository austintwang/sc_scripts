import os
import pickle
import gzip
import glob
import re
import numpy as np
import pandas as pd
# import scipy.stats

def add_data(res_path, res_name, data_lst):
    # print(res_path) ####
    # print(res_name) ####
    study = res_name.split(".")[0]
    with open(res_path) as res_file:
        for line in res_file:
            cols = line.strip().split("\t")
            line_id = cols[1] 
            cluster = line_id.split(".")[-3]
            gene = cols[2]
            model = cols[15]
            z = float(cols[18])
            p = float(cols[19])

            data_lst.append([study, model, cluster, gene, z, p, 1.])

def load_ldsc_out(name, res_dir_base, out_dir):
    res_dir = os.path.join(res_dir_base, name)
    data_lst = []
    for i in glob.glob(os.path.join(res_dir, "*.ase.bonf_top")):
        basename = os.path.basename(i)
        add_data(i, basename, data_lst)
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
    df.to_csv(csv_path, index=False, na_rep="None")

if __name__ == '__main__':
    res_dir_base = "/agusevlab/awang/sc_kellis/twas_res_2/"
    out_dir = os.path.join(res_dir_base, "agg")

    load_ldsc_out("SCTWAS_RESULTS", res_dir_base, out_dir)
    
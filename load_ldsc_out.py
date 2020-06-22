import os
import pickle
import gzip
import glob
import re
import numpy as np
import pandas as pd
# import scipy.stats

def add_data(res_path, data_lst):
    fields = ["CT", "Prop._SNPs", "Prop._h2", "Prop._h2_std_error", "Enrichment", "Enrichment_std_error", "Enrichment_p"]
    with open(res_path) as res_file:
        header = next(res_file)
        cols = {val: ind for ind, val in enumerate(header.strip().split())}
        # print(cols) ####
        for line in res_file:
            vals = line.strip().split()
            info = vals[cols["GWAS"]]
            study, paramstr = info.split(".", 1)
            threshold_match = re.search(r"joint_sc_(.*?)(_|\.results)", paramstr)   
            if threshold_match is None:
                threshold = 0.1
            else:
                threshold = float(threshold_match.group(1))
            window_match = re.search(r"pm(.*?)kb", paramstr)
            if window_match is None:
                window = 10
            else:
                window = int(window_match.group(1))
            print(paramstr, threshold, window) ####
            data_lst.append([study, threshold, window] + [vals[cols[i]] for i in fields])

def load_ldsc_out(name, res_dir_base, out_dir):
    res_dir = os.path.join(res_dir_base, name)
    data_lst = []
    for i in glob.glob(os.path.join(res_dir, "*.tab")):
        res_path = os.path.join(res_dir, i)
        add_data(res_path, data_lst)
    cols = [
        "Study", 
        "Threshold", 
        "Window", 
        "Cluster", 
        "SNPs", 
        "H2", 
        "H2StdError", 
        "Enrichment", 
        "EnrichmentStdError",
        "EnrichmentP"
    ]
    df = pd.DataFrame(data_lst, columns=cols)

    csv_path = os.path.join(out_dir, f"{name}.csv")
    df.to_csv(csv_path, index=False, na_rep="None")

if __name__ == '__main__':
    res_dir_base = "/agusevlab/awang/sc_kellis/ldsc_res/"
    out_dir = os.path.join(res_dir_base, "agg")

    load_ldsc_out("results_chisq", res_dir_base, out_dir)
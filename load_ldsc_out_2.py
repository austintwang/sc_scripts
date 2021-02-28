import os
import pickle
import gzip
import glob
# import re
import numpy as np
import pandas as pd
# import scipy.stats

def add_data(res_path, res_name, data_lst):
    fields = ["Prop._SNPs", "Prop._h2", "Prop._h2_std_error", "Enrichment", "Enrichment_std_error", "Enrichment_p"]
    field_map = {v: i for i, v in enumerate(fields)}
    field_vals = [None for i in fields]
    with open(res_path) as res_file:
        # header = next(res_file)
        # cols = {val: ind for ind, val in enumerate(header.strip().split())}
        # print(cols) ####
        study, paramstr = res_name.split(".", 1)
        paramstr = paramstr.rsplit(".", 1)[0]
        paramstr, windowstr = paramstr.rsplit("_", 1)
        # print(paramstr, windowstr) ####
        window = int(windowstr[2:-2])
        rems = paramstr.rsplit("_", 1)
        if len(rems) == 1 or rems[0] == "":
            return
        cluster = rems[0]
        threshold = float(rems[1])
        name_info = [study, threshold, window, cluster]

        for line in res_file:
            field, val = line.strip().split()
            idx = field_map.get(field)
            if idx is None:
                continue
            field_vals[idx] = float(val)

        print(name_info, field_vals) ####
        data_lst.append(name_info + field_vals)

def load_ldsc_out(name, res_dir_base, out_dir):
    res_dir = os.path.join(res_dir_base, name)
    data_lst = []
    for i in os.listdir(res_dir):
        res_path = os.path.join(res_dir, i)
        add_data(res_path, i, data_lst)
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

    load_ldsc_out("scqtl_2020_09_28_ldsc_enrichment", res_dir_base, out_dir)
    # load_ldsc_out("results_total_expression", res_dir_base, out_dir)
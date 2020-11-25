import os
import pickle
import operator
import numpy as np
import pandas as pd

OP_MAP = {
    "lt": operator.lt,
    "le": operator.le,
    "eq": operator.eq,
    "ne": operator.ne,
    "ge": operator.ge,
    "gt": operator.gt,
}

def build_sets(data, categories):
    clinical_sets = {}
    for name, criteria in categories:
        evals = [OP_MAP[op](data[trait], val) for trait, val, op in criteria]
        idx = np.all(evals, axis=0)
        ids = set(data.loc[idx, "projid"].astype(str))
        # print(ids) ####
        # print(data.loc[data[trait] == val]["projid"].astype(str)) ####
        # print(type(ids.pop())) ####
        print(name, len(ids)) ####
        clinical_sets[name] = ids
    return clinical_sets

def cast_num_no_plus(val):
    try: 
        return float(val.rstrip("+"))
    except ValueError:
        return np.nan

def get_clinical_data(data_path, categories, out_path):
    converters = {
        "age_first_ad_dx": cast_num_no_plus,
        "age_death": cast_num_no_plus,
        "age_at_visit_max": cast_num_no_plus,
    }
    data = pd.read_csv(data_path, converters=converters)
    sets = build_sets(data, categories)
    with open(out_path, "wb") as out_file:
        pickle.dump(sets, out_file)
            
if __name__ == '__main__':
    data_path = "/agusevlab/awang/sc_kellis/rosmap_clinical/ROSMAP_clinical.csv"
    categories = [
        ("Female", (("msex", 0, "eq"),)),
        ("Male", (("msex", 1, "eq"),)),
        ("AgeUnder80", (("age_death", 80, "lt"),)),
        ("Age80To90", (("age_death", 80, "ge"), ("age_death", 90, "lt"),)),
        ("AgeOver90", (("age_death", 90, "ge"),)),
        ("ReaganNeg", (("ad_reagan", 0, "eq"),)),
        ("ReaganPos", (("ad_reagan", 1, "eq"),)),
        ("CeradNCI", (("ceradsc", 5, "eq"),)),
        ("CeradMCI", (("ceradsc", 3, "ge"), ("ceradsc", 4, "le"),)),
        ("CeradAD", (("ceradsc", 1, "ge"), ("ceradsc", 2, "le"),)),
    ]
    out_path = "/agusevlab/awang/sc_kellis/rosmap_clinical/clinical_sets.pickle"
    get_clinical_data(data_path, categories, out_path)
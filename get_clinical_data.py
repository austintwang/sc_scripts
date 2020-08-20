import os
import pickle
import operator
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
        # print(data.loc[data[trait] == val]["projid"].astype(str)) ####
        # print(type(ids.pop())) ####
        clinical_sets[name] = ids
    return clinical_sets

def get_clinical_data(data_path, categories, out_path):
    data = pd.read_csv(data_path)
    sets = build_sets(data, categories)
    with open(out_path, "wb") as out_file:
        pickle.dump(sets, out_file)
            
if __name__ == '__main__':
    data_path = "/agusevlab/awang/sc_kellis/rosmap_clinical/ROSMAP_clinical.csv"
    categories = [
        ("Female", (("msex", 0, "eq"),)),
        ("Male", (("msex", 1, "eq"),)),
    ]
    out_path = "/agusevlab/awang/sc_kellis/rosmap_clinical/clinical_sets.pickle"
    get_clinical_data(data_path, categories, out_path)
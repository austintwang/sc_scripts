import os
import pickle
import pandas as pd

def build_sets(data, categories):
    clinical_sets = {}
    for name, trait, val in categories:
        ids = set(data.loc[data[trait] == val]["projid"].astype(str))
        # print(data.loc[data[trait] == val]["projid"].astype(str)) ####
        print(type(set.pop()))
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
        ("Female", "msex", 0),
        ("Male", "msex", 0)
    ]
    out_path = "/agusevlab/awang/sc_kellis/rosmap_clinical/clinical_sets.pickle"
    get_clinical_data(data_path, categories, out_path)
import os
import sys
import pickle
import numpy as np

def load_names(names_path, out_dir):
    id_to_name = {}
    name_to_id = {}
    with open(names_path, "r") as names_file:
        for line in names_file:
            cols = line.strip("\n").split("\t")
            eid = line[0]
            name = line[1]
            id_to_name[eid] = name
            name_to_id[name] = eid

    with open(os.path.join(out_dir, "id_to_name.pickle"), "wb") as out_file:
        pickle.dump(id_to_name, out_file)
    with open(os.path.join(out_dir, "name_to_id.pickle"), "wb") as out_file:
        pickle.dump(name_to_id, out_file)

if __name__ == '__main__':
    ensembl_dir = "/agusevlab/awang/ensembl"
    names_path = os.path.join(ensembl_dir, "ensembl_names.txt")
    load_names(names_path, ensembl_dir)
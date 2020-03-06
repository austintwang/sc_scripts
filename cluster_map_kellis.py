import os
import sys
import pickle
import numpy as np

def get_cluster_map(metadata_path, out_path):
    cluster_map = {}
    with open(metadata_path, "r") as metadata_file:
        next(metadata_file)
        for line in metadata_file:
            cols = line.strip("\n").split("\t")
            barcode = cols[0].split(".")[0]
            well = cols[1] 
            clusters = cols[4:]
            for c in clusters:
                cells = cluster_map.setdefault(c, [])
                cells.append((barcode, well),)

    with open(out_path, "wb") as out_file:
        pickle.dump(cluster_map, out_file)

if __name__ == '__main__':
    base_dir = "/agusevlab/awang/sc_kellis"
    metadata_path = os.path.join(base_dir, "filtered_column_metadata.txt")
    out_path = os.path.join(base_dir, "cluster_map.pickle")
    get_cluster_map(metadata_path, out_path)
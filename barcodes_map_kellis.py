import os
import sys
import pickle
import numpy as np

def load_metadata(metadata_path, out_path):
    barcodes_map = {}
    with open(metadata_path, "r") as metadata_file:
        next(metadata_file)
        for line in metadata_file:
            cols = line.strip("\n").split("\t")
            barcode = cols[0].split(".")[0]
            well = cols[1] 
            barcodes_map[(barcode, well)] = well

    return barcodes_map

if __name__ == '__main__':
    base_dir = "/agusevlab/awang/sc_kellis"
    metadata_path = os.path.join(base_dir, "filtered_column_metadata.txt")
    out_path = os.path.join(base_dir, "metadata.pickle")
    load_metadata(metadata_path, out_path)
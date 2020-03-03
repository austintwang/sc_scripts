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
    load_gene(*sys.argv[1:])
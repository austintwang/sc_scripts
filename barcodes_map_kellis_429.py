import os
import sys
import pickle
import numpy as np

def load_metadata(bam_info_path, out_path):
    well_map = {}
    with open(bam_info_path, "r") as bam_info_file:
        next(bam_info_file)
        for line in bam_info_file:
            cols = line.strip("\n").split(",")
            sample = cols[0]
            well = cols[1]
            well_map[well] = sample

    with open(out_path, "wb") as out_file:
        pickle.dump(well_map, out_file)

if __name__ == '__main__':
    base_dir = "/agusevlab/awang/sc_kellis"
    bam_info_path = os.path.join(base_dir, "Bam_paths_432_PFC_HM_Austin.csv")
    out_path = os.path.join(base_dir, "metadata_429.pickle")
    load_metadata(bam_info_path, out_path)
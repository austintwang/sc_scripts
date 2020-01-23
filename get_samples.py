import os
import pickle

def load_samples(in_path):
    samples = {}
    with open(in_path, "r") as in_file:
        header = next(in_file)
        cols = header.split(",")
        sample_idx = cols.index("ind_cov")
        well_idx = cols.index("well")
        for line in in_file:
            entries = line.split(",")
            bc = entries[sample_idx].split("-")[0]
            samples.setdefault(entries[well_idx], []).append(bc)

    return samples


def get_samples(in_path, out_path):
    samples = load_samples(in_path)
    for v in samples.values():
        v.sort()
    with open(out_path, "wb") as out_file:
        pickle.dump(samples, out_file)

if __name__ == '__main__':
    data_dir = "/agusevlab/awang/sc_le/"
    in_path = os.path.join(data_dir, "adata_obs.V6.1.txt")
    out_path = os.path.join(data_dir, "samples.pickle")
    get_samples(in_path, out_path)
import os
import pickle

def agg_counts(agg_dir, out_path):
    aggs = {}
    for name in os.listdir(agg_dir):
        path = os.path.join(agg_dir, name)
        with open(path, "rb") as agg_file:
            agg = pickle.load(agg_file)
        aggs.update(agg)

    with open(out_path, "wb") as out_file:
        pickle.dump(aggs, out_file)

if __name__ == '__main__':
    base_dir = "/agusevlab/awang/sc_kellis"
    agg_dir = os.path.join(base_dir, "agg_counts")
    out_path = os.path.join(base_dir, "agg_counts.pickle")
    agg_counts(agg_dir, out_path)


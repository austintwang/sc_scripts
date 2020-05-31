import os
import pickle

def agg_counts(agg_dir, out_path):
    aggs = {}
    for name in os.listdir(agg_dir):
        # ctype = name.split("_")[-1]
        ctype = name
        path = os.path.join(agg_dir, name)
        with open(path, "rb") as agg_file:
            agg = pickle.load(agg_file)
        aggs.setdefault(ctype, {}).update(agg)
        # for k, v in agg.items():
        #     aggs["_all"].setdefault(k, 0) 
        #     aggs["_all"][k] += v
        # print(aggs.keys()) ####

    with open(out_path, "wb") as out_file:
        pickle.dump(aggs, out_file)

if __name__ == '__main__':
    base_dir = "/agusevlab/awang/sc_kellis"
    agg_dir = os.path.join(base_dir, "agg_counts_processed")
    out_path = os.path.join(base_dir, "agg_counts.pickle")
    agg_counts(agg_dir, out_path)


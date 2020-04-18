import os
import pickle
import gzip

def get_cluster_map(clusters_path, out_path):
    cluster_map = {"_all": []}
    with gzip.open(clusters_path, "rt", encoding='utf-8') as clusters_file:
        for line in clusters_file:
            cols = line.strip("\n").split(" ")
            barcode = cols[0].split("_")[0]
            cluster = cols[1]
            cluster_map["_all"].append(barcode)
            cells = cluster_map.setdefault(cluster, [])
            cells.append(barcode)

    with open(out_path, "wb") as out_file:
        pickle.dump(cluster_map, out_file)

if __name__ == '__main__':
    base_dir = "/agusevlab/awang/sc_kellis"
    metadata_path = os.path.join(base_dir, "broad.annot.gz")
    out_path = os.path.join(base_dir, "cluster_map_429.pickle")
    get_cluster_map(metadata_path, out_path)
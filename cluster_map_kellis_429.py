import os
import pickle
import gzip

def get_cluster_map(clusters_path, well_map_path, out_path, metas):
    with open(well_map_path, "rb") as well_map_file:
        well_map = pickle.load(well_map_file)
    sample_map = {}
    for k, v in well_map.items():
        sample_map.setdefault(v, []).append(k)
    meta_map = {}
    for k, v in metas.items():
        for c in v:
            meta_map.setdefault(c, []).append(k)
    # print(meta_map) ####

    cluster_map = {}
    with gzip.open(clusters_path, "rt", encoding='utf-8') as clusters_file:
        for line in clusters_file:
            cols = line.strip("\n").split(" ")
            barcode, sample = cols[0].split("_")
            cluster = cols[1]
            clust_metas = meta_map.get(cluster, [])
            wells = sample_map[sample]
            cells = cluster_map.setdefault(cluster, [])
            for w in wells:
                cell = (barcode, w)
                # cluster_map["_all"].append(cell)
                cells.append(cell)
                for m in clust_metas:
                    cluster_map.setdefault(m, []).append(cell)

    print(cluster_map.keys()) ####
    with open(out_path, "wb") as out_file:
        pickle.dump(cluster_map, out_file)

if __name__ == '__main__':
    base_dir = "/agusevlab/awang/sc_kellis"
    metadata_path = os.path.join(base_dir, "broad.annot.gz")
    out_path = os.path.join(base_dir, "cluster_map_429.pickle")
    well_map_path = os.path.join(base_dir, "metadata_429.pickle")
    metas = {
        "_neur": ["Ex", "In"],
        "_glia": ["Oligo", "Astro", "Microglia", "OPC"],
        "_all": ["Ex", "Oligo", "Astro", "In", "Endo", "Microglia", "OPC", "Per"]
    }
    get_cluster_map(metadata_path, well_map_path, out_path, metas)
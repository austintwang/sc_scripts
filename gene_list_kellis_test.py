import os
import pickle

def build_list(data_dir, contig, out_dir):
    gene_list = []
    for gene_name in os.listdir(data_dir):
        gene_dir = os.path.join(data_dir, gene_name)
        gene_path = os.path.join(gene_dir, "gene_data.pickle")
        try:
            with open(gene_path, "rb") as gene_file:
                gene_data = pickle.load(gene_file)
                if gene_data["markers"][0][0] == contig:
                    print(gene_data["markers"][0][0]) ####
                    gene_list.append(gene_name)
        except Exception as e:
            continue

    print(len(gene_list)) ####
    with open(os.path.join(out_dir, f"list_429_test_{contig}.pickle"), "wb") as out_file:
        pickle.dump(gene_list, out_file)
            
if __name__ == '__main__':
    # Kellis 429
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")

    build_list(genes_dir_kellis, "22", data_path_kellis)
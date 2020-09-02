import os
import pickle

def build_list(data_dir, out_dir):
    gene_list = []
    for gene_name in os.listdir(data_dir):
        gene_dir = os.path.join(data_dir, gene_name)
        gene_path = os.path.join(gene_dir, "gene_data.pickle")
        try:
            with open(gene_path, "rb") as gene_file:
                gene_data = pickle.load(gene_file)
                # print(gene_name) ####
                gene_list.append(gene_name)
        except Exception as e:
            print(e)
            continue

    print(len(gene_list)) ####
    with open(os.path.join(out_dir, "list_429_all.pickle"), "wb") as out_file:
        pickle.dump(gene_list, out_file)
            
if __name__ == '__main__':
    # Kellis 429
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")

    build_list(genes_dir_kellis, data_path_kellis)
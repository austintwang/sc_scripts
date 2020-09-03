import os
import pickle

def add_genes(list_path, gene_set):
    try:
        with open(list_path) as list_file:
            for line in list_file:
                gene_set.add(line.strip())
    except FileNotFoundError as e:
        print(e)

def build_list(gwas_dir, res_dir, out_dir):
    gene_set = set()
    studies = os.listdir(gwas_dir)
    for study in studies
        gwas_path = os.path.join(gwas_dir, study)
        gwas_name = study.split(".")[0]
        list_path = os.path.join(res_dir, f"{gwas_name}_sig", "genes.txt")
        add_genes(list_path, gene_set)

    gene_list = list(gene_set)
    print(len(gene_list)) ####
    with open(os.path.join(out_dir, "list_429_sig.pickle"), "wb") as out_file:
        pickle.dump(gene_list, out_file) 
            
if __name__ == '__main__':
    # Kellis 429
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")

    build_list(genes_dir_kellis, data_path_kellis)
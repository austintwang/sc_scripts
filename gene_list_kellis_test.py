def build_list(data_dir, contig):
    gene_list = []
    for gene_name in os.listdir(data_dir):
        gene_dir = os.path.join(data_dir, gene_name)
        gene_path = os.path.join(gene_dir, "gene_data.pickle")
        try:
            with open(gene_path, "rb") as gene_file:
                gene_data = pickle.load(gene_file)
                if gene_data["markers"][0] == contig:
                    gene_list.append(gene_name)
        except Exception as e:
            continue
            
if __name__ == '__main__':
    # Kellis 429
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")

    build_list(data_dir, "22")
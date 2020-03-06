import os
import pickle

def get_tss(genes_path, out_path):
    genes = {}
    with open(genes_path, "r") as genes_file:
        for line in genes_file:
            data = line.split("\t")
            contig = data[0]
            pos = int(data[4])
            gene = data[3].split(".")[0]
            genes[gene] = (contig, pos)

    with open(out_path, "wb") as out_file:
        pickle.dump(genes, out_file)

if __name__ == '__main__':
    genes_path = "/agusevlab/DATA/ANNOTATIONS/gencode.protein_coding.transcripts.bed"
    out_path = "/agusevlab/awang/gen_data/tss.pickle"

    get_tss(genes_path, out_path)
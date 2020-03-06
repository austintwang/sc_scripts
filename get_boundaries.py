import os
import pickle

def get_boundaries(boundaries_path, out_path):
    genes = {}
    with open(boundaries_path, "r") as boundaries_file:
        for line in boundaries_file:
            if line.startswith("##"):
                continue
            data = line.split("\t")
            if data[2] == "gene":
                contig = data[0]
                start = int(data[3])
                end = int(data[4])
                gene = data[-1].split(";")[0].split(" ")[1].strip("\"")
                genes[gene] = (contig, start, end)

    with open(out_path, "wb") as out_file:
        pickle.dump(genes, out_file)

if __name__ == '__main__':
    boundaries_path = "/agusevlab/DATA/ANNOTATIONS/gencode.v26lift37.annotation.patched_contigs.gtf"
    out_path = "/agusevlab/awang/gen_data/boundaries.pickle"

    get_boundaries(boundaries_path, out_path)
import os
import pickle
import gzip
import numpy as np

def get_gwas_data(gwas_path, out_path):
    markers = {}
    sample_sizes = []
    with gzip.open(gwas_path, 'rb') as gwas_file:
        next(gwas_file)
        for line in gwas_file:
            data = line.decode('utf-8').split("\t")
            if data[3] == 'Inf':
                continue
            marker = data[0]
            zscr = float(data[3])
            markers[marker] = zscr
            try:
                sample_sizes.append(int(float(data[4].rstrip())))
            except OverflowError as e:
                print(line.decode('utf-8'))
                print(gwas_path)
                raise e

    markers["_size"] = np.mean(sample_sizes)

    with open(out_path, "wb") as out_file:
        pickle.dump(markers, out_file)

if __name__ == '__main__':
    # gwas_dir = "/agusevlab/DATA/GWAS/"
    
    # neur_dir = os.path.join(gwas_dir, "INTERNAL")
    # gwas_files = os.listdir(neur_dir)

    # for i in gwas_files:
    #     name = i.split("_")[0]
    #     path = os.path.join(neur_dir, i)
    #     out_path = "/agusevlab/awang/gwas_data/{0}.pickle".format(name)
    #     get_gwas_data(path, out_path)

    summ_dir = "/agusevlab/awang/gwas/SUMM"
    gwas_files = os.listdir(summ_dir)

    for i in gwas_files:
        name = i.split("_")[0]
        path = os.path.join(summ_dir, i)
        out_path = "/agusevlab/awang/gwas_data/{0}.pickle".format(name)
        get_gwas_data(path, out_path)


    # gwas_path_alz = os.path.join(gwas_dir, "INTERNAL", "AlzheimersProxyMetaIGAP_Marioni2018.sumstats.gz")
    # out_path = "/agusevlab/awang/gwas_data/alz.pickle"
    # get_gwas_data(gwas_path_alz, out_path)
import os
import sys
import pickle
import subprocess
import numpy as np

def dispatch(script_path, names, dataset_name, radius, data_dir, vcf_path, barcodes_map_path, boundaries_map_path, tss_map_path, memory, fails_only=False):
    jobs = []
    for name in names:
        status_path = os.path.join(data_dir, name, "load_status.txt")
        if fails_only and os.path.isfile(status_path):
            with open(status_path) as status_file:
                if status_file.read() == "Complete":
                    continue
        if not fails_only:
            with open(status_path, "w") as status_file:
                status_file.write("")

        err_name = os.path.join(data_dir, name, "load_%j.out")
        cmd = [
            "sbatch", "--mem={0}".format(memory), "-J", name, "-o", err_name, "-x", "node02,node13",
            script_path, name, dataset_name, str(radius), data_dir, vcf_path, barcodes_map_path, boundaries_map_path, tss_map_path, status_path
        ]
        print(" ".join(cmd))
        jobs.append(cmd)

    timeout = "sbatch: error: Batch job submission failed: Socket timed out on send/recv operation"
    for i in jobs:
        while True:
            try:
                submission = subprocess.run(i, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(str(submission.stdout, 'utf-8').rstrip())
                break
            except subprocess.CalledProcessError as e:
                # print(e.stdout) ####
                err = str(e.stderr, 'utf-8').rstrip()
                print(err)
                if err == timeout:
                    print("Retrying Submit")
                    continue
                else:
                    raise e

if __name__ == '__main__':
    curr_path = os.path.abspath(os.path.dirname(__file__))
    script_path = os.path.join(curr_path, "load_gene.py")

    gen_data_path = "/agusevlab/awang/gen_data"
    boundaries_map_path = os.path.join(gen_data_path, "boundaries.pickle") 
    tss_map_path = os.path.join(gen_data_path, "tss.pickle") 
    radius = 100000

    # Kellis 48
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    barcodes_map_path_kellis = os.path.join(data_path_kellis, "metadata.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes")
    names_kellis = os.listdir(genes_dir_kellis)
    vcf_path_kellis = os.path.join(data_path_kellis, "gen", "impute", "rosmap_phased.vcf.gz")

    # dispatch(
    #     script_path, 
    #     names_kellis,
    #     "Kellis",
    #     radius, 
    #     genes_dir_kellis, 
    #     vcf_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     boundaries_map_path, 
    #     tss_map_path, 
    #     2000, 
    #     fails_only=False
    # )

    dispatch(
        script_path, 
        names_kellis,
        "Kellis",
        radius, 
        genes_dir_kellis, 
        vcf_path_kellis, 
        barcodes_map_path_kellis, 
        boundaries_map_path, 
        tss_map_path, 
        2000, 
        fails_only=True
    )








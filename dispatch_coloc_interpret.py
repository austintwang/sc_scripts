import os
import sys
import pickle
import subprocess
import time
import numpy as np

def dispatch(script_path, genes_dir, genes_list_path, status_dir, gwas_names, plasma_run_name, coloc_run_name, gene_map_path, cluster_map_path, out_path, memory, fails_only=False):
    jobs = []
    for name in gwas_names:
        os.makedirs(os.path.join(status_dir, name), exist_ok=True)
        status_path = os.path.join(status_dir, name, "coloc_interpret_status.txt")
        if fails_only and os.path.isfile(status_path):
            with open(status_path) as status_file:
                if status_file.read() == "Complete":
                    continue
        if not fails_only:
            with open(status_path, "w") as status_file:
                status_file.write("")

        err_name = os.path.join(status_dir, name, "coloc_%j.out")
        cmd = [
            "sbatch", "--mem={0}".format(memory), "-J", name, "-o", err_name, "-x", "node12,node13",
            script_path, genes_dir, genes_list_path, gene_map_path, name, plasma_run_name, coloc_run_name, cluster_map_path, out_path, status_path
        ]
        print(" ".join(cmd))
        jobs.append(cmd)

    timeout = "sbatch: error: Batch job submission failed: Socket timed out on send/recv operation"
    limit = "sbatch: error: QOSMaxSubmitJobPerUserLimit"
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
                elif err.startswith(limit):
                    print("Waiting for queue to clear...")
                    time.sleep(600)
                else:
                    raise e

if __name__ == '__main__':
    curr_path = os.path.abspath(os.path.dirname(__file__))
    script_path = os.path.join(curr_path, "coloc_interpret.py")

    gene_map_path = "/agusevlab/awang/ensembl/id_to_name.pickle"

    data_path_kellis = "/agusevlab/awang/sc_kellis"

    cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")
    out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/colocalization"
    gwas_dir = "/agusevlab/awang/gwas_data"
    gwas_names = [i.split(".")[0] for i in os.listdir(gwas_dir)]
    status_dir_kellis = os.path.join(data_path_kellis, "statuses/coloc_interpret")

    # dispatch(
    #     script_path, 
    #     genes_dir_kellis, 
    #     status_dir_kellis,
    #     gwas_names,
    #     "combined",
    #     "coloc",
    #     gene_map_path,
    #     cluster_map_path_kellis, 
    #     out_path_kellis, 
    #     7000, 
    #     fails_only=False
    # )

    # dispatch(
    #     script_path, 
    #     genes_dir_kellis, 
    #     status_dir_kellis,
    #     gwas_names,
    #     "combined",
    #     "coloc",
    #     gene_map_path,
    #     cluster_map_path_kellis, 
    #     out_path_kellis, 
    #     10000, 
    #     fails_only=True
    # )

    groups = [
        "Female",
        "Male",
        "AgeUnder80",
        "Age80To90",
        "AgeOver90",
        "ReaganNeg",
        "ReaganPos",
        "CeradNCI",
        "CeradMCI",
        "CeradAD"
    ]
    # groups = groups[:4] ####
    groups = groups[4:] ####

    genes_list_path = os.path.join(data_path_kellis, "list_429_test_1.pickle")

    for group in groups:
        out_path_group = os.path.join("/agusevlab/awang/ase_finemap_results/sc_results/kellis_429/colocalization_clinical", group)
        dispatch(
            script_path, 
            genes_dir_kellis,
            genes_list_path, 
            status_dir_kellis,
            gwas_names,
            f"clinical_coloc_{group}",
            f"clinical_coloc_res_{group}",
            gene_map_path,
            cluster_map_path_kellis, 
            out_path_group, 
            7000, 
            fails_only=False
        )




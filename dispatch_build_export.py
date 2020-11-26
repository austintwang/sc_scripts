import os
import sys
import pickle
import subprocess
import time
import numpy as np

def dispatch(script_path, names, run_name, data_dir, out_dir, barcodes_map_path, memory, fails_only=False):
    jobs = []
    for name in names:
        status_path = os.path.join(data_dir, name, "twas_in_status.txt")
        if fails_only and os.path.isfile(status_path):
            with open(status_path) as status_file:
                if status_file.read() == "Complete":
                    continue
        if not fails_only:
            with open(status_path, "w") as status_file:
                status_file.write("")

        err_name = os.path.join(data_dir, name, "twas_in_%j.out")
        cmd = [
            "sbatch", "--mem={0}".format(memory), "-J", name, "-o", err_name, "-x", "node12,node13",
            script_path, name, run_name, data_dir, out_dir, barcodes_map_path, status_path
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
    script_path = os.path.join(curr_path, "build_export.py")

    data_path_kellis = "/agusevlab/awang/sc_kellis"
    cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    barcodes_map_path_kellis = os.path.join(data_path_kellis, "metadata_429.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")
    names = os.listdir(genes_dir_kellis)
    gwas_path = "/agusevlab/awang/gwas_data"

    run_name = "combined"
    run_name_coloc = "coloc"
    out_dir_kellis = os.path.join(data_path_kellis, "export_429", "_everyone")
    os.makedirs(out_dir_kellis, exist_ok=True)
    dispatch(script_path, names, run_name, run_name_coloc, gwas_path, genes_dir_kellis, out_dir_kellis, barcodes_map_path_kellis, 2000, fails_only=False)

    groups = [
        "Female",
        "Male",
        "AgeUnder80",
        "Age80To90",
        "AgeOver90",
        "ReaganNeg",
        "ReaganPos",
        "Cerad1",
        "Cerad2",
        "Cerad3",
        "Cerad4",
    ]

    for group in groups:
        run_name = f"clinical_{group}"
        run_name_coloc = f"clinical_coloc_res_{group}"
        out_dir_kellis = os.path.join(data_path_kellis, "export_429", group)
        os.makedirs(out_dir_kellis, exist_ok=True)
        dispatch(script_path, names, run_name, run_name_coloc, gwas_path, genes_dir_kellis, out_dir_kellis, barcodes_map_path_kellis, 2000, fails_only=False)




import os
import sys
import pickle
import subprocess
import time
import numpy as np

def dispatch(script_path, names, bulk_name, data_dir, params, params_path, filter_path, boundaries_map_path, memory, fails_only=False):
    with open(params_path, "wb") as params_file:
        pickle.dump(params, params_file)

    jobs = []
    for name in names:
        status_path = os.path.join(data_dir, name, f"{bulk_name}_fm_status.txt")
        if fails_only and os.path.isfile(status_path):
            with open(status_path) as status_file:
                if status_file.read() == "Complete":
                    continue
        if not fails_only:
            with open(status_path, "w") as status_file:
                status_file.write("")

        err_name = os.path.join(data_dir, name, "coloc_%j.out")
        cmd = [
            "sbatch", "--mem={0}".format(memory), "-J", name, "-o", err_name, "-x", "node12,node13",
            script_path, name, bulk_name, data_dir, params_path, filter_path, boundaries_map_path, status_path
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
    script_path = os.path.join(curr_path, "run_plasma_bulk.py")
    gen_data_path = "/agusevlab/awang/gen_data"
    boundaries_map_path = os.path.join(gen_data_path, "boundaries.pickle") 

    # Kellis 429
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    params_path_kellis = os.path.join(data_path_kellis, "coloc_params.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")
    names_kellis = os.listdir(genes_dir_kellis)

    params_kellis = {
        "total_exp_herit_prior": 0.05,
        "imbalance_herit_prior": 0.40,
        "num_ppl_total_exp": 400, 
        "cross_corr_prior": 0.9,
        "min_causal": 1,
        "num_causal": 1.,
        "search_mode": "exhaustive",
        "max_causal": 1,
        "confidence": 0.95, 
        "model_flavors_qtl": "all",
        "model_flavors_gwas": "all",
    }

    params_kellis_test = params_kellis.copy()
    params_path_kellis_test = os.path.join(data_path_kellis, "rosmap_fm_params.pickle")

    dispatch(
        script_path, 
        names_kellis, 
        "rosmap",
        genes_dir_kellis, 
        params_kellis_test, 
        params_path_kellis, 
        "all", 
        boundaries_map_path,
        2000, 
        fails_only=False
    )






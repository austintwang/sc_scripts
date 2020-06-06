import os
import sys
import pickle
import subprocess
import time
import numpy as np

def dispatch(script_path, names, data_dir, params, params_path, filter_path, gwas_dir, gwas_gen_path, boundaries_map_path, memory, fails_only=False):
    with open(params_path, "wb") as params_file:
        pickle.dump(params, params_file)

    jobs = []
    for name in names:
        status_path = os.path.join(data_dir, name, "coloc_status.txt")
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
            script_path, name, data_dir, params_path, filter_path, gwas_dir, gwas_gen_path, boundaries_map_path, status_path
        ]
        # print(" ".join(cmd))
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
    script_path = os.path.join(curr_path, "colocalize.py")

    gwas_gen_path = "/agusevlab/awang/gwas_data/gen/LDREF/1000G.EUR."
    gen_data_path = "/agusevlab/awang/gen_data"
    boundaries_map_path = os.path.join(gen_data_path, "boundaries.pickle") 

    # # Kellis 48
    # data_path_kellis = "/agusevlab/awang/sc_kellis"
    # params_path_kellis = os.path.join(data_path_kellis, "coloc_params.pickle")
    # cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map.pickle")
    # barcodes_map_path_kellis = os.path.join(data_path_kellis, "metadata.pickle")
    # overdispersion_path_kellis = os.path.join(data_path_kellis, "overdispersions.pickle")
    # genes_dir_kellis = os.path.join(data_path_kellis, "genes")
    # names_kellis = os.listdir(genes_dir_kellis)
    # vcf_path_kellis = os.path.join(data_path_kellis, "gen", "")

    # params_kellis = {
    #     "total_exp_herit_prior": 0.05,
    #     "imbalance_herit_prior": 0.40,
    #     "cross_corr_prior": 0.9,
    #     "min_causal": 1,
    #     "num_causal": 1.,
    #     "search_mode": "exhaustive",
    #     "max_causal": 1,
    #     "confidence": 0.95, 
    #     "model_flavors_qtl": "all",
    #     "model_flavors_gwas": "all",
    # }

    # # Alzheimers
    # alz_path = "/agusevlab/awang/gwas_data/alz.pickle"
    # params_kellis_alz = params_kellis.copy()
    # params_kellis_alz["num_ppl"] = 388324
    # params_path_kellis_alz = os.path.join(data_path_kellis, "plasma_c_alz_params.pickle")

    # # dispatch(
    # #     script_path, 
    # #     names_kellis, 
    # #     genes_dir_kellis, 
    # #     "alz",
    # #     params_kellis_alz, 
    # #     params_path_kellis, 
    # #     "all", 
    # #     alz_path,
    # #     gwas_gen_path,
    # #     boundaries_map_path,
    # #     2000, 
    # #     fails_only=False
    # # )

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     "alz",
    #     params_kellis_alz, 
    #     params_path_kellis, 
    #     "all", 
    #     alz_path,
    #     gwas_gen_path,
    #     boundaries_map_path,
    #     2000, 
    #     fails_only=True
    # )

    # Kellis 429
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    params_path_kellis = os.path.join(data_path_kellis, "coloc_params.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")
    names_kellis = os.listdir(genes_dir_kellis)

    params_kellis = {
        "total_exp_herit_prior": 0.05,
        "imbalance_herit_prior": 0.40,
        "cross_corr_prior": 0.9,
        "min_causal": 1,
        "num_causal": 1.,
        "search_mode": "exhaustive",
        "max_causal": 1,
        "confidence": 0.95, 
        "model_flavors_qtl": "all",
        "model_flavors_gwas": "all",
    }

    gwas_dir = "/agusevlab/awang/gwas_data"
    # for i in os.listdir(gwas_dir):
    #     gwas_path = os.path.join(gwas_dir, i)
    #     gwas_name = i.split(".")[0]
    #     params_kellis_test = params_kellis.copy()
    #     # params_kellis_test["num_ppl"] = 388324
    #     params_path_kellis_test = os.path.join(data_path_kellis, "plasma_c_{0}_params.pickle".format(gwas_name))

    #     dispatch(
    #         script_path, 
    #         names_kellis, 
    #         genes_dir_kellis, 
    #         gwas_name,
    #         params_kellis_test, 
    #         params_path_kellis, 
    #         "all", 
    #         gwas_path,
    #         gwas_gen_path,
    #         boundaries_map_path,
    #         2000, 
    #         fails_only=(gwas_name == "Intelligence")
    #     )

    params_kellis_test = params_kellis.copy()
    # params_kellis_test["num_ppl"] = 388324
    params_path_kellis_test = os.path.join(data_path_kellis, "plasma_c_params.pickle")

    dispatch(
        script_path, 
        names_kellis, 
        genes_dir_kellis, 
        params_kellis_test, 
        params_path_kellis, 
        "all", 
        gwas_dir,
        gwas_gen_path,
        boundaries_map_path,
        2000, 
        fails_only=False
    )



    # gwas_gen_path_ros = "/agusevlab/awang/sc_kellis/gen/impute/rosmap1709-chr"
    # alz_path = "/agusevlab/awang/gwas_data/alz.pickle"
    # params_kellis_alz = params_kellis.copy()
    # params_kellis_alz["num_ppl"] = 388324
    # params_path_kellis_ros = os.path.join(data_path_kellis, "plasma_c_ros_params.pickle")






import os
import sys
import pickle
import subprocess
import time
import numpy as np

def dispatch(script_path, names, data_dir, params, params_path, filter_path, cluster_map_path, barcodes_map_path, overdispersion_path, memory, fails_only=False):
    with open(params_path, "wb") as params_file:
        pickle.dump(params, params_file)

    jobs = []
    for name in names:
        status_path = os.path.join(data_dir, name, "plasma_status.txt")
        if fails_only and os.path.isfile(status_path):
            with open(status_path) as status_file:
                if status_file.read() == "Complete":
                    continue
        if not fails_only:
            with open(status_path, "w") as status_file:
                status_file.write("")

        err_name = os.path.join(data_dir, name, "plasma_%j.out")
        cmd = [
            "sbatch", "--mem={0}".format(memory), "-J", name, "-o", err_name, "-x", "node12,node13",
            script_path, name, data_dir, params_path, filter_path, cluster_map_path, barcodes_map_path, overdispersion_path, status_path
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
                    time.sleep(1800)
                else:
                    raise e

if __name__ == '__main__':
    curr_path = os.path.abspath(os.path.dirname(__file__))
    script_path = os.path.join(curr_path, "run_plasma.py")

    # # Kellis 48
    # data_path_kellis = "/agusevlab/awang/sc_kellis"
    # params_path_kellis = os.path.join(data_path_kellis, "plasma_params.pickle")
    # cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map.pickle")
    # barcodes_map_path_kellis = os.path.join(data_path_kellis, "metadata.pickle")
    # overdispersion_path_kellis = os.path.join(data_path_kellis, "overdispersions.pickle")
    # genes_dir_kellis = os.path.join(data_path_kellis, "genes")
    # names_kellis = os.listdir(genes_dir_kellis)
    # aliases_kellis = {
    #     "Oligo": "Oli",
    #     "OPC": "Opc",
    #     "Endo": "End",
    #     "Astro": "Ast",
    #     "Microglia": "Mic"
    # }

    # params_kellis = {
    #     "total_exp_herit_prior": 0.05,
    #     "imbalance_herit_prior": 0.40,
    #     "cross_corr_prior": 0.9,
    #     "min_causal": 1,
    #     "num_causal": 1.,
    #     "search_mode": "exhaustive",
    #     "max_causal": 1,
    #     "confidence": 0.95, 
    #     "model_flavors": "all",
    #     "cell_type_aliases": aliases_kellis
    # }

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     params_kellis, 
    #     params_path_kellis, 
    #     "all", 
    #     cluster_map_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     overdispersion_path_kellis, 
    #     2000, 
    #     fails_only=False
    # )

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     params_kellis, 
    #     params_path_kellis, 
    #     "all", 
    #     cluster_map_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     overdispersion_path_kellis, 
    #     2000, 
    #     fails_only=True
    # )


    # Kellis 429
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    # params_path_kellis = os.path.join(data_path_kellis, "plasma_params.pickle")
    cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    barcodes_map_path_kellis = os.path.join(data_path_kellis, "metadata_429.pickle")
    overdispersion_path_kellis = os.path.join(data_path_kellis, "overdispersions_429.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")
    names_kellis = os.listdir(genes_dir_kellis)
    vcf_path_kellis = os.path.join(data_path_kellis, "gen", "")
    aliases_kellis = {
        "Oligo": "Oli",
        "OPC": "Opc",
        "Endo": "End",
        "Astro": "Ast",
        "Microglia": "Mic"
    }

    params_kellis = {
        "run_name": "combined",
        "total_exp_herit_prior": 0.05,
        "imbalance_herit_prior": 0.40,
        "cross_corr_prior": 0.9,
        "min_causal": 1,
        "num_causal": 1.,
        "search_mode": "exhaustive",
        "max_causal": 1,
        "confidence": 0.95, 
        "model_flavors": "all",
        "cell_type_aliases": aliases_kellis,
        "splits": [1.],
        "min_carriers": 5,
        "pre_flags": "clmt",
    }
    params_path_kellis_combined = os.path.join(data_path_kellis, "plasma_params_combined.pickle")

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     params_kellis, 
    #     params_path_kellis_combined, 
    #     "all", 
    #     cluster_map_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     overdispersion_path_kellis, 
    #     2000, 
    #     fails_only=False
    # )

    params_kellis_coloc = params_kellis.copy()
    params_kellis_coloc.update({
        "run_name": "combined_coloc",
        "min_causal": 0,
    })
    params_path_kellis_coloc = os.path.join(data_path_kellis, "plasma_params_combined_coloc.pickle")

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     params_kellis_coloc, 
    #     params_path_kellis_coloc, 
    #     "all", 
    #     cluster_map_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     overdispersion_path_kellis, 
    #     2000, 
    #     fails_only=False
    # )

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     params_kellis_coloc, 
    #     params_path_kellis_coloc, 
    #     "all", 
    #     cluster_map_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     overdispersion_path_kellis, 
    #     5000, 
    #     fails_only=True
    # )

    params_kellis_xval = params_kellis.copy()
    params_kellis_xval.update({
        "run_name": "split",
        "splits": [0.5, 0.5],
    })
    params_path_kellis_xval = os.path.join(data_path_kellis, "plasma_params_split.pickle")

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     params_kellis_xval, 
    #     params_path_kellis_xval, 
    #     "all", 
    #     cluster_map_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     overdispersion_path_kellis, 
    #     2000, 
    #     fails_only=False
    # )

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     params_kellis, 
    #     params_path_kellis, 
    #     "all", 
    #     cluster_map_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     overdispersion_path_kellis, 
    #     5000, 
    #     fails_only=True
    # )

    params_kellis_xval = params_kellis.copy()
    params_kellis_xval.update({
        "run_name": "split_strict",
        "splits": [0.5, 0.5],
        "min_carriers": 10
    })
    params_path_kellis_xval = os.path.join(data_path_kellis, "plasma_params_split_strict.pickle")

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     params_kellis_xval, 
    #     params_path_kellis_xval, 
    #     "all", 
    #     cluster_map_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     overdispersion_path_kellis, 
    #     2000, 
    #     fails_only=False
    # )

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     params_kellis_xval, 
    #     params_path_kellis_xval, 
    #     "all", 
    #     cluster_map_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     overdispersion_path_kellis, 
    #     2000, 
    #     fails_only=False
    # )

    params_kellis_xval = params_kellis.copy()
    params_kellis_xval.update({
        "run_name": "split5",
        "splits": [0.2, 0.2, 0.2, 0.2, 0.2],
    })
    params_path_kellis_xval = os.path.join(data_path_kellis, "plasma_params_split.pickle")

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     params_kellis_xval, 
    #     params_path_kellis_xval, 
    #     "all", 
    #     cluster_map_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     overdispersion_path_kellis, 
    #     2000, 
    #     fails_only=False
    # )

    # dispatch(
    #     script_path, 
    #     names_kellis, 
    #     genes_dir_kellis, 
    #     params_kellis_xval, 
    #     params_path_kellis_xval, 
    #     "all", 
    #     cluster_map_path_kellis, 
    #     barcodes_map_path_kellis, 
    #     overdispersion_path_kellis,
    #     5000, 
    #     fails_only=True
    # )

    
    flags_lst = []
    for c1 in ["", "c"]:
        for gn in ["", "r", "l"]:
            for pc in ["", "f", "t"]:
                for c2 in ["", "c", "n"]:
                    flags_lst.append(f"{c1}{gn}m{pc}{c2}")

    # flags_lst = flags_lst[flags_lst.index("crmtc"):] ####

    names_test_path = os.path.join(data_path_kellis, "list_429_test_1.pickle")
    with open(names_test_path, "rb") as names_test_file:
        names_test = pickle.load(names_test_file)

    for flags in flags_lst:
        params_kellis_test = params_kellis.copy()
        params_kellis_test.update({
            "run_name": f"test_{flags}",
            "pre_flags": flags,
        })
        params_path_kellis_test = os.path.join(data_path_kellis, "test_429_params", f"plasma_params_test_{flags}.pickle")

        # dispatch(
        #     script_path, 
        #     names_test, 
        #     genes_dir_kellis, 
        #     params_kellis_test, 
        #     params_path_kellis_test, 
        #     "all", 
        #     cluster_map_path_kellis, 
        #     barcodes_map_path_kellis, 
        #     overdispersion_path_kellis, 
        #     2000, 
        #     fails_only=False
        # )

    # names_test_path = os.path.join(data_path_kellis, "list_429_test_1.pickle")
    # with open(names_test_path, "rb") as names_test_file:
    #     names_test = pickle.load(names_test_file)

    for flags in flags_lst:
        params_kellis_test = params_kellis.copy()
        params_kellis_test.update({
            "run_name": f"test_split_{flags}",
            "pre_flags": flags,
            "splits": [0.5, 0.5],
        })
        params_path_kellis_test = os.path.join(data_path_kellis, "test_429_params", f"plasma_params_test_split_{flags}.pickle")

        # dispatch(
        #     script_path, 
        #     names_test, 
        #     genes_dir_kellis, 
        #     params_kellis_test, 
        #     params_path_kellis_test, 
        #     "all", 
        #     cluster_map_path_kellis, 
        #     barcodes_map_path_kellis, 
        #     overdispersion_path_kellis, 
        #     2000, 
        #     fails_only=False
        # )

        # dispatch(
        #     script_path, 
        #     names_test, 
        #     genes_dir_kellis, 
        #     params_kellis_test, 
        #     params_path_kellis_test, 
        #     "all", 
        #     cluster_map_path_kellis, 
        #     barcodes_map_path_kellis, 
        #     overdispersion_path_kellis, 
        #     5000, 
        #     fails_only=True
        # )

    groups = [
        # "Female",
        # "Male",
        # "AgeUnder80",
        # "Age80To90",
        # "AgeOver90",
        # "ReaganNeg",
        # "ReaganPos",
        "Cerad1",
        "Cerad2",
        "Cerad3",
        "Cerad4",
    ]

    # names_test_path = os.path.join(data_path_kellis, "list_429_test_1.pickle")
    names_test_path = os.path.join(data_path_kellis, "list_429_sig.pickle")
    with open(names_test_path, "rb") as names_test_file:
        names_test = pickle.load(names_test_file)

    for group in groups:
        params_kellis_test = params_kellis.copy()
        params_kellis_test.update({
            "run_name": f"clinical_{group}",
            "clinical_group": group,
        })
        params_path_kellis_test = os.path.join(data_path_kellis, "clinical_429_params", f"plasma_params_clinical_{group}.pickle")

        dispatch(
            script_path, 
            names_test, 
            genes_dir_kellis, 
            params_kellis_test, 
            params_path_kellis_test, 
            "all", 
            cluster_map_path_kellis, 
            barcodes_map_path_kellis, 
            overdispersion_path_kellis, 
            2000, 
            fails_only=False
        )

    names_test_path = os.path.join(data_path_kellis, "list_429_sig.pickle")
    with open(names_test_path, "rb") as names_test_file:
        names_test = pickle.load(names_test_file)

    groups = [
        # "Female",
        # "Male",
        # "AgeUnder80",
        # "Age80To90",
        # "AgeOver90",
        # "ReaganNeg",
        # "ReaganPos",
        "Cerad1",
        "Cerad2",
        "Cerad3",
        "Cerad4",
    ]

    for group in groups:
        params_kellis_test = params_kellis.copy()
        params_kellis_test.update({
            "run_name": f"clinical_coloc_{group}",
            "min_causal": 0,
            "clinical_group": group,
        })
        params_path_kellis_test = os.path.join(data_path_kellis, "clinical_429_params", f"plasma_params_clinical_coloc_{group}.pickle")

        dispatch(
            script_path, 
            names_test, 
            genes_dir_kellis, 
            params_kellis_test, 
            params_path_kellis_test, 
            "all", 
            cluster_map_path_kellis, 
            barcodes_map_path_kellis, 
            overdispersion_path_kellis, 
            2000, 
            fails_only=False
        )






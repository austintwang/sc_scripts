import os
import sys
import pickle
import gzip
import glob
import subprocess
import numpy as np


def dispatch(script_path, clusters, base_path, rows_path, genes_dir, agg_out_dir, job_data_dir, flags, memory):
    jobs = []
    for file_name, patterns in clusters.items():
        os.makedirs(os.path.join(job_data_dir, file_name), exist_ok=True)
        err_name = os.path.join(job_data_dir, file_name, "load_%j.out")
        cmd = [
            "sbatch", "--mem={0}".format(memory), "-J", file_name, "-o", err_name, "-x", "node02,node13",
            script_path, file_name, ",".join(patterns), base_path, rows_path, genes_dir, agg_out_dir, *flags
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
    script_path = os.path.join(curr_path, "total_counts_kellis.py")

    base_dir = "/agusevlab/awang/sc_kellis"
    genes_dir = os.path.join(base_dir, "genes")
    agg_out_dir = os.path.join(base_dir, "agg_counts_processed")

    data_dir = os.path.join(base_dir, "snRNAseq_PFC_eQTL")
    counts_dir = os.path.join(data_dir, "result", "aggregate", "merged", "broad")
    names = {
        (os.path.splitext(os.path.splitext(os.path.basename(i))[0])[0]).split("_")[-1] 
        for i in glob.glob(os.path.join(counts_dir, "*.s1.gz"))
    }
    clusters = {
        "Ex": ["*_Ex"],
        "Oligo": ["*_Oligo"],
        "Astro": ["*_Astro"],
        "In": ["*_In"],
        "Endo": ["*_Endo"],
        "Microglia": ["*_Microglia"],
        "OPC": ["*_OPC"],
        "Per": ["*_Per"],
        "_neur": [f"*_{i}" for i in ["Ex", "In"]],
        "_glia": [f"*_{i}" for i in ["Oligo", "Astro", "Microglia", "OPC"]],
        "_all": [f"*_{i}" for i in ["Ex", "Oligo", "Astro", "In", "Endo", "Microglia", "OPC", "Per"]],
    }
    # clusters["_all"] = "*"
    rows_path = os.path.join(data_dir, "auxiliary", "all_features.gz")

    job_data_dir = os.path.join(data_dir, "job_data_processed")

    # dispatch(script_path, names, counts_dir, rows_path, genes_dir, agg_out_dir, job_data_dir, 2000)

    genes_dir = os.path.join(base_dir, "genes_429")

    flags = ["c"]
    for c1 in ["", "c"]:
        for gn in ["", "r", "l"]:
            for pc in ["", "f", "t"]:
                for c2 in ["", "c", "n"]:
                    flags.append(f"{c1}{gn}m{pc}{c2}")

    dispatch(script_path, clusters, counts_dir, rows_path, genes_dir, agg_out_dir, job_data_dir, flags, 50000)



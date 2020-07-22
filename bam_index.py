import os
import sys
import pickle
import subprocess
import numpy as np

def format_command(job_name, bam_path, out_prefix, memory):
    idx_cmd = [
        "samtools",
        "index",
        bam_path,
    ]

    err_name = out_prefix + "_%j.out"
    cmd = [
        "sbatch",
        "--mem={0}".format(memory),
        "-J",
        job_name,
        "-o",
        err_name,
        "-x", "node09,node10,node11,node19",
        "--wrap='{0}'".format(" ".join(idx_cmd)) 
    ]

    print(" ".join(cmd))

    return cmd

def dispatch_star(bam_map, out_path_base, memory, paired=False, selection=None):
    if selection is not None:
        bam_map = {k: v for k, v in bam_map.items() if k in selection}

    jobs = []
    for k, v in bam_map.items():
        out_path = os.path.join(out_path_base, k)
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        out_prefix = os.path.join(out_path, f"idx_{k}")
        cmd = format_command(k, v, out_prefix, memory)
        jobs.append(cmd)

    # print(" & ".join([" ".join(cmd) for cmd in jobs])) ####
    with open("exec.sh", "w") as script_file:
        script_file.write("#!/bin/bash\n") ####
        script_file.writelines([(" ".join(cmd) + "\n") for cmd in jobs]) ####

    subprocess.run("./exec.sh", stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # timeout = "sbatch: error: Batch job submission failed: Socket timed out on send/recv operation"
    # for i in jobs:
    #     while True:
    #         try:
    #             submission = subprocess.run(i, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #             print(str(submission.stdout, 'utf-8').rstrip())
    #             break
    #         except subprocess.CalledProcessError as e:
    #             # print(e.stdout) ####
    #             err = str(e.stderr, 'utf-8').rstrip()
    #             print(err)
    #             if err == timeout:
    #                 print("Retrying Submit")
    #                 continue
    #             else:
    #                 raise e

def get_failed_jobs(bam_map):
    fails = set()
    for k, v in bam_map.items():
        print(v + ".bai") ####
        if not os.path.isfile(v + ".bai"):
            fails.add(k)
    # for i in names:
    #     out_bam_path = os.path.join(out_path_base, i, i + "Aligned.sortedByCoord.out.bam")
    #     if not os.path.isfile(out_bam_path) or os.path.getsize(out_bam_path) < 1e5:
    #         fails.add(i)
    return fails

if __name__ == '__main__':
    # Kellis 429
    kellis_path_base = "/agusevlab/awang/sc_kellis"
    bam_path_kellis = os.path.join(kellis_path_base, "PFC_bam_files")
    bam_map_kellis_429 = {}
    with open(os.path.join(kellis_path_base, "Bam_paths_432_PFC_HM_Austin.csv")) as bam_data:
        next(bam_data)
        for line in bam_data:
            cols = line.strip().split(",")
            bam_map_kellis_429[cols[1]] = os.path.join(bam_path_kellis, cols[2].lstrip("/"), cols[3].lstrip("/"))

    out_path_base_kellis_429 = os.path.join(kellis_path_base, "processed_429")
    # print(bam_map_kellis_429) ####

    # dispatch_star(
    #     bam_map_kellis_429, out_path_base_kellis_429, 5000
    # )

    fail_kellis_429 = get_failed_jobs(bam_map_kellis_429)
    dispatch_star(
        bam_map_kellis_429, out_path_base_kellis_429, 5000, selection=fail_kellis_429
    )

    # fail_kellis_429 = get_failed_jobs(bam_map_kellis_429.keys(), out_path_base_kellis_429)
    # dispatch_star(
    #     bam_map_kellis_429, out_path_base_kellis_429, 260000, selection=fail_kellis_429
    # )



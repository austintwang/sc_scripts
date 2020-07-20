import os
import sys
import pickle
import subprocess
import numpy as np

def format_command(job_name, bam_path, contig, res_dir, err_dir, memory):
    idx_cmd = [
        "samtools",
        "view",
        "-h",
        "-b",
        "-o",
        os.path.join(res_dir, f"{job_name}_{contig}.bam"),
        bam_path,
        contig
    ]

    err_name = os.path.join(err_dir, "partition_%j.out")
    cmd = [
        "sbatch",
        "--mem={0}".format(memory),
        "-J",
        job_name,
        "-o",
        err_name,
        # "-x", "node03,node06,node07,node11,node13",
        "--wrap='{0}'".format(" ".join(idx_cmd)) 
    ]

    print(" ".join(cmd))

    return cmd

def dispatch(bam_map, out_path_base, memory, contigs, selection=None):
    if selection is not None:
        bam_map = {k: v for k, v in bam_map.items() if k in selection}

    jobs = []
    for k, v in bam_map.items():
        res_dir = os.path.join(out_path_base, k, "raw")
        err_dir = os.path.join(out_path_base, k, "status")
        os.makedirs(res_dir, exist_ok=True)
        os.makedirs(err_dir, exist_ok=True)
        for c in contigs:
            cmd = format_command(k, v, c, res_dir, err_dir, memory)
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

def get_failed_jobs(names, out_path_base):
    fails = set()
    for i in names:
        out_bam_path = os.path.join(out_path_base, i, i + "Aligned.sortedByCoord.out.bam")
        if not os.path.isfile(out_bam_path) or os.path.getsize(out_bam_path) < 1e5:
            fails.add(i)
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

    out_path_base_kellis_429 = os.path.join(kellis_path_base, "partitioned_429")
    contigs = [str(i) for i in range(1, 23)]
    # print(bam_map_kellis_429) ####

    dispatch(
        bam_map_kellis_429, out_path_base_kellis_429, contigs, 5000
    )

    # fail_kellis_429 = get_failed_jobs(bam_map_kellis_429.keys(), out_path_base_kellis_429)
    # dispatch_star(
    #     bam_map_kellis_429, out_path_base_kellis_429, 260000, selection=fail_kellis_429
    # )



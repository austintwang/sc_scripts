import os
import sys
import pickle
import subprocess
import numpy as np
import glob ####

def dispatch(script_path, dataset_name, data_dir, boundaries_path, names, contigs, out_pattern_base, memory, fails_only=False):
    jobs = []
    # print(data_dir) ####
    for name in names:
        for contig in contigs:
            split_name = f"{name}_{contig}"
            bam_path = os.path.join(data_dir, name, contig, f"{split_name}Aligned.sortedByCoord.out.bam")
            if not os.path.isfile(bam_path) or os.path.getsize(bam_path) < 1e5:
                # print(bam_path) ####
                continue

            status_path = os.path.join(data_dir, name, contig, "countstatus.txt")
            if fails_only and os.path.isfile(status_path):
                with open(status_path) as status_file:
                    # print(repr(status_file.read())) ####
                    # continue ####
                    if status_file.read() == "Complete":
                        # print("complete") ####
                        # outs = glob.glob(os.path.join(data_dir, name, "count_*.out")) ####
                        # with open(max(outs)) as of: ####
                        #     ol = of.readlines()
                        # print(ol) ####
                        # print(len(ol)) ####
                        # if len(ol) == 1:
                        continue
            if not fails_only:
                with open(status_path, "w") as status_file:
                    status_file.write("")

            err_name = os.path.join(data_dir, name, contig, "count_%j.out")
            out_pattern = out_pattern_base.format(split_name)
            cmd = [
                "sbatch", "--mem={0}".format(memory), "-J", split_name, "-o", err_name, "-x", "node02",
                script_path, dataset_name, bam_path, boundaries_path, out_pattern, status_path
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
    script_path = os.path.join(curr_path, "count_reads.py")
    # print("start") ####

    boundaries_path = "/agusevlab/DATA/ANNOTATIONS/gencode.v26lift37.annotation.patched_contigs.gtf"

    # Kellis 429 Partitioned
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    contigs = ["9", "10", "11", "12", "13", "14", "15", "17"]
    bam_path_kellis = os.path.join(data_path_kellis, "partitioned_429")
    names_kellis = os.listdir(bam_path_kellis)
    # print(names_kellis) ####
    out_pattern_base_kellis = os.path.join(data_path_kellis, "genes_429/{{0}}/bamdata/{{0}}_{0}.pickle")
    dispatch(script_path, "Kellis_429", bam_path_kellis, boundaries_path, names_kellis, contigs, out_pattern_base_kellis, 2000)
    # dispatch(script_path, "Kellis_429", bam_path_kellis, boundaries_path, names_kellis, contigs, out_pattern_base_kellis, 5000, fails_only=True)











import os
import sys
import pickle
import subprocess
import time
import numpy as np

def dispatch(script_path, names, data_dir, out_dir, memory, barcodes_map_path, fails_only=False):
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
            script_path, name, data_dir, out_dir, barcodes_map_path, status_path
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
    script_path = os.path.join(curr_path, "get_twas_inputs.py")

    data_path_kellis = "/agusevlab/awang/sc_kellis"
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")
    out_dir_kellis = os.path.join(data_path_kellis, "export_429")
    names = os.listdir(genes_dir_kellis)

    dispatch(script_path, names, genes_dir_kellis, out_dir_kellis, 2000, fails_only=False)

    # dispatch(script_path, names, genes_dir_kellis, out_dir_kellis, 2000, fails_only=True)






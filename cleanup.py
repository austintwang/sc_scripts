import os
import glob

def cleanup(genes_dir, wildcards):
	for i in os.listdir(genes_dir):
		for w in wildcards:
			rm = os.path.join(genes_dir, i, w)
			matches = glob.glob(rm)
			for m in matches:
				print(m)
				os.remove(m)

if __name__ == '__main__':
    # Kellis 429
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")
    wildcards = ["*.out"]

    cleanup(genes_dir_kellis, wildcards)
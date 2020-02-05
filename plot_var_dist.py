import os
import pickle
import numpy as np
import seaborn as sns

def plot_dist(in_path, out_path):
    with open(in_path, "rb") as in_file:
        dist_data = pickle.load(in_file)
    x = np.array(dist_data.values())
    print("\n".join(sorted(x))) ####

    sns.set()
    sns.distplot(x, rug=True)
    plt.savefig(out_path)
    plt.clf()

if __name__ == '__main__':
    data_dir = "/agusevlab/awang/sc_le/"
    in_path = os.path.join(data_dir, "variant_dist.pickle")
    in_path = os.path.join(data_dir, "dist.svg")
    bam_dist(in_path, out_path, 1000000)
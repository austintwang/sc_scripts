import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_dist(in_path, out_path):
    with open(in_path, "rb") as in_file:
        dist_data = pickle.load(in_file)
    x = np.log10(np.array(list(dist_data.values())))
    # print(dist_data) ##s##
    # print(x) ####
    # print("\n".join(map(str, sorted(x)))) ####
    print(np.count_nonzero(x >= 1)) ####
    print(np.count_nonzero(x >= 2)) ####

    sns.set()
    sns.distplot(x, norm_hist=False, kde=False)
    plt.savefig(out_path)
    plt.clf()

if __name__ == '__main__':
    data_dir = "/agusevlab/awang/sc_le/"
    in_path = os.path.join(data_dir, "variant_dist.pickle")
    out_path = os.path.join(data_dir, "dist.svg")
    plot_dist(in_path, out_path)
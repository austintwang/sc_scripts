import os
import pickle
import h5py

def load_hdf5(in_file):
    with h5py.File(in_file, "r") as hdf_in:
        barcodes = []
        for i in hdf_in["obs"]:
            if isinstance(i[0], bytes):
                entry = i[0].decode("utf-8")
            else:
                entry = i[0]

            # print(i[0].decode("utf-8").split("-")) ####
            try:
                barcodes.append(entry.split("-")[0])
            except Exception as e:
                print(entry) ####
                raise e

    return barcodes

def get_barcodes(in_dir, out_path):
    in_files = [os.path.join(in_dir, i) for i in os.listdir(in_dir) if i.endswith(".h5ad")]
    barcodes = []
    for i in in_files:
        barcodes.extend(load_hdf5(i))
    barcodes.sort()
    for ind, val in enumerate(barcodes): ####
        if val == barcodes[ind-1]:
            print(val)
    with open(out_path, "wb") as out_file:
        pickle.dump(barcodes, out_file)

if __name__ == '__main__':
    data_dir = "/agusevlab/awang/sc_le/"
    in_dir = os.path.join(data_dir, "genotypes")
    out_path = os.path.join(data_dir, "barcodes.pickle")
    get_barcodes(in_dir, out_path)
import os
import pickle
import h5py

# def load_hdf5(in_file):
#     with h5py.File(in_file, "r") as hdf_in:
#         barcodes = []
#         for i in hdf_in["obs"]:
#             if isinstance(i[0], bytes):
#                 entry = i[0].decode("utf-8")
#             else:
#                 entry = i[0]

#             # print(i[0].decode("utf-8").split("-")) ####
#             try:
#                 barcodes.append(entry.split("-")[0])
#             except Exception as e:
#                 print(entry) ####
#                 raise e

#     return barcodes

def load_barcodes(in_path):
    barcodes = {}
    with open(in_path, "r") as in_file:
        header = next(in_file)
        cols = header.split(",")
        barcode_idx = cols.index("index")
        well_idx = cols.index("well")
        for line in in_file:
            entries = line.split(",")
            bc = entries[barcode_idx].split("-")[0]
            barcodes.setdefault(entries[well_idx], []).append(bc)

    return barcodes


def get_barcodes(in_path, out_path):
    barcodes = load_barcodes(in_path)
    for v in barcodes.values():
        v.sort()
    with open(out_path, "wb") as out_file:
        pickle.dump(barcodes, out_file)

if __name__ == '__main__':
    data_dir = "/agusevlab/awang/sc_le/"
    in_path = os.path.join(data_dir, "adata_obs.V6.1.txt")
    out_path = os.path.join(data_dir, "barcodes.pickle")
    get_barcodes(in_path, out_path)
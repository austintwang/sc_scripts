import os
import pickle
import pysam
import numpy as np

def bam_dist(in_path, out_path, max_rows):
    dist = {}
    with pysam.AlignmentFile(in_path, "rb") as in_file:
        for ind, line in enumerate(in_file):
            if ind >= max_rows:
                break

            try:
                wasp_pass = line.get_tag("vW")
                if wasp_pass != 1:
                    continue
            except KeyError:
                continue

            try:
                intersects = line.get_tag("vG")
                # print(var) ####
                for var in intersects:
                    dist.setdefault(var, 0)
                    dist[var] += 1
            except KeyError:
                continue

            # if line.startswith("@"):
            #     continue

            # cols = line.split()
            # wasp_pass = False
            # var = None
            # for col in cols:
            #     if col.startswith("vW"):
            #         if col.split(":")[-1] == "1":
            #             wasp_pass = True
            #     if col.startswith("vG"):
            #         var = col
            #     if (var is not None) and wasp_pass:
            #         dist.setdefault(var, 0)
            #         dist[var] += 1
            #         break


    with open(out_path, "wb") as out_file:
        pickle.dump(dist, out_file)


if __name__ == '__main__':
    data_dir = "/agusevlab/awang/sc_le/"
    in_path = os.path.join(data_dir, "processed", "YE_7-19-1", "YE_7-19-1Aligned.sortedByCoord.out.bam")
    out_path = os.path.join(data_dir, "variant_dist.pickle")
    bam_dist(in_path, out_path, np.inf)


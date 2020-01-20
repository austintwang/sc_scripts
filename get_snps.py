import os
import pickle
import pysam

def load_vcf(in_path):
    vcf_in = pysam.VariantFile(in_path)
    snps = []
    for rec in vcf_in.fetch():
        entry = (rec.contig, rec.pos,)
        snps.append(entry)
    return snps

def get_snps(in_path, out_path):
    snps = load_vcf(in_path)
    snps.sort(key=lambda x: x[1])
    snps.sort(key=lambda x: x[0])
    for ind, val in enumerate(snps): ####
        if val == snps[ind-1]:
            print(val)
    with open(out_path, "wb") as out_file:
        pickle.dump(snps, out_file)

if __name__ == '__main__':
    data_dir = "/agusevlab/awang/sc_le/"
    in_path = os.path.join(data_dir, "HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz")
    out_path = os.path.join(data_dir, "markers.pickle")
    get_snps(path, out_path)
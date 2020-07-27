import os
import pysam

def make_star_vcf(in_path, out_path, min_af):
    with open(in_path, "r") as in_file, open(out_path, "w") as out_file:
        for line in in_file:
            if line.startswith("##"):
                continue
            cols = line.strip().split()
            if cols[0] == "#CHROM":
                joined = '\t'.join(cols)
                out_file.write(f"{joined}\tFORMAT\tVCF\n")
            else:
                fields = cols[-1].split(";")
                for f in fields:
                    if f.startswith("AF="):
                        af = float(f.split("=")[1])
                if af >= min_af:
                    joined = '\t'.join(cols)
                    out_file.write(f"{joined}\tGT\t0/1\n")

if __name__ == '__main__':
    genotypes_dir = "/agusevlab/awang/sc_le/genotypes/"
    if not os.path.exists(genotypes_dir):
        os.makedirs(genotypes_dir)
    in_path = os.path.join(genotypes_dir, "HRC.r1-1.GRCh37.wgs.mac5.sites.vcf")
    out_path = "/agusevlab/awang/sc_data/HRC.r1-1.GRCh37.wgs.mac5.maf05.sites.vcf"
    make_star_vcf(in_path, out_path, 0.05)

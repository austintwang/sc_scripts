import os
import pysam

def convert(in_path, out_path):
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    vcf_in = pysam.VariantFile(in_path)
    with open(out_path, "w") as bam_out:
        for rec in vcf_in.fetch():
            line = "{0}\t{1}\t{2}\n".format(rec.contig, rec.pos, rec.pos+1)
            bam_out.write(line)

if __name__ == '__main__':
    genotypes_dir = "/agusevlab/awang/sc_le/genotypes/"
    in_path = os.path.join(genotypes_dir, "HRC.r1-1.GRCh37.wgs.mac5.sites.vcf.gz")
    out_path = os.path.join(genotypes_dir, "hrc_sites.bed")
    convert(in_path, out_path)

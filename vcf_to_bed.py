import os
import pysam

def convert(in_path, out_path):
    contig_freqs = {} ####
    vcf_in = pysam.VariantFile(in_path)
    with open(out_path, "w") as bam_out:
        for rec in vcf_in.fetch():
            line = "{0}\t{1}\t{2}\n".format(rec.contig, rec.pos, rec.pos+1)
            bam_out.write(line)
            contig_freqs.setdefault(rec.contig, 0) ####
            contig_freqs[rec.contig] += 1 ####

if __name__ == '__main__':
    genotypes_dir = "/agusevlab/awang/sc_le/genotypes/"
    if not os.path.exists(genotypes_dir):
        os.makedirs(genotypes_dir)
    in_path = os.path.join(genotypes_dir, "HRC.r1-1.GRCh37.wgs.mac5.sites.vcf")
    out_path = os.path.join(genotypes_dir, "hrc_sites.bed")
    convert(in_path, out_path)

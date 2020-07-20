#!/usr/bin/env python3

import os
import sys
import pickle
import subprocess
import numpy as np

def count_bam(bam_path):
    contig_counts = {} ####
    req_tags = set(["vW", "vA", "vG", "CB", "RG"])
    args = ["/agusevlab/awang/samtools/bin/samtools", "view", bam_path]
    # print(" ".join(args)) ####
    with subprocess.Popen(" ".join(args), shell=True, stdout=subprocess.PIPE, bufsize=1, text=True) as p:
        for line in p.stdout:
            # print(line) ####
            cols = line.split("\t")
            chromosome = cols[2]
            start = int(cols[3])
            tag_data = {}
            for tag in cols[11:]:
                tag_name = tag[:2]
                if tag_name in req_tags:
                    tag_data[tag_name] = tag[2:]

            # contig_counts.setdefault(chromosome, 0) ####
            # contig_counts[chromosome] += 1 ####
            # if start % 10000 == 0:
            #     print(contig_counts) ####

            # wasp_pass = tag_data.get("vW")
            # if (wasp_pass is None) or int(wasp_pass[-1]) != 1:
            #     continue

            if len(chromosome) <= 2:
                contig_counts.setdefault(chromosome, 0) ####
                contig_counts[chromosome] += 1 ####
                if start % 100000 == 0:
                    print(contig_counts) ####

            # # print(line.get_tag("vA")) ####
            # genotype_raw = tag_data.get("vA", None)
            # if genotype_raw is None or len(genotype_raw) < 6:
            #     continue
            # genotype = int(genotype_raw[5]) - 1
            # # print(genotype, genotype_raw) ####
            # if not (genotype == 0 or genotype == 1):
            #     continue

            # barcode_raw = tag_data.get("CB")
            # if barcode_raw is None or len(barcode_raw) < 4:
            #     continue
            # barcode = barcode_raw[3:].split("-")[0]
            # # print(barcode, barcode_raw) ####

            # well_raw = tag_data.get("RG")
            # if well_raw is None or len(well_raw) < 4:
            #     continue
            # well = well_raw[3:].split(":")[0]
            # cell = barcode, well
            # # print(well, well_raw) ####

            # intersects_raw = tag_data.get("vG")
            # if intersects_raw is None or len(intersects_raw) < 6:
            #     continue
            # try:
            #     intersects = list(map(int, intersects_raw[5:].split(",")))
            #     # print(intersects, intersects_raw) ####
            # except ValueError:
            #     continue 

if __name__ == '__main__':
    bam_path = "/agusevlab/awang/sc_kellis/PFC_bam_files/SM-J2CJJ_SK-3TQI_pfc/D19-8353/outs/possorted_genome_bam.bam"
    count_bam(bam_path)
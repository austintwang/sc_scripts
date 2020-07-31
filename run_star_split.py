import os
import sys
import pickle
import subprocess
import numpy as np

def format_command(job_name, contig, readcmd, bam_path, bed_path, vcf_path, genome_path, boundaries_path, whitelist_path, out_prefix, paired, memory):
    # threads = str(min(64, 400 // (1400000 // memory)))
    # threads = str(24) ####
    threads = str(1)
    star_cmd = [
        "STAR",
        "--runMode", "alignReads",
        "--readFilesType", "SAM {0}".format("PE" if paired else "SE"),
        "--readFilesCommand", readcmd, str(contig), bed_path,
        "--outFilterMultimapNmax", "1",
        "--outFilterMatchNmin", "35",
        "--limitBAMsortRAM", str(int((memory - 6000) * 1e6)),
        "--runThreadN", threads,
        "--quantMode", "GeneCounts",
        "--twopassMode", "Basic",
        "--outFileNamePrefix", out_prefix,
        "--genomeDir", genome_path,
        "--sjdbGTFfile", boundaries_path,
        "--waspOutputMode", "SAMtag",
        "--varVCFfile", vcf_path,
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--soloCBwhitelist", whitelist_path,
        "--soloType", "Droplet",
        "--readFilesIn", bam_path,
        "--outSAMattributes", "NH", "HI", "AS", "nM", "vW", "vG", "vA", "CR", "CY", "UR", "UY",
        # "--outStd", "SAM"
    ]

    err_name = out_prefix + "_%j.out"
    cmd = [
        "sbatch",
        "--mem={0}".format(memory),
        "-c", threads,
        "-J",
        job_name,
        "-o",
        err_name,
        "-x", "node02,node04,node06,node07,node12,node13,node14,node21",
        "--wrap='{0}'".format(" ".join(star_cmd)) 
    ]

    print(" ".join(cmd))

    return cmd

def dispatch_star(bam_map, vcf_map, bed_map, contigs, readcmd, genome_path, boundaries_path, whitelist_path, out_path_base, memory, paired=False, selection=None):
    if selection is not None:
        bam_map = {k: v for k, v in bam_map.items() if k in selection}
        vcf_map = {k: v for k, v in vcf_map.items() if k in selection}
        bed_map = {k: v for k, v in bed_map.items() if k in selection}

    jobs = []
    for k, v in bam_map.items():
        vcf_path = vcf_map[k]
        bed_path = bed_map[k]
        for c in contigs:
            out_path = os.path.join(out_path_base, k, c)
            if not os.path.exists(out_path):
                os.makedirs(out_path)
            out_prefix = os.path.join(out_path, f"{k}_{c}")
            cmd = format_command(k, c, readcmd, v, bed_path, vcf_path, genome_path, boundaries_path, whitelist_path, out_prefix, paired, memory)
            jobs.append(cmd)

    # print(" & ".join([" ".join(cmd) for cmd in jobs])) ####
    with open("exec.sh", "w") as script_file:
        script_file.write("#!/bin/bash\n") ####
        script_file.writelines([(" ".join(cmd) + "\n") for cmd in jobs]) ####

    subprocess.run("./exec.sh", stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # timeout = "sbatch: error: Batch job submission failed: Socket timed out on send/recv operation"
    # for i in jobs:
    #     while True:
    #         try:
    #             submission = subprocess.run(i, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #             print(str(submission.stdout, 'utf-8').rstrip())
    #             break
    #         except subprocess.CalledProcessError as e:
    #             # print(e.stdout) ####
    #             err = str(e.stderr, 'utf-8').rstrip()
    #             print(err)
    #             if err == timeout:
    #                 print("Retrying Submit")
    #                 continue
    #             else:
    #                 raise e

def get_failed_jobs(names, out_path_base):
    fails = set()
    for i in names:
        out_bam_path = os.path.join(out_path_base, i, i + "Aligned.sortedByCoord.out.bam")
        if not os.path.isfile(out_bam_path) or os.path.getsize(out_bam_path) < 1e5:
            fails.add(i)
    return fails

if __name__ == '__main__':
    genome_path = "/cluster/agusevlab/awang/STAR_hg19/"
    boundaries_path = "/agusevlab/DATA/ANNOTATIONS/gencode.v26lift37.annotation.patched_contigs.gtf"
    whitelist_path = "/agusevlab/awang/sc_data/737K-august-2016.txt"
    vcf_hrc = "/agusevlab/DATA/ANNOTATIONS/HRC.r1-1.GRCh37.wgs.mac5.maf05.sites.vcf"
    bed_hrc = "/agusevlab/awang/sc_le/genotypes/hrc_sites.bed"
    contigs = [str(i) for i in range(1, 23)]
    curr_path = os.path.abspath(os.path.dirname(__file__))
    readcmd = os.path.join(curr_path, "bam_stream.py")


    # Kellis 429
    kellis_path_base = "/agusevlab/awang/sc_kellis"
    bam_path_kellis = os.path.join(kellis_path_base, "PFC_bam_files")
    bam_map_kellis_429 = {}
    vcf_map_kellis_429 = {}
    bed_map_kellis_429 = {}
    with open(os.path.join(kellis_path_base, "Bam_paths_432_PFC_HM_Austin.csv")) as bam_data:
        next(bam_data)
        for line in bam_data:
            cols = line.strip().split(",")
            bam_map_kellis_429[cols[1]] = os.path.join(bam_path_kellis, cols[2].lstrip("/"), cols[3].lstrip("/"))
            vcf_map_kellis_429[cols[1]] = vcf_hrc
            bed_map_kellis_429[cols[1]] = bed_hrc

    out_path_base_kellis_429 = os.path.join(kellis_path_base, "partitioned_429")
    # print(bam_map_kellis_429) ####

    dispatch_star(
        bam_map_kellis_429, vcf_map_kellis_429, bed_map_kellis_429, contigs, readcmd, genome_path, boundaries_path, whitelist_path, out_path_base_kellis_429, 20000
    )

    # fail_kellis_429 = get_failed_jobs(bam_map_kellis_429.keys(), out_path_base_kellis_429)
    # dispatch_star(
    #     bam_map_kellis_429, vcf_map_kellis_429, bed_map_kellis_429, genome_path, boundaries_path, whitelist_path, out_path_base_kellis_429, 260000, selection=fail_kellis_429
    # )



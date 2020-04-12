import os
import sys
import pickle
import subprocess
import numpy as np

def format_command(job_name, bam_path, bed_path, vcf_path, genome_path, boundaries_path, whitelist_path, out_prefix, paired, memory):
    threads = str(min(64, 400 // (1400000 // memory)))
    star_cmd = [
        "STAR",
        "--runMode", "alignReads",
        "--readFilesType", "SAM {0}".format("PE" if paired else "SE"),
        "--readFilesCommand", "samtools", "view", "-h", "-L", vcf_path,
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
        "-x", "node02",
        "--wrap='{0}'".format(" ".join(star_cmd)) 
    ]

    # print(" ".join(cmd))

    return cmd

def dispatch_star(bam_map, vcf_map, bed_map, genome_path, boundaries_path, whitelist_path, out_path_base, memory, paired=False, selection=None):
    if selection is not None:
        bam_map = {k: v for k, v in bam_map.items() if k in selection}
        vcf_map = {k: v for k, v in vcf_map.items() if k in selection}
        bed_map = {k: v for k, v in bed_map.items() if k in selection}

    jobs = []
    for k, v in bam_map.items():
        out_path = os.path.join(out_path_base, k)
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        out_prefix = os.path.join(out_path, k)
        vcf_path = vcf_map[k]
        bed_path = bed_map[k]
        cmd = format_command(k, v, bed_path, vcf_path, genome_path, boundaries_path, whitelist_path, out_prefix, paired, memory)
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
    whitelist_path = "/agusevlab/DATA/SCRNA/737K-august-2016.txt"
    vcf_hrc = "/agusevlab/DATA/ANNOTATIONS/HRC.r1-1.GRCh37.wgs.mac5.maf05.sites.vcf"
    bed_hrc = "/agusevlab/awang/sc_le/genotypes/hrc_sites.bed"

    # Ye lab (except "flare" bams)
    bam_path_ye = "/agusevlab/awang/sc_le/bam/"
    ye_non_flare = {
        "immvar_8_31-1" : "immvarYE_0831_1.bam.1",
        "immvar_8_31-2" : "immvarYE_0831_2.bam.1",
        "immvar_8_31-3" : "immvarYE_0831_3.bam.1",
        "immvar_8_31-4" : "immvarYE_0831_4.bam.1",
        "immvar_9_07-1" : "immvarYE_0907_1.bam.1",
        "immvar_9_07-2" : "immvarYE_0907_2.bam.1",
        "immvar_9_07-3" : "immvarYE_0907_3.bam.1",
        "immvar_9_07-4" : "immvarYE_0907_4.bam.1",
        "immvar_8_30-1" : "immvarYE_8_30_1.bam.1",
        "immvar_8_30-2" : "immvarYE_8_30_2.bam.1",
        "immvar_8_30-3" : "immvarYE_8_30_3.bam.1",
        "immvar_8_30-4" : "immvarYE_8_30_4.bam.1",
        "YE110-1" : "YE110_1.bam.1",
        "YE110-2" : "YE110_2.bam.1",
        "YE110-3" : "YE110_3.bam.1",
        "YE110-4" : "YE110_4.bam.1",
        "YE_7-13-1" : "YE_7_13_1.bam.1",
        "YE_7-13-2" : "YE_7_13_2.bam.1",
        "YE_7-13-3" : "YE_7_13_3.bam.1",
        "YE_7-13-4" : "YE_7_13_4.bam.1",
        "YE_7-19-1" : "YE_7_19_1.bam.1",
        "YE_7-19-2" : "YE_7_19_2.bam.1",
        "YE_7-19-3" : "YE_7_19_3.bam.1",
        "YE_7-19-4" : "YE_7_19_4.bam.1",
        "YE_7-20-1" : "YE_7_20_1.bam.1",
        "YE_7-20-2" : "YE_7_20_2.bam.1",
        "YE_7-20-3" : "YE_7_20_3.bam.1",
        "YE_7-20-4" : "YE_7_20_4.bam.1",
        "YE_7-26-1" : "YE_7_26_1.bam.1",
        "YE_7-26-2" : "YE_7_26_2.bam.1",
        "YE_7-26-3" : "YE_7_26_3.bam.1",
        "YE_7-26-4" : "YE_7_26_4.bam.1",
        "YE_8-16-1" : "YE_8_16_1.bam.1",
        "YE_8-16-2" : "YE_8_16_2.bam.1",
        "YE_8-16-3" : "YE_8_16_3.bam.1",
        "YE_8-16-4" : "YE_8_16_4.bam.1",
        "YE_8-17-1" : "YE_8_17_1.bam.1",
        "YE_8-17-2" : "YE_8_17_2.bam.1",
        "YE_8-17-3" : "YE_8_17_3.bam.1",
        "YE_8-17-4" : "YE_8_17_4.bam.1",
        "YE_8-2-1" : "YE_8_2_1.bam.1",
        "YE_8-2-2" : "YE_8_2_2.bam.1",
        "YE_8-23-1" : "YE_8_23_1.bam.1",
        "YE_8-23-2" : "YE_8_23_2.bam.1",
        "YE_8-23-3" : "YE_8_23_3.bam.1",
        "YE_8-23-4" : "YE_8_23_4.bam.1",
        "YE_8-2-3" : "YE_8_2_3.bam.1",
        "YE_8-2-4" : "YE_8_2_4.bam.1",
        "YE_8-3-1" : "YE_8_3_1.bam.1",
        "YE_8-3-2" : "YE_8_3_2.bam.1",
        "YE_8-3-3" : "YE_8_3_3.bam.1",
        "YE_8-3-4" : "YE_8_3_4.bam.1",
        "YE_8-9-3" : "YE_8_9_3.bam.1",
        "YE_8-9-4" : "YE_8_9_4.bam.1",
    }
    bam_map_ye_nf = {k: os.path.join(bam_path_ye, v) for k, v in ye_non_flare.items()}
    vcf_map_ye_nf = {k: vcf_hrc for k in ye_non_flare.keys()}
    bed_map_ye_nf = {k: bed_hrc for k in ye_non_flare.keys()}
    out_path_base_ye_nf = "/agusevlab/awang/sc_le/processed"
    # dispatch_star(bam_map_ye_nf, vcf_map_ye_nf, bed_map_ye_nf, genome_path, boundaries_path, whitelist_path, out_path_base_ye_nf, 10000)

    # Clean up Ye
    # fail_ye_nf = get_failed_jobs(ye_non_flare.keys(), out_path_base_ye_nf)
    # dispatch_star(
    #     bam_map_ye_nf, vcf_map_ye_nf, bed_map_ye_nf, genome_path, boundaries_path, whitelist_path, out_path_base_ye_nf, 300000, selection=fail_ye_nf
    # )

     # : "flare1_1.bam.1",
     # : "flare1_2.bam.1",
     # : "flare2_1.bam.1",
     # : "flare2_2.bam.1",
     # : "flare3_1.bam.1",
     # : "flare3_2.bam.1",
     # : "flare3_3.bam.1",
     # : "flare3_4.bam.1",
     # : "flare4_1.bam.1",
     # : "flare4_2.bam.1",
     # : "flare4_3.bam.1",
     # : "flare4_4.bam.1",


    # Kellis 48
    # kellis_path_base = "/agusevlab/awang/sc_kellis"
    # bam_path_kellis = os.path.join(kellis_path_base, "121719_10xdata")
    # kellis_48 = {i: "{0}/{0}.bam".format(i) for i in os.listdir(bam_path_kellis)}
    # bam_map_kellis_48 = {k: os.path.join(bam_path_kellis, v) for k, v in kellis_48.items()}
    # vcf_map_kellis_48 = {k: vcf_hrc for k in kellis_48.keys()}
    # bed_map_kellis_48 = {k: bed_hrc for k in kellis_48.keys()}
    # out_path_base_kellis_48 = os.path.join(kellis_path_base, "processed")
    # dispatch_star(
    #     bam_map_kellis_48, vcf_map_kellis_48, bed_map_kellis_48, genome_path, boundaries_path, whitelist_path, out_path_base_kellis_48, 20000
    # )

    # # Clean up Kellis
    # fail_kellis_48 = get_failed_jobs(kellis_48.keys(), out_path_base_kellis_48)
    # dispatch_star(
    #     bam_map_kellis_48, vcf_map_kellis_48, bed_map_kellis_48, genome_path, boundaries_path, whitelist_path, out_path_base_kellis_48, 60000, selection=fail_kellis_48
    # )

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

    out_path_base_kellis_429 = os.path.join(kellis_path_base, "processed_429")
    # print(bam_map_kellis_429) ####

    # dispatch_star(
    #     bam_map_kellis_429, vcf_map_kellis_429, bed_map_kellis_429, genome_path, boundaries_path, whitelist_path, out_path_base_kellis_429, 20000
    # )

    fail_kellis_429 = get_failed_jobs(bam_map_kellis_429.keys(), out_path_base_kellis_429)
    dispatch_star(
        bam_map_kellis_429, vcf_map_kellis_429, bed_map_kellis_429, genome_path, boundaries_path, whitelist_path, out_path_base_kellis_429, 40000, selection=fail_kellis_429
    )



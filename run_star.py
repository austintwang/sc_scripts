import os
import sys
import pickle
import subprocess
import numpy as np

def format_command(bam_path, bed_path, vcf_path, genome_path, boundaries_path, whitelist_path, out_path):
	star_cmd = [
		"STAR",
		"--runMode", "alignReads",
		"--readFilesCommand", "samtools", "view", "-h", "-L", vcf_path,
		"--outFilterMultimapNmax", "1",
		"--outFilterMatchNmin", "35",
		"--quantMode", "GeneCounts",
		"--twopassMode", "Basic",
		"--outFileNamePrefix", out_path,
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

	cmd = [
		"sbatch",
		"--mem=50000",
		"--wrap='{0}'".format(" ".join(star_cmd))
	]

	return cmd

def dispatch_star(bam_map, vcf_map, bed_map, genome_path, boundaries_path, whitelist_path, out_path_base):
	jobs = []
	for k, v in bam_map.items():
		out_path = os.path.join(out_path_base, k, k)
		vcf_path = vcf_map[k]
		cmd = format_command(v, bed_path, vcf_path, genome_path, boundaries_path, whitelist_path, out_path)
		jobs.append(cmd)

	for i in jobs:
		while True:
			try:
				submission = subprocess.run(i, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				print(str(submission.stdout, 'utf-8').rstrip())
				break
			except subprocess.CalledProcessError as e:
				# print(e.stdout) ####
				err = str(e.stderr, 'utf-8').rstrip()
				print(err)
				if err == timeout:
					print("Retrying Submit")
					continue
				else:
					raise e

if __name__ == '__main__':
	genome_path = "/agusevlab/DATA/GENOMES/STAR_hg19/"
	boundaries_path = "/agusevlab/DATA/ANNOTATIONS/gencode.v26lift37.annotation.patched_contigs.gtf"
	whitelist_path = "/agusevlab/DATA/SCRNA/737K-august-2016.txt"

	# Ye lab (except "flare" bams)
	bam_path_ye = "/agusevlab/awang/sc_le/bam/"
	vcf_hrc = "/agusevlab/awang/sc_le/genotypes/HRC.r1-1.GRCh37.wgs.mac5.sites.vcf"
	bed_hrc = "/agusevlab/awang/sc_le/genotypes/hrc_sites.bed"
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
	dispatch_star(bam_map_ye_nf, vcf_map_ye_nf, bed_map_ye_nf, genome_path, boundaries_path, whitelist_path, out_path_base_ye_nf)

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




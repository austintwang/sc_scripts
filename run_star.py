#!/usr/bin/env python3

import os
import sys
import pickle
import subprocess
import numpy as np

def format_command(bam_path, bed_path, vcf_path, genome_path, boundaries_path, whitelist_path, log_path):
	cmd = [
		"STAR",
		"--runMode", "alignReads",
		"--readFilesCommand", "samtools", "view", "-h", "-L", vcf_path, bed_path,
		"--outFilterMultimapNmax", "1",
		"--outFilterMatchNmin", "35",
		"--quantMode", "GeneCounts",
		"--twopassMode", "Basic",
		"--outFileNamePrefix", log_path,
		"--genomeDir", genome_path,
		"--sjdbGTFfile", boundaries_path,
		"--waspOutputMode", "SAMtag",
		"--varVCFfile", vcf_path,
		"--outSAMtype", "BAM", "SortedByCoordinate",
		"--soloCBwhitelist", whitelist_path,
		"--soloType", "Droplet",
		"--readFilesIn", bam_path,
		"--outSAMattributes", "NH", "HI", "AS", "nM", "vW", "vG", "vA", "CR", "CY", "UR", "UY",
		"--outStd", "SAM"
	]
	return cmd

def star(
		well_name,
		markers_path,
		barcodes_path,
		bam_path, 
		bed_path, 
		vcf_path, 
		genome_path, 
		boundaries_path, 
		whitelist_path, 
		log_path,
		out_path
	):

	with open(markers_path, "rb") as markers_file:
		markers = pickle.load(markers_file)

	with open(barcodes_path, "rb") as barcodes_file:
		barcodes_all = pickle.load(barcodes_file)
	# print(barcodes_all.keys()) ####
	barcodes = barcodes_all[well_name]

	marker_map = {val: ind for ind, val in enumerate(markers)}
	barcode_map = {val: ind for ind, val in enumerate(barcodes)}

	# ref = np.zeros((len(markers), len(barcodes),), dtype='uint8')
	# alt = np.zeros((len(markers), len(barcodes),), dtype='uint8')

	cmd = format_command(bam_path, bed_path, vcf_path, genome_path, boundaries_path, whitelist_path, log_path)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	for line in proc.stdout:
		if line.startswith("@"):
			continue
		print(line) ####

if __name__ == '__main__':
	star(*sys.argv[1:])







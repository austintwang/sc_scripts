#!/usr/bin/env python3

import os
import sys
import pickle
import subprocess
import numpy as np

class ReadBuffer(object):
	def __init__(self, depth, cov_threshold, num_cells):
		self.depth = depth
		self.num_cells = num_cells

		self.buffer = depth * [None]
		self.buffer_data = {}
		self.markers = []
		self.marker_data = []
		self.pos = 0

	def add_read(self, marker, cell_idx, var_idx):
		if marker not in self.buffer_data:
			retire_marker = self.buffer[self.pos]
			if retire_marker is not None:
				retire_data = self.buffer_data.pop(retire_marker)
				if np.sum(retire_data) >= cov_threshold:
					self.markers.append(retire_marker)
					self.markers.append(retire_data)
			
			self.buffer[self.pos] = marker
			self.buffer_data[marker] = np.zeros((self.num_cells, 2,), dtype='uint16')

			self.pos = (self.pos + 1) % self.depth

		self.buffer_data[marker][cell_idx] += 1

	def purge_buffer(self):
		for i in range(self.pos, self.pos + self.depth):
			idx = i % self.depth
			marker = self.buffer[idx]
			if marker is not None:
				self.markers.append(marker)
				self.marker_data.append(buffer_data[marker][var_idx])

		self.buffer = depth * [None]
		self.buffer_data = {}
		self.markers = []
		self.marker_data = []
		self.pos = 0

def format_command(bam_path, bed_path, vcf_path, genome_path, boundaries_path, whitelist_path, log_path):
	cmd = [
		"STAR",
		"--runMode", "alignReads",
		"--readFilesCommand", "samtools", "view", "-h", "-L", vcf_path,
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
		samples_path,
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

	with open(samples_path, "rb") as samples_file:
		samples_all = pickle.load(samples_file)
	samples = samples_all[well_name]

	with open(barcodes_path, "rb") as barcodes_file:
		barcodes_all = pickle.load(barcodes_file)
	# print(barcodes_all.keys()) ####
	barcodes = barcodes_all[well_name]

	barcode_map = {val: ind for ind, val in enumerate(barcodes)}

	# ref = np.zeros((len(markers), len(barcodes),), dtype='uint8')
	# alt = np.zeros((len(markers), len(barcodes),), dtype='uint8')
	rbuf = ReadBuffer(10, 100 * len(samples), len(barcodes))

	cmd = format_command(bam_path, bed_path, vcf_path, genome_path, boundaries_path, whitelist_path, log_path)
	proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
	for line in proc.stdout:
		if line.startswith("@"):
			continue
		print(line) ####
		cols = line.split("\t")
		wasp_pass = True
		for c in cols:
			if c.startswith("vW"):
				wasp_status = int(c.split(":")[-1])
				if wasp_status != 1:
					wasp_pass = False
					break
			elif c.startswith("vG"):
				marker = c.split(":")[-1]
			elif c.startswith("CB"):
				cell_idx = barcode_map[c.split(":")[-1].split("-")[0]]
			elif c.startswith("vA"):
				var_idx = int(c.split(":")[-1]) - 1

		if wasp_pass == True:
			rbuf.add_read(marker, cell_idx, var_idx)

	rbuf.purge_buffer()

	with open(out_path, "wb") as out_file:
		pickle.dump(out_file, {"markers": rbuf.markers, "data": rbuf.marker_data})

if __name__ == '__main__':
	star(*sys.argv[1:])







#!/usr/bin/env python3

import os
import sys
import pickle
import subprocess
import numpy as np
import pysam

class ReadBuffer(object):
	def __init__(self, depth, cov_threshold):
		self.depth = depth
		self.buffer = depth * [None]
		self.buffer_data = {}
		self.markers = []
		self.marker_data = []
		self.out_data = {}
		self.pos = 0

	def add_read(self, chrm, pos, cell, genotype):
		marker = (chrm, pos)
		if marker not in self.buffer_data:
			retire_marker = self.buffer[self.pos]
			if retire_marker is not None:
				retire_data = self.buffer_data.pop(retire_marker)
				if np.sum(retire_data) >= cov_threshold:
					self.markers.append(retire_marker)
					self.marker_data.append(retire_data)
			
			self.buffer[self.pos] = marker
			self.buffer_data[marker] = {}
			self.pos = (self.pos + 1) % self.depth

		self.buffer_data[marker].setdefault(cell, np.zeros(2, dtype='uint16'))[genotype] += 1

	def purge(self):
		for i in range(self.pos, self.pos + self.depth):
			idx = i % self.depth
			marker = self.buffer[idx]
			if marker is not None:
				self.markers.append(marker)
				self.marker_data.append(buffer_data[marker][var_idx])

		self.buffer = depth * [None]
		self.buffer_data = {}
		self.pos = 0

	def get_data(self):
		for m, d in zip(markers, marker_data):
			chrm, pos = m
			entry = self.out_data.setdefault(chrm, [])
			entry.append([pos, d])

		return self.out_data


def count_bam(bam_path, readbuf):
	with pysam.AlignmentFile(bam_path, "rb") as bam_file:
		for ind, line in enumerate(in_file):
            try:
            	wasp_pass = line.get_tag("vW")
                if wasp_pass != 1:
                    continue

            	genotype = line.get_tag("vA") - 1
            	if not (genotype == 0 or genotype == 1):
            		continue

            	barcode = line.get_tag("CB")
            	well = line.get_tag("RG").split(":")[0]
            	cell = (barcode, well)
            	chromosome = line.reference_name
                intersects = line.get_tag("vG")

                for var in intersects:
                    readbuf.add_read(chromosome, var, cell, genotype)
            
            except KeyError:
                continue

    readbuf.purge() 

def count_reads(data_dir, name, buffer_depth, cov_threshold):
	bam_path = os.path.join(data_dir, "{0}Aligned.sortedByCoord.out.bam".format(name))
	out_path = os.path.join(data_dir, "{0}_counts.pickle")

	readbuf = ReadBuffer(buffer_depth, cov_threshold)
	count_bam(bam_path, readbuf)
	out_data = readbuf.get_data()

	with open(out_path, "rb") as out_file:
		pickle.dump(out_data, out_file)

if __name__ == '__main__':
	count_reads(*sys.argv[1:])







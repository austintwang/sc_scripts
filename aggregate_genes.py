#!/usr/bin/env python3

import os
import sys
import pickle
import subprocess
import numpy as np
import pysam
# import ncls

class IntervalBuffer(object):
	def __init__(self, intervals, keys):
		self.intervals = sorted(intervals)
		self.keymap = dict(zip(intervals, keys))

		self.idx = 0
		self.window = set([])

	def query(query_pos):
		# if self.intervals[self.idx][0] > query_pos:
		# 	raise ValueError("Query's position behind previous query's")

		while self.intervals[self.idx][0] <= query_pos:
			curr_interval = self.intervals[self.idx]
			if curr_interval[1] < query_pos:
				self.idx += 1
			else:
				self.window.add(curr_interval)
				self.idx += 1

		intersects = []
		retires = []
		for i in self.window:
			if i[1] >= query_pos:
				intersects.append(i)
			else:
				retires.append(i)

		for i in retires:
			self.window.remove(i)

		intersects_names = [self.keymap[i] for i in intersects]

		return intersects_names

	def reset(self):
		self.idx = 0
		self.window.clear()

class MarkerBuffer(object):
	def __init__(self, intervals, keys):
		self.intervals = sorted(intervals)
		self.keymap = dict(zip(intervals, keys))

		self.idx = 0
		self.window = set([])

	def query(self, query_pos):
		# if self.intervals[self.idx][0] > query_pos:
		# 	raise ValueError("Query's position behind previous query's")

		while self.intervals[self.idx][0] <= query_pos:
			curr_interval = self.intervals[self.idx]
			if curr_interval[1] < query_pos:
				self.idx += 1
			else:
				self.window.add(curr_interval)
				self.idx += 1

		intersects = []
		retires = []
		for i in self.window:
			if i[1] >= query_pos:
				intersects.append(i)
			else:
				retires.append(i)

		for i in retires:
			self.window.remove(i)

		intersects_names = [self.keymap[i] for i in intersects]

		return intersects_names

	def reset(self):
		self.idx = 0
		self.window.clear()

class GeneAggregator(object):
	def __init__(self, exons, barcode_map):
		self.buffer_inputs = {}
		for name, chrm, start, end in exons:
			intervals, keys = self.buffer_data.setdefault(chrm, [[],[]])
			intervals.append((start, end),)
			keys.append(name)
		self.buffers = {k, IntervalBuffer(*v) for k, v in self.buffer_inputs.items()}
		self.gene_data = {}
		self.barcode_map = barcode_map

	def add_data(bam_data):
		for k, v in bam_data.items():
			for pos, read_data in v:
				genes = self.buffers[k].query(pos)
				for g in genes:

					samples = self.gene_data.setdefault(g, {})




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







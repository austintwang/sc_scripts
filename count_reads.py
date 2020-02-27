#!/usr/bin/env python3

import os
import sys
import pickle
import subprocess
import numpy as np
import pysam

class ReadBuffer(object):
    def __init__(self, depth, marker_buf):
        self.depth = depth
        self.buffer = depth * [None]
        self.buffer_data = {}
        self.marker_buf = marker_buf
        self.pos = 0

    def add_read(self, chrm, posns, cell, genotype):
        markers = [(chrm, pos) for pos in posns]
        for marker in markers:
            if marker not in self.buffer_data:
                retire_marker = self.buffer[self.pos]
                if retire_marker is not None:
                    retire_data = self.buffer_data.pop(retire_marker)
                    self.marker_buf.add_marker(retire_marker, retire_data)
                
                self.buffer[self.pos] = marker
                self.buffer_data[marker] = {}
                self.pos = (self.pos + 1) % self.depth

            marker_data = self.buffer_data[marker].setdefault(cell, np.zeros(2, dtype='uint16'))
            marker_data += 1

    def purge(self):
        for i in range(self.pos, self.pos + self.depth):
            idx = i % self.depth
            marker = self.buffer[idx]
            if marker is not None:
                retire_data = self.buffer_data.pop(marker)
                if np.sum(retire_data[0]) >= cov_threshold:
                    self.marker_buf.add_marker(retire_marker, retire_data)

        self.marker_buf.purge()
        self.buffer = depth * [None]
        self.buffer_data = {}
        self.pos = 0


class MarkerBuffer(object):
    def __init__(self, depth, out_pattern, gene_finder):
        self.out_pattern = out_pattern
        self.depth = depth
        self.buffer = depth * [None]
        self.buffer_data = {}
        self.pos = 0
        self.gene_finder = gene_finder

    def add_marker(self, marker, data):
        genes = self.gene_finder.query()
        for g in genes:
            if g not in self.buffer_data:
                retire_gene = self.buffer[self.pos]
                if retire_gene is not None:
                    retire_data = self.buffer_data.pop(retire_gene)
                    self._retire(retire_gene, retire_data)
            
                self.buffer[self.pos] = g
                self.buffer_data[g] = {}
                self.pos = (self.pos + 1) % self.depth

            gene_data = self.buffer_data[g].setdefault(marker, {})
            gene_data[marker] = data

    def purge(self):
        for i in range(self.pos, self.pos + self.depth):
            idx = i % self.depth
            retire_gene = self.buffer[idx]
            if retire_gene is not None:
                retire_data = self.buffer_data.pop(retire_gene)
                self._retire(retire_gene, retire_data)
        
            self.buffer[self.pos] = g
            self.buffer_data[g] = {}
            self.pos = (self.pos + 1) % self.depth

        gene_data = self.buffer_data[g].setdefault(marker, {})
        gene_data[marker] = data


    def _retire(self, gene, data):
        out_path = self.out_pattern.format(gene)
        out_dir = os.path.dirname(out_path)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        with open(out_path, "wb") as out_file:
            pickle.dump(data, out_file)


class GeneFinder(object):
    def __init__(self, exons, contig_order):
        self.contig_map = {val: ind for ind, val in enumerate(contig_order)}
        self.exons = [[self.contig_map[i[0]]] + i[1:] for i in exons]
        self.intervals = sorted(self.exons, lambda x: (contig_map[x[0]], x[1]))
        self.idx = 0
        self.window = set([])

    def query(query_pos):
        while self.intervals[self.idx][1] <= query_pos[1]:
            curr_interval = self.intervals[self.idx]
            if curr_interval[0] > self.contig_map[query_pos[0]]:
                break
            if curr_interval[0] < self.contig_map[query_pos[0]]:
                self.idx += 1
            elif curr_interval[2] < query_pos[1]: 
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

        return [i[2] for i in intersects]

    def reset(self):
        self.idx = 0
        self.window.clear()

def count_bam(bam_path, exons, out_pattern):
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        contig_data = bam_file.header["SQ"]
        contig_order = [i["SN"] for i in contig_data]
        gene_finder = GeneFinder(exons, contig_order)
        markerbuf = MarkerBuffer(10, out_pattern, gene_finder)
        readbuf = ReadBuffer(10, markerbuf)

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
                readbuf.add_read(chromosome, intersects, cell, genotype)
            
            except KeyError:
                continue

    readbuf.purge()

def load_exons(boundaries_path):
    exons = []
    with open(boundaries_path, "r") as boundaries_file:
        for line in boundaries_file:
            if line.startswith("##"):
                continue
            data = line.split("\t")
            # print(data) ####
            if data[2] == "exon":
                contig = data[0]
                start = data[3]
                end = data[4]
                gene = data[-1].split(";")[0].split(" ")[1].strip("\"")
                exons.append([contig, start, end, gene])

    return exons

def count_reads(bam_path, boundaries_path, out_pattern):
    exons = load_exons(boundaries_path)
    count_bam(bam_path, exons, out_pattern)

if __name__ == '__main__':
    count_reads(*sys.argv[1:])







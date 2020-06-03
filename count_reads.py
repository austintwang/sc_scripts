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

    def add_read(self, chrm, start, posns, cell, genotype):
        if len(posns) == 0:
            return

        markers = [(chrm, start, True),]
        markers.extend((chrm, pos, False) for pos in posns)
        # print(markers) ####
        for i, marker in enumerate(markers):
            if marker not in self.buffer_data:
                retire_marker = self.buffer[self.pos]
                if retire_marker is not None:
                    retire_data = self.buffer_data.pop(retire_marker)
                    self.marker_buf.add_marker(retire_marker, retire_data, i==0)
                
                self.buffer[self.pos] = marker
                self.buffer_data[marker] = {}
                self.pos = (self.pos + 1) % self.depth

            marker_data = self.buffer_data[marker].setdefault(cell, np.zeros(2, dtype='uint16'))
            marker_data[genotype] += 1

    def purge(self):
        for i in range(self.pos, self.pos + self.depth):
            idx = i % self.depth
            marker = self.buffer[idx]
            if marker is not None:
                retire_data = self.buffer_data.pop(marker)
                self.marker_buf.add_marker(marker, retire_data, i==0)

        self.marker_buf.purge()
        self.buffer = self.depth * [None]
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

    def add_marker(self, marker, data, checkpoint):
        # print(marker) ####
        # print(data) ####
        genes = self.gene_finder.query(marker)
        for g in genes:
            if g not in self.buffer_data:
                retire_gene = self.buffer[self.pos]
                if retire_gene is not None:
                    retire_data = self.buffer_data.pop(retire_gene)
                    self._retire(retire_gene, retire_data)
            
                self.buffer[self.pos] = g
                self.buffer_data[g] = {}
                self.pos = (self.pos + 1) % self.depth

            if not marker[2]:
                self.buffer_data[g][marker[:2]] = data

    def purge(self):
        for i in range(self.pos, self.pos + self.depth):
            idx = i % self.depth
            retire_gene = self.buffer[idx]
            if retire_gene is not None:
                retire_data = self.buffer_data.pop(retire_gene)
                self._retire(retire_gene, retire_data)
        
        self.buffer[self.pos] = self.depth * [None]
        self.buffer_data = {}
        self.pos = 0


    def _retire(self, gene, data):
        out_path = self.out_pattern.format(gene)
        out_dir = os.path.dirname(out_path)
        try:
            os.makedirs(out_dir)
        except FileExistsError:
            pass

        with open(out_path, "wb") as out_file:
            pickle.dump(data, out_file)
        # print([(k, np.sum(np.stack(v.values()), axis=0)) for k, v in data.items()]) ####


class GeneFinder(object):
    def __init__(self, exons, contig_order):
        # print(exons[0]) ####
        self.contig_map = {val: ind for ind, val in enumerate(contig_order)}
        # print(self.contig_map) ####
        self.exons = [tuple([self.contig_map[i[0]]] + i[1:]) for i in exons]
        self.intervals = sorted(self.exons)
        self.idx = 0
        self.window = set([])
        self.idx_checkpoint = 0
        self.window_checkpoint = set([])
        # print(self.intervals[0], self.intervals[-1]) ####

    def query(self, query_pos):
        checkpoint = query_pos[2]
        if (
            (
                (self.intervals[self.idx][1] > query_pos[1]) 
                and (self.intervals[self.idx][0] == self.contig_map[query_pos[0]])
            )
            or (self.intervals[self.idx][0] > self.contig_map[query_pos[0]])
        ):
            print(self.intervals[self.idx], self.intervals[self.idx_checkpoint]) ####
            self.idx = self.idx_checkpoint
            self.window = self.window_checkpoint

        while True:
            if (self.idx + 1) == len(self.intervals):
                break
            next_interval = self.intervals[self.idx + 1]
            if next_interval[0] > self.contig_map[query_pos[0]]:
                break
            if next_interval[1] > query_pos[1]:
                break
            if next_interval[0] < self.contig_map[query_pos[0]]:
                self.idx += 1
            elif next_interval[2] < query_pos[1]: 
                self.idx += 1
            else:
                self.window.add(next_interval)
                self.idx += 1

        # while (self.intervals[self.idx][1] <= query_pos[1]) or (self.intervals[self.idx][2] < self.contig_map[query_pos[0]]):
        #     curr_interval = self.intervals[self.idx]
        #     if curr_interval[0] > self.contig_map[query_pos[0]]:
        #         break
        #     if curr_interval[0] < self.contig_map[query_pos[0]]:
        #         self.idx += 1
        #     elif curr_interval[2] < query_pos[1]: 
        #         self.idx += 1
        #     else:
        #         self.window.add(curr_interval)
        #         self.idx += 1

        intersects = []
        retires = []
        for i in self.window:
            if i[2] >= query_pos[1]:
                intersects.append(i)
            else:
                retires.append(i)

        for i in retires:
            self.window.remove(i)

        print(query_pos, checkpoint, intersects) ####

        if checkpoint:
            self.idx_checkpoint = self.idx
            self.window_checkpoint = self.window.copy()

        return [i[3] for i in intersects]

    def reset(self):
        self.idx = 0
        self.window.clear()

def get_readdata_ye(line):
    barcode = line.get_tag("CB").split("-")[0]
    well = line.get_tag("RG").split(":")[0]
    return barcode, well

def get_readdata_kellis(line):
    # print(line.get_tag("RG")) ####
    data = line.get_tag("RG").split(":")
    barcode = data[0].split("-")[0]
    well = data[1].split("_")[0]
    return barcode, well

def get_readdata_kellis_429(line):
    barcode = line.get_tag("CB").split("-")[0]
    well = line.get_tag("RG").split(":")[0]
    return barcode, well   

def count_bam(bam_path, exons, readdata_fn, out_pattern, parse_manual):
    with pysam.AlignmentFile(bam_path, "rb") as bam_file:
        contig_data = bam_file.header["SQ"]
        # contigs = {i["SN"]: i["LN"] for i in contig_data}
        # contig_order = sorted(contigs.keys(), key=contigs.get)
        # print(contig_order) ####
        contig_order = [i["SN"] for i in contig_data]
        gene_finder = GeneFinder(exons, contig_order)
        markerbuf = MarkerBuffer(10, out_pattern, gene_finder)
        readbuf = ReadBuffer(10, markerbuf)

        if not parse_manual:
            for line in bam_file:
                # print(line) ####
                # print(line.reference_name) ####
                try:
                    wasp_pass = line.get_tag("vW")
                    if wasp_pass != 1:
                        continue

                    # print(line.get_tag("vA")) ####
                    genotype = line.get_tag("vA")[0] - 1
                    if not (genotype == 0 or genotype == 1):
                        continue

                    cell = readdata_fn(line)

                    chromosome = line.reference_name
                    intersects = line.get_tag("vG")
                    readbuf.add_read(chromosome, intersects, cell, genotype)
                
                except KeyError:
                    continue

    if parse_manual:
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

                wasp_pass = tag_data.get("vW")
                if (wasp_pass is None) or int(wasp_pass[-1]) != 1:
                    continue

                # print(line.get_tag("vA")) ####
                genotype_raw = tag_data.get("vA", None)
                if genotype_raw is None or len(genotype_raw) < 6:
                    continue
                genotype = int(genotype_raw[5]) - 1
                # print(genotype, genotype_raw) ####
                if not (genotype == 0 or genotype == 1):
                    continue

                barcode_raw = tag_data.get("CB")
                if barcode_raw is None or len(barcode_raw) < 4:
                    continue
                barcode = barcode_raw[3:].split("-")[0]
                # print(barcode, barcode_raw) ####

                well_raw = tag_data.get("RG")
                if well_raw is None or len(well_raw) < 4:
                    continue
                well = well_raw[3:].split(":")[0]
                cell = barcode, well
                # print(well, well_raw) ####

                intersects_raw = tag_data.get("vG")
                if intersects_raw is None or len(intersects_raw) < 6:
                    continue
                try:
                    intersects = list(map(int, intersects_raw[5:].split(",")))
                    # print(intersects, intersects_raw) ####
                except ValueError:
                    continue 

                readbuf.add_read(chromosome, start, intersects, cell, genotype)

    readbuf.purge()

    # print(chromosome, intersects) ####

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
                start = int(data[3])
                end = int(data[4])
                gene = data[-1].split(";")[0].split(" ")[1].strip("\"")
                exons.append([contig, start, end, gene])

    return exons

def count_reads(dataset_name, bam_path, boundaries_path, out_pattern, status_path):
    with open(status_path, "w") as status_file:
        status_file.write("")
    exons = load_exons(boundaries_path)
    if dataset_name == "Ye":
        readdata_fn = get_readdata_ye
        parse_manual = False
    elif dataset_name == "Kellis":
        readdata_fn = get_readdata_kellis
        parse_manual = False
    elif dataset_name == "Kellis_429":
        readdata_fn = get_readdata_kellis_429
        parse_manual = True
    count_bam(bam_path, exons, readdata_fn, out_pattern, parse_manual)
    with open(status_path, "w") as status_file:
        status_file.write("Complete")

if __name__ == '__main__':
    count_reads(*sys.argv[1:])







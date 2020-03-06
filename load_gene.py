#!/usr/bin/env python3

import os
import sys
import pickle
import numpy as np
import pysam

def read_vcf(vcf_path, contig, start, end):
    genotypes_list = []
    markers = [] 
    marker_ids = []
    vcf = pysam.VariantFile(vcf_path)
    samples = vcf.samples
    for record in vcf.fetch(contig, start, end):
        add_marker = True
        record_gens = np.zeros((len(samples), 2,),)
        for ind, sample in enumerate(samples)
            val = record.samples[sample]
            if not val.phased:
                add_marker = False
                break
            record_gens[ind,:] = np.array(val.alleles)
        if add_marker:
            genotypes_list.append(record_gens)
            markers.append((chrom, record.pos,),)
            marker_ids.append(record.id) 

    genotypes = np.stack(genotypes_list, axis=1)
    return genotypes, samples, markers, marker_ids

def add_data(agg_counts, var_data, cell_map, genotypes, sample_gen_map, marker_gen_map):
    for var, cells in var_data.items():
        for cell, counts in cells.items():
            cell_agg = agg_counts.setdefault(cell, np.array([0,0,0]))
            cell_agg[2] += np.sum(counts)
            if not (var in markers):
                continue
            sample_idx = sample_gen_map[cell_map[k]]
            marker_idx = marker_gen_map[var]
            gen = genotypes[sample, var, :]
            if np.sum(gen) == 1:
                cell_agg += counts[gen]

def load_gene(gene_name, radius, gene_dir, vcf_path, barcodes_map_path, boundaries_map_path, tss_map_path):
    with open(barcodes_map_path, "rb") as barcodes_map_file:
        barcodes_map = pickle.load(barcodes_map_file)
    with open(boundaries_map_path, "rb") as boundaries_map_file:
        boundaries_map = pickle.load(boundaries_map_file)
    with open(tss_map_path, "rb") as tss_map_file:
        tss_map = pickle.load(tss_map_file)

    contig, start, end = boundaries_map[gene_name]
    genotypes, samples, markers, marker_ids = read_vcf(vcf_path, contig, start, end)
    sample_gen_map = dict([(val, ind) for ind, val in enumerate(samples)])
    markers = dict([(val, ind) for ind, val in enumerate(markers)])

    contig, tss_pos = tss_map[gene_name.split(".")[0]]
    genotypes_nc, samples_nc, markers_nc, marker_ids_nc = read_vcf(
        vcf_path, contig, tss_pos - radius, tss_pos + radius + 1
    )
    sample_gen_map_nc = dict([(val, ind) for ind, val in enumerate(samples_nc)])
    markers_nc = dict([(val, ind) for ind, val in enumerate(markers_nc)])
    
    agg_counts = {}
    var_data_paths = os.listdir(os.path.join(gene_dir, "bamdata")) 
    for path in var_data_paths:
        with open(path, "rb") as var_file:
            var_data = pickle.load(var_file)
            add_data(agg_counts, var_data, cell_map, genotypes, sample_gen_map, marker_gen_map)

    out_data = {
        "name": gene_name, 
        "genotypes": genotypes_nc, 
        "samples": samples_nc, 
        "markers": markers_nc, 
        "marker_ids": marker_ids_nc,
        "cell_counts": agg_counts
    }
    out_path = os.path.join(gene_dir, "gene_data.pickle")
    with open(out_path, "wb") as out_file:
        pickle.dump(out_data, out_file)

if __name__ == '__main__':
    load_gene(*sys.argv[1:])







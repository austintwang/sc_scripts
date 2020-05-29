#!/usr/bin/env python3

import os
import sys
import pickle
import numpy as np
import pysam

def read_vcf(vcf_path, contig, start, end, min_maf=0., min_info=0.):
    genotypes_list = []
    markers = [] 
    marker_ids = []
    vcf = pysam.VariantFile(vcf_path)
    samples = vcf.header.samples
    for record in vcf.fetch(contig, start, end):
        maf = record.info["RefPanelAF"][0]
        info = record.info["INFO"]
        # print(maf) ####
        # print(info) ####
        if (maf < min_maf) or (info < min_info):
            continue

        add_marker = True
        record_gens = np.zeros((len(samples), 2,), dtype=int)
        for ind, sample in enumerate(samples):
            val = record.samples[sample]
            if not val.phased:
                add_marker = False
                break
            record_gens[ind,:] = np.array(val.allele_indices, dtype=int)
            # print(val.allele_indices) ####

        if add_marker:
            genotypes_list.append(record_gens)
            markers.append((contig, record.pos-1,),)
            marker_ids.append(record.id) 

    if len(genotypes_list) == 0:
        genotypes = np.array([])
    else:
        genotypes = np.stack(genotypes_list, axis=1)
    return genotypes, samples, markers, marker_ids

def add_data(agg_counts, var_data, cell_map, genotypes, sample_gen_map, marker_gen_map, well_only):
    # print(var_data) ####
    for var, cells in var_data.items():
        for cell, counts in cells.items():
            cell_full = cell
            if well_only:
                cell = cell[1]
            if not (cell in cell_map):
                # print(cell) ####
                continue
            if not (cell_map[cell] in sample_gen_map):
                # print(cell_map[cell]) ####
                continue
            counts_all = np.sum(counts)
            cell_agg = agg_counts.setdefault(cell_full, np.array([0,0,0]))
            cell_agg[2] += counts_all
            # print(var, var in marker_gen_map) ####
            if not (var in marker_gen_map):
                # print(var) ####
                continue
            sample_idx = sample_gen_map[cell_map[cell]]
            marker_idx = marker_gen_map[var]
            gen = genotypes[sample_idx, marker_idx, :]
            # print(gen) ####
            if np.sum(gen) == 1:
                cell_agg[:2] += counts[gen]

def process_samplename_kellis(sample_names):
    return [i[-8:] for i in sample_names]

def process_countsname_kellis(counts_names):
    return [i.split("_")[-1] for i in counts_names]

def load_gene(gene_name, dataset_name, radius, min_maf, min_info, well_only, ignore_total, data_dir, vcf_path, barcodes_map_path, boundaries_map_path, tss_map_path, total_counts_norm_path, status_path):
    with open(status_path, "w") as status_file:
        status_file.write("")

    gene_dir = os.path.join(data_dir, gene_name)
    if dataset_name == "Kellis":
        sample_process_fn = process_samplename_kellis
        counts_process_fn = process_countsname_kellis

    with open(barcodes_map_path, "rb") as barcodes_map_file:
        barcodes_map = pickle.load(barcodes_map_file)
    with open(boundaries_map_path, "rb") as boundaries_map_file:
        boundaries_map = pickle.load(boundaries_map_file)
    with open(tss_map_path, "rb") as tss_map_file:
        tss_map = pickle.load(tss_map_file)
    ignore_total = ignore_total == "True" 
    if ignore_total:
        total_counts_norm = None
    else:
        with open(total_counts_norm_path, "rb") as total_counts_norm_file:
            total_counts_norm = pickle.load(total_counts_norm_file)

    radius = int(radius)
    min_maf = float(min_maf)
    min_info = float(min_info)
    well_only = well_only == "True"
    # print(barcodes_map) ####

    if gene_name.split(".")[0] in tss_map:
        contig, start, end = boundaries_map[gene_name]
        genotypes, samples, markers, marker_ids = read_vcf(vcf_path, contig, start, end)
        samples = sample_process_fn(samples)
        sample_gen_map = dict([(val, ind) for ind, val in enumerate(samples)])
        marker_gen_map = dict([(val, ind) for ind, val in enumerate(markers)])
        # print(marker_gen_map.keys()) ####

        contig, tss_pos = tss_map[gene_name.split(".")[0]]
        contig = contig[3:]
        genotypes_nc, samples_nc, markers_nc, marker_ids_nc = read_vcf(
            vcf_path, contig, max(0, tss_pos - radius), tss_pos + radius + 1, min_maf=min_maf, min_info=min_info
        )
        samples_nc = sample_process_fn(samples_nc)
        # sample_gen_map_nc = dict([(val, ind) for ind, val in enumerate(samples_nc)])
        # marker_gen_map_nc = dict([(val, ind) for ind, val in enumerate(markers_nc)])

        total_counts_dir = os.path.join(gene_dir, "processed_counts")
        # print(ignore_total) ####
        if ignore_total or (not os.path.isdir(total_counts_dir)):
            # print(ignore_total or (not os.isdir(total_counts_dir))) ####
            total_counts = False
        else:
            # total_counts = {"_all": {}}
            total_counts = {}
            fnames = os.listdir(total_counts_dir)
            cell_types = counts_process_fn(fnames)
            for fname, cell_type in zip(fnames, cell_types):
                path = os.path.join(total_counts_dir, fname)
                with open(path, "rb") as count_file:
                    counts = pickle.load(count_file)
                    total_counts.setdefault(cell_type, {}).update(counts)
                    # for k, v in counts.items():
                    #     # print(k) ####
                    #     total_counts["_all"].setdefault(k, 0)
                    #     total_counts["_all"][k] += v
            # print(total_counts) ####


        agg_counts = {}
        var_data_paths = os.listdir(os.path.join(gene_dir, "bamdata")) 
        for path in var_data_paths:
            with open(os.path.join(gene_dir, "bamdata", path), "rb") as var_file:
                var_data = pickle.load(var_file)
                add_data(agg_counts, var_data, barcodes_map, genotypes, sample_gen_map, marker_gen_map, well_only)

        # print(agg_counts) ####
        # print(genotypes_nc) ####
        out_data = {
            "name": gene_name, 
            "genotypes": genotypes_nc, 
            "samples": samples_nc, 
            "markers": markers_nc, 
            "marker_ids": marker_ids_nc,
            "cell_counts": agg_counts,
            "total_counts": total_counts,
            "counts_norm": total_counts_norm
        }
        out_path = os.path.join(gene_dir, "gene_data.pickle")
        with open(out_path, "wb") as out_file:
            pickle.dump(out_data, out_file)

    with open(status_path, "w") as status_file:
        status_file.write("Complete")

if __name__ == '__main__':
    load_gene(*sys.argv[1:])







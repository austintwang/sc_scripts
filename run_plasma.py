#!/usr/bin/env python3

import numpy as np
import os
import sys
import traceback
import pickle
import gc
import shutil

if __name__ == '__main__' and __package__ is None:
    __package__ = 'run'
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    sys.path.insert(0, "/agusevlab/awang/plasma")
    
from . import Finemap, FmBenner

def restore_informative(shape, values, informative_snps, default):
    # print(shape) ####
    # print(values) ####
    # print(informative_snps) ####
    vals_all = np.full(shape, default)
    vals_all[informative_snps] = values
    return vals_all

def run_model(model_cls, inputs, input_updates, informative_snps, return_stats=False):
    model_inputs = inputs.copy()
    model_inputs.update(input_updates)
    # print(model_inputs) ####

    model = model_cls(**model_inputs)
    model.initialize()

    if inputs["search_mode"] == "exhaustive":
        model.search_exhaustive(inputs["min_causal"], inputs["max_causal"])
    elif inputs["search_mode"] == "shotgun":
        model.search_shotgun(
            inputs["min_causal"], 
            inputs["max_causal"], 
            inputs["prob_threshold"], 
            inputs["streak_threshold"], 
            inputs["search_iterations"]
        )

    shape_orig = np.shape(inputs["snp_ids"])

    causal_set_inf = model.get_causal_set(inputs["confidence"])
    causal_set = restore_informative(shape_orig, causal_set_inf, informative_snps, 1)
    ppas_inf = model.get_ppas()
    ppas = restore_informative(shape_orig, ppas_inf, informative_snps, np.nan)
    size_probs = model.get_size_probs()

    if return_stats:
        z_phi = restore_informative(shape_orig, model.imbalance_stats, informative_snps, np.nan)
        z_beta = restore_informative(shape_orig, model.total_exp_stats, informative_snps, np.nan)
        phi = restore_informative(shape_orig, model.phi, informative_snps, np.nan)
        beta = restore_informative(shape_orig, model.beta, informative_snps, np.nan)
        imbalance_errors = model.imbalance_errors 
        imbalance = model.imbalance

    gc.collect()

    # print(causal_set) ####
    if return_stats:
        return causal_set, ppas, size_probs, z_phi, z_beta, phi, beta, imbalance_errors, imbalance
    else:
        return causal_set, ppas, size_probs

def calc_reads(cell_counts, barcodes, barcodes_map, sample_names):
    sample_counts = {}
    sample_counts_norm = {}
    sample_num_cells = {}
    if isinstance(next(iter(barcodes_map.keys())), str):
        well_only = True
    else:
        well_only = False
    for i in barcodes:
        if well_only:
            bar_key = i[1]
        else:
            bar_key = i
        if (i in cell_counts) and (bar_key in barcodes_map):
            counts = cell_counts[i]
            sample = barcodes_map[bar_key]
            sc = sample_counts.setdefault(sample, np.array([0,0,0])) 
            sc += counts
            sample_num_cells.setdefault(sample, 0)
            sample_num_cells[sample] += 1

    counts_all = np.stack([sample_counts.get(i, np.array([0,0,0])) for i in sample_names], axis=0)
    num_cells_all = np.array([sample_num_cells.get(i, 0) for i in sample_names])
    return counts_all, num_cells_all

def load_clusters(gene_data, cluster_map_path, barcodes_map_path, overdispersion_path):
    # with open(gene_path, "rb") as gene_file:
    #     gene_data = pickle.load(gene_file)
    with open(cluster_map_path, "rb") as cluster_map_file:
        cluster_map = pickle.load(cluster_map_file)
    with open(barcodes_map_path, "rb") as barcodes_map_file:
        barcodes_map = pickle.load(barcodes_map_file)
    with open(overdispersion_path, "rb") as overdispersion_file:
        overdispersion = pickle.load(overdispersion_file)

    cell_counts = gene_data["cell_counts"]
    sample_names = gene_data["samples"]
    # sample_map = {val: ind for ind, val in enumerate(sample_names)}

    cluster_inputs = {}
    for cluster, barcodes in cluster_map.items():
        counts, num_cells = calc_reads(cell_counts, barcodes, barcodes_map, sample_names)
        overdispersion_clust = np.array([overdispersion[cluster].get(i, np.nan) for i in sample_names])
        cluster_inputs[cluster] = {
            "counts1": counts[:,0],
            "counts2": counts[:,1],
            "counts_total": counts[:,2],
            "overdispersion": overdispersion_clust,
            "num_cells": num_cells
        }

    return cluster_inputs

def write_output(output_path, result):
    # if not os.path.exists(output_path):
    #     os.makedirs(output_path)

    # output_return = os.path.join(output_path, "output.pickle")
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)
    with open(output_path, "wb") as output_file:
        pickle.dump(result, output_file)

    gc.collect()

def run_plasma(name, data_dir, params_path, filter_path, cluster_map_path, barcodes_map_path, overdispersion_path, status_path):
    with open(status_path, "w") as status_file:
        status_file.write("")

    try:
        gene_dir = os.path.join(data_dir, name)
        gene_path = os.path.join(gene_dir, "gene_data.pickle")

        with open(params_path, "rb") as params_file:
            params = pickle.load(params_file)

        output_dir = os.path.join(gene_dir, params["run_name"])
        os.makedirs(output_dir, exist_ok=True)
        output_path_base = os.path.join(gene_dir, params["run_name"], "plasma_{0}.pickle")

        with open(gene_path, "rb") as gene_file:
            gene_data = pickle.load(gene_file)

        if filter_path == "all":
            snp_filter = False
        else:
            with open(filter_path, "rb") as filter_file:
                snp_filter = pickle.load(filter_file)

        inputs_all = {
            "hap1": gene_data["genotypes"][:,:,0],
            "hap2": gene_data["genotypes"][:,:,1],
            "sample_names": gene_data["samples"],
            "snp_ids": gene_data["marker_ids"],
            "snp_pos": gene_data["markers"],
            "total_counts": gene_data.get("total_counts", False),
            "agg_counts": gene_data.get("counts_norm", False)
        }
        inputs_all.update(params)

        if inputs_all["total_counts"]:
            for k, v in inputs_all["cell_type_aliases"].items():
                inputs_all["total_counts"][v] = inputs_all["total_counts"][k]
        if inputs_all["agg_counts"]:
            for k, v in inputs_all["cell_type_aliases"].items():
                inputs_all["agg_counts"][v] = inputs_all["agg_counts"][k]

        clusters = load_clusters(gene_data, cluster_map_path, barcodes_map_path, overdispersion_path)

    except Exception as e:
        results = {}
        all_complete = False
        trace = traceback.format_exc()
        print(trace, file=sys.stderr)
        message = repr(e)
        results["run_error_all"] = message
        results["traceback_all"] = trace
        write_output(output_path_base.format(0), results)

        with open(status_path, "w") as status_file:
            status_file.write("Fail")
        return

    # print(clusters.keys()) ####
    # print(inputs_all["total_counts"].keys()) ####
    # print(inputs_all["total_counts"]) ####
    # print(inputs_all["agg_counts"])  ####

    splits = np.array(inputs_all.get("splits", [1.]))
    num_samples = len(inputs_all["sample_names"])
    allocs_raw = num_samples * splits
    cumu = np.cumsum(allocs_raw)
    rems = 1 - (-cumu % 1)
    # floors = (cumu - rems).astype(int)
    adds = np.random.binomial(1, rems)
    # cumu_int = floors + adds + (1 - np.roll(adds, 1))
    # allocs = np.copy(cumu_int)
    # allocs[1:] -= cumu_int[:-1]
    allocs_int = allocs_raw - rems
    allocs = (allocs_int + adds + (1 - np.roll(adds, 1))).astype(int)
    part_idxs = np.concatenate([np.full(val, ind) for ind, val in enumerate(allocs)])
    partitions = np.random.permutation(part_idxs)
    # print(rems) ####
    # print(adds) ####
    # print(allocs) ####
    # print(partitions) ####
    # print(np.sum(allocs), num_samples) ####

    all_complete = True
    for split in range(len(splits)):
        results = {}
        output_path = output_path_base.format(split)
        for cluster, inputs in clusters.items():
            result = results.setdefault(cluster, {})
            try:
                inputs = inputs.copy()
                inputs.update(inputs_all)
                print(cluster) ####
                # print(inputs["total_counts"].keys()) ####
                # print(inputs["total_counts"]) ####
                if inputs["total_counts"] and inputs["total_counts"].get(cluster, False):
                    processed_counts = True
                    # print(inputs["total_counts"][cluster]) ####
                    # print(inputs["total_counts"][cluster]) ####
                    # print(inputs["agg_counts"][cluster]) ####
                    inputs["counts_total"] = np.array([inputs["total_counts"][cluster].get(i, np.nan) for i in inputs["sample_names"]])
                    inputs["counts_norm"] = np.array([inputs["agg_counts"][cluster].get(i, np.nan) for i in inputs["sample_names"]])
                    # print(inputs["counts_total"]) ####
                    # print(inputs["counts_norm"]) ####
                    # print(np.count_nonzero(~np.isnan(inputs["counts_total"]))) ####
                    # print(np.count_nonzero(~np.isnan(inputs["counts_norm"]))) ####
                else:
                    processed_counts = False

                select_counts = np.logical_and.reduce([
                    partitions == split,
                    inputs["counts1"] >= 1, 
                    inputs["counts2"] >= 1, 
                    np.logical_not(np.isnan(inputs["overdispersion"]))
                ], axis=0)
                result["split"] = split
                result["effective_sample_size"] = np.sum(select_counts)
                result["sample_size"] = select_counts.size

                inputs["hap1"] = inputs["hap1"][select_counts]
                inputs["hap2"] = inputs["hap2"][select_counts]
                inputs["counts1"] = inputs["counts1"][select_counts]
                inputs["counts2"] = inputs["counts2"][select_counts]
                inputs["counts_total"] = inputs["counts_total"][select_counts]
                inputs["overdispersion"] = inputs["overdispersion"][select_counts]
                inputs["sample_names"] = np.array(inputs["sample_names"])[select_counts]
                inputs["num_cells"] = inputs["num_cells"][select_counts]
                if processed_counts:
                    inputs["counts_norm"] = inputs["counts_norm"][select_counts]

                result["avg_counts_total"] = np.nanmean(inputs["counts_total"])
                result["avg_counts_mapped"] = np.nanmean(inputs["counts1"] + inputs["counts2"])
                result["overdispersion"] = inputs["overdispersion"]
                # result["avg_overdispersion"] = np.nanmean(inputs["overdispersion"])
                result["avg_num_cells"] = np.nanmean(inputs["num_cells"])

                if processed_counts:
                    # print(inputs["counts_total"]) ####
                    # print(inputs["counts_norm"]) ####
                    inputs["counts_total"] = inputs["counts_total"] * np.mean(inputs["counts_norm"]) / inputs["counts_norm"]
                    result["avg_counts_total_scaled"] = np.nanmean(inputs["counts_total"])
                else:
                    result["avg_counts_total_scaled"] = None

                if snp_filter:
                    snps_in_filter = [ind for ind, val in enumerate(inputs["snp_ids"]) if val in snp_filter]
                    inputs["snp_ids"] = inputs["snp_ids"][snps_in_filter]
                    inputs["snp_pos"] = inputs["snp_pos"][snps_in_filter]
                    inputs["hap1"] = inputs["hap1"][:, snps_in_filter]
                    inputs["hap2"] = inputs["hap2"][:, snps_in_filter]

                haps_comb = inputs["hap1"] + inputs["hap2"]

                if np.size(inputs["counts1"]) <= 1:
                    result["data_error"] = "Insufficient Read Counts"
                    continue

                informative_snps = np.nonzero(np.logical_not(np.all(haps_comb == haps_comb[0,:], axis=0)))[0]
                result["informative_snps"] = informative_snps
                result["num_snps_total"] = np.size(inputs["snp_ids"])
                result["snp_ids"] = inputs["snp_ids"]
                result["num_snps_informative"] = np.count_nonzero(informative_snps)

                inputs["hap1"] = inputs["hap1"][:, informative_snps]
                inputs["hap2"] = inputs["hap2"][:, informative_snps]

                inputs["num_causal_prior"] = inputs["num_causal"]

                if inputs["hap1"].size == 0:
                    result["data_error"] = "Insufficient Markers"
                    continue

                inputs["hap_A"] = inputs["hap1"].astype(np.int)
                inputs["hap_B"] = inputs["hap2"].astype(np.int)

                inputs["counts_A"] = inputs["counts1"].astype(np.int)
                inputs["counts_B"] = inputs["counts2"].astype(np.int)
                inputs["total_exp"] = inputs["counts_total"].astype(np.int)

                results["hap_A"] = inputs["hap_A"]
                results["hap_B"] = inputs["hap_B"]
                results["counts_A"] = inputs["counts_A"]
                results["counts_B"] = inputs["counts_B"]
                results["total_exp"] = inputs["total_exp"]

                if inputs["model_flavors"] == "all":
                    model_flavors = set(["full", "indep", "eqtl", "ase", "acav", "fmb"])
                else:
                    model_flavors = inputs["model_flavors"]

                if "full" in model_flavors:
                    updates_full = {"num_ppl": None}
                    result["causal_set_full"], result["ppas_full"], result["size_probs_full"] = run_model(
                        Finemap, inputs, updates_full, informative_snps
                    )

                if "indep" in model_flavors:
                    updates_indep = {"cross_corr_prior": 0., "num_ppl": None}
                    result["causal_set_indep"], result["ppas_indep"], result["size_probs_indep"], result["z_phi"], result["z_beta"] , result["phi"], result["beta"], result["imbalance_errors"], result["imbalance"] = run_model(
                        Finemap, inputs, updates_indep, informative_snps, return_stats=True
                    )
                    
                if "eqtl" in model_flavors:
                    updates_eqtl = {"qtl_only": True, "num_ppl": None}
                    result["causal_set_eqtl"], result["ppas_eqtl"], result["size_probs_eqtl"] = run_model(
                        Finemap, inputs, updates_eqtl, informative_snps
                    )

                if "ase" in model_flavors:
                    updates_ase = {"as_only": True, "num_ppl": None}
                    result["causal_set_ase"], result["ppas_ase"], result["size_probs_ase"] = run_model(
                        Finemap, inputs, updates_ase, informative_snps
                    )

                if "fmb" in model_flavors:
                    updates_fmb = {"qtl_only": True, "num_ppl": None}
                    result["causal_set_fmb"], result["ppas_fmb"], result["size_probs_fmb"] = run_model(
                        FmBenner, inputs, updates_fmb, informative_snps
                    )

            except Exception as e:
                all_complete = False
                trace = traceback.format_exc()
                print(trace, file=sys.stderr)
                message = repr(e)
                result["run_error"] = message
                result["traceback"] = trace
                write_output(output_path, results)

        write_output(output_path, results)

    with open(status_path, "w") as status_file:
        if all_complete:
            status_file.write("Complete")
        else:
            status_file.write("Fail")

if __name__ == '__main__':
    run_plasma(*sys.argv[1:])


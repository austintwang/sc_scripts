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
    # print(shape_orig, len(causal_set_inf), informative_snps.shape) ####
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

def calc_reads(cell_counts, cell_moments, barcodes, barcodes_map, sample_names):
    sample_counts = {}
    sample_num_cells = {}
    sample_counts_norm = {}
    sample_moments = {}
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
            moments = cell_moments.get(i, np.zeros((4,2),))
            sample = barcodes_map[bar_key]
            sc = sample_counts.setdefault(sample, np.array([0,0,0])) 
            sm = sample_moments.setdefault(sample, np.zeros((4,2),))
            sc += counts
            sm += moments
            sample_num_cells.setdefault(sample, 0)
            sample_num_cells[sample] += 1

    counts_all = np.stack([sample_counts.get(i, np.array([0,0,0])) for i in sample_names], axis=0)
    moments_all = np.stack([sample_moments.get(i, np.zeros((4,2),)) for i in sample_names], axis=0)
    num_cells_all = np.array([sample_num_cells.get(i, 0) for i in sample_names])

    r1 = moments_all[:,1,:] / moments_all[:,0,:]
    r2 = moments_all[:,2,:] / moments_all[:,1,:]
    r3 = moments_all[:,3,:] / moments_all[:,2,:]
    r1r3 = r1 * r3
    r2r3 = r2 * r3
    r1r2 = r1 * r2
    k_on = 2*r1*(r3-r2)/(r1r2-2*r1r3+r2r3)
    k_off = 2*(r2-r1)*(r1-r3)*(r3-r2)/((r1r2-2*r1r3+r2r3)*(r1-2*r2+r3))
    k_syn = (-r1r2+2*r1r3-r2r3)/(r1-2*r2+r3)

    return counts_all, num_cells_all, moments_all, k_on, k_off, k_syn

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
    cell_moments = gene_data["burst_data"]
    sample_names = gene_data["samples"]
    # sample_map = {val: ind for ind, val in enumerate(sample_names)}

    cluster_inputs = {}
    for cluster, barcodes in cluster_map.items():
        # print(cluster_map.keys()) ####
        counts, num_cells, moments, k_on, k_off, k_syn = calc_reads(cell_counts, cell_moments, barcodes, barcodes_map, sample_names)
        overdispersion_clust = np.array([overdispersion[cluster].get(i, np.nan) for i in sample_names])
        cluster_inputs[cluster] = {
            "counts1": counts[:,0],
            "counts2": counts[:,1],
            "counts_total": counts[:,2],
            "overdispersion": overdispersion_clust,
            "num_cells": num_cells,
            "moments1": moments[:,:,0],
            "moments2": moments[:,:,1],
            "kon1": k_on[:,0],
            "kon2": k_on[:,1],
            "koff1": k_off[:,0],
            "koff2": k_off[:,1],
            "ksyn1": k_syn[:,0],
            "ksyn2": k_syn[:,1],
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
            "snp_ids": np.array(gene_data["marker_ids"]),
            "snp_pos": np.array(gene_data["markers"]),
            "snp_alleles": np.array(gene_data["marker_alleles"]),
            "total_counts": gene_data.get("total_counts", False),
            "agg_counts": gene_data.get("counts_norm", False),
            "tss": gene_data.get("tss"),
            "sample_masks": gene_data.get("sample_masks", {})
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
    adds = np.random.binomial(1, rems)
    allocs_int = allocs_raw - rems
    allocs = (allocs_int + adds + (1 - np.roll(adds, 1))).astype(int)
    part_idxs = np.concatenate([np.full(val, ind) for ind, val in enumerate(allocs)])
    partitions = np.random.permutation(part_idxs)
    # print(partitions) ####

    all_complete = True
    for all_but in [False, True]:
        for split in range(len(splits)):
            select_counts = np.logical_xor(partitions == split, all_but)
            results = {"_gen": {}}
            results["_gen"]["hap_A"] = inputs_all["hap1"][select_counts]
            results["_gen"]["hap_B"] = inputs_all["hap2"][select_counts]
            results["_gen"]["snp_ids"] = inputs_all["snp_ids"]
            results["_gen"]["snp_alleles"] = inputs_all["snp_alleles"]
            results["_gen"]["snp_pos"] = inputs_all["snp_pos"]
            results["_gen"]["tss"] = inputs_all["tss"]
            out_prefix = "x" if all_but else "i"
            output_path = output_path_base.format(out_prefix + str(split))
            for cluster, inputs in clusters.items():
                result = results.setdefault(cluster, {})
                try:
                    inputs = inputs.copy()
                    inputs.update(inputs_all)
                    print(cluster, split, all_but) ####
                    # print(inputs["total_counts"].keys()) ####
                    # print(inputs["total_counts"].keys()) ####
                    if inputs["total_counts"] and inputs["total_counts"].get(cluster, False):
                        processed_counts = True
                        # print(inputs["total_counts"][cluster]) ####
                        # print(inputs["total_counts"][cluster]) ####
                        # print(inputs["agg_counts"][cluster]) ####
                        # print(inputs["total_counts"][cluster].keys()) ####
                        inputs["counts_total"] = np.array([inputs["total_counts"][cluster][inputs["pre_flags"]].get(i, np.nan) for i in inputs["sample_names"]])
                        # inputs["counts_norm"] = np.array([inputs["agg_counts"][cluster].get(i, np.nan) for i in inputs["sample_names"]])
                        # counts_raw_total = np.array([inputs["total_counts"][f"{cluster}_raw"].get(i, np.nan) for i in inputs["sample_names"]])
                        # counts_raw_norm = np.array([inputs["agg_counts"][f"{cluster}_raw"].get(i, np.nan) for i in inputs["sample_names"]])
                        # results["counts_raw"] = np.array([inputs["total_counts"][cluster]["cm"].get(i, np.nan) for i in inputs["sample_names"]])
                        # print(inputs["counts_total"]) ####
                        # print(inputs["counts_norm"]) ####
                        # print(np.count_nonzero(~np.isnan(inputs["counts_total"]))) ####
                        # print(np.count_nonzero(~np.isnan(inputs["counts_norm"]))) ####
                    else:
                        processed_counts = False

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

                    result["moments_A"] =result["moments1"] = inputs["moments1"][select_counts]
                    result["moments_B"] = result["moments2"] = inputs["moments2"][select_counts]
                    result["kon_A"] = result["kon1"] = inputs["kon1"][select_counts]
                    result["kon_B"] = result["kon2"] = inputs["kon2"][select_counts]
                    result["koff_A"] = result["koff1"] = inputs["koff1"][select_counts]
                    result["koff_B"] = result["koff2"] = inputs["koff2"][select_counts]
                    result["ksyn_A"] = result["ksyn1"] = inputs["ksyn1"][select_counts]
                    result["ksyn_B"] = result["ksyn2"] = inputs["ksyn2"][select_counts]
                    # if processed_counts:
                    #     inputs["counts_norm"] = inputs["counts_norm"][select_counts]

                    # print(inputs.get("clinical_group", True)) ####
                    clinical_mask = inputs.get("sample_masks", {}).get(inputs.get("clinical_group"), True)
                    print(clinical_mask) ####
                    # inputs["mask_imbalance"] = mask_imbalance = np.logical_and.reduce([
                    #     clinical_mask,
                    #     inputs["counts1"] >= 1, 
                    #     inputs["counts2"] >= 1, 
                    #     np.logical_not(np.isnan(inputs["overdispersion"]))
                    # ], axis=0)

                    inputs["mask_imbalance"] = mask_imbalance = (
                        clinical_mask
                        & (inputs["counts1"] >= 1)
                        & (inputs["counts2"] >= 1)
                        & ~np.isnan(inputs["overdispersion"])
                    )
                    # print(inputs.get("clinical_group", True)) ####
                    # print(np.logical_not(np.isnan(inputs["counts_total"]))) ####
                    inputs["mask_total_exp"] = mask_total_exp = clinical_mask & ~np.isnan(inputs["counts_total"])
                    print(np.sum(mask_total_exp)) ####
                    # print(inputs["counts_total"][mask_total_exp]) ####

                    result["avg_counts_total"] = np.nanmean(inputs["counts_total"])
                    result["avg_counts_mapped"] = np.nanmean(inputs["counts1"] + inputs["counts2"])
                    result["overdispersion"] = inputs["overdispersion"]
                    # result["avg_overdispersion"] = np.nanmean(inputs["overdispersion"])
                    result["avg_num_cells"] = np.nanmean(inputs["num_cells"])

                    if processed_counts:
                        # print(inputs["counts_total"]) ####
                        # print(inputs["counts_norm"]) ####
                        # print(processed_counts) ####
                        # print(inputs["counts_total"][inputs["mask_total_exp"]]) ####
                        # print(np.mean(inputs["counts_norm"]) / inputs["counts_norm"]) ####
                        # inputs["counts_total"] = inputs["counts_total"] * 1e3 / inputs["counts_norm"]
                        result["avg_counts_total_scaled"] = np.nanmean(inputs["counts_total"])
                        # print(inputs["counts_total"]) ####
                        # print(result["avg_counts_total_scaled"]) ####

                    else:
                        result["avg_counts_total_scaled"] = None

                    if snp_filter:
                        snps_in_filter = [ind for ind, val in enumerate(inputs["snp_ids"]) if val in snp_filter]
                        inputs["snp_ids"] = inputs["snp_ids"][snps_in_filter]
                        inputs["snp_pos"] = inputs["snp_pos"][snps_in_filter]
                        inputs["hap1"] = inputs["hap1"][:, snps_in_filter]
                        inputs["hap2"] = inputs["hap2"][:, snps_in_filter]

                    # haps_comb = (inputs["hap1"] + inputs["hap2"])[mask_total_exp, :]
                    # haps_diff = (inputs["hap1"] - inputs["hap2"])[mask_imbalance, :]
                    hap_c1 = inputs["hap1"][mask_total_exp, :]
                    hap_c2 = inputs["hap2"][mask_total_exp, :]
                    hap_d1 = inputs["hap1"][mask_imbalance, :]
                    hap_d2 = inputs["hap2"][mask_imbalance, :]

                    if np.size(inputs["counts1"][mask_imbalance]) <= 1:
                        result["data_error"] = "Insufficient Read Counts"
                        print("data_error") ####
                        continue

                    mc = inputs["min_carriers"]
                    informative_snps = np.nonzero(np.logical_and.reduce([
                        # np.logical_not(np.all(haps_comb == haps_comb[0,:], axis=0)),
                        # np.logical_not(np.all(haps_diff == haps_diff[0,:], axis=0)),
                        np.sum(hap_c1, axis=0) >= mc,
                        np.sum(1 - hap_c1, axis=0) >= mc,
                        np.sum(hap_c2, axis=0) >= mc,
                        np.sum(1 - hap_c2, axis=0) >= mc,
                        np.sum(hap_d1, axis=0) >= mc,
                        np.sum(1 - hap_d1, axis=0) >= mc,
                        np.sum(hap_d2, axis=0) >= mc,
                        np.sum(1 - hap_d2, axis=0) >= mc,
                    ]))[0]
                    informative_snps_weak = np.nonzero(np.logical_and.reduce([
                        np.sum(hap_c1, axis=0) >= mc,
                        np.sum(1 - hap_c1, axis=0) >= mc,
                        np.sum(hap_c2, axis=0) >= mc,
                        np.sum(1 - hap_c2, axis=0) >= mc,
                    ]))[0]
                    result["informative_snps"] = informative_snps
                    result["informative_snps_weak"] = informative_snps_weak
                    result["num_snps_total"] = np.size(inputs["snp_ids"])
                    result["snp_ids"] = inputs["snp_ids"]
                    result["num_snps_informative"] = np.count_nonzero(informative_snps)
                    result["num_snps_informative_weak"] = np.count_nonzero(informative_snps_weak)

                    inputs["hap1"] = inputs["hap1"][:, informative_snps]
                    inputs["hap2"] = inputs["hap2"][:, informative_snps]

                    inputs["hap_A"] = inputs["hap1"].astype(np.int)
                    inputs["hap_B"] = inputs["hap2"].astype(np.int)

                    inputs["num_causal_prior"] = inputs["num_causal"]

                    if inputs["hap1"].size == 0:
                        result["data_error"] = "Insufficient Markers"
                        continue

                    inputs["counts_A"] = inputs["counts1"].astype(np.int)
                    inputs["counts_B"] = inputs["counts2"].astype(np.int)
                    inputs["total_exp"] = inputs["counts_total"].astype(float)
                    # print(list(inputs["total_exp"])) ####
                    # print(inputs["counts_norm"][inputs["mask_total_exp"]]) ####
                    # print(inputs["counts_total"][inputs["mask_total_exp"]]) ####
                    # print(inputs["total_exp"][inputs["mask_total_exp"]]) ####

                    result["counts_A"] = inputs["counts_A"]
                    result["counts_B"] = inputs["counts_B"]
                    result["total_exp"] = inputs["total_exp"]

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


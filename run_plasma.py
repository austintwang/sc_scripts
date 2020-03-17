#!/usr/bin/env python3

import numpy as np
import os
import sys
import traceback
import pickle
import gc

if __name__ == '__main__' and __package__ is None:
    __package__ = 'run'
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    sys.path.insert(0, "/agusevlab/awang/plasma")
    
from . import Finemap, FmBenner

def restore_informative(shape, values, informative_snps, default):
    vals_all = np.full(shape, default)
    np.put(vals_all, informative_snps, values)
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

    shape_orig = np.ones(np.shape(inputs["snp_ids"]))

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

    gc.collect()

    # print(causal_set) ####
    if return_stats:
        return causal_set, ppas, size_probs, z_phi, z_beta, phi, beta
    else:
        return causal_set, ppas, size_probs

def calc_reads(cell_counts, barcodes, barcodes_map, sample_names):
    sample_counts = {}
    for i in barcodes:
        if (i in cell_counts) and (i in barcodes_map):
            counts = cell_counts[i]
            sample = barcodes_map[i]
            sc = sample_counts.setdefault(sample, np.array([0,0,0])) 
            sc += counts

    counts_all = np.stack([sample_counts.get(i, np.array([0,0,0])) for i in sample_names], axis=0)
    return counts_all

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
        counts = calc_reads(cell_counts, barcodes, barcodes_map, sample_names)
        overdispersion_clust = np.array([overdispersion[cluster].get(i, np.nan) for i in sample_names])
        cluster_inputs[cluster] = {
            "counts1": counts[:,0],
            "counts2": counts[:,1],
            "counts_total": counts[:,2],
            "overdispersion": overdispersion_clust
        }

    return cluster_inputs

def write_output(output_path, result):
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    output_return = os.path.join(output_path, "output.pickle")
    with open(output_return, "wb") as output_file:
        pickle.dump(result, output_file)

    gc.collect()

def run_plasma(name, data_dir, params_path, filter_path, cluster_map_path, barcodes_map_path, overdispersion_path, status_path):
    with open(status_path, "w") as status_file:
        status_file.write("")

    try:
        gene_dir = os.path.join(data_dir, name)
        gene_path = os.path.join(gene_dir, "gene_data.pickle")
        output_path = os.path.join(gene_dir, "plasma.pickle")

        with open(params_path, "rb") as params_file:
            params = pickle.load(params_file)
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
            "snp_pos": gene_data["markers"]
        }
        inputs_all.update(params)

        clusters = load_clusters(gene_data, cluster_map_path, barcodes_map_path, overdispersion_path)

    except Exception as e:
        all_complete = False
        trace = traceback.format_exc()
        print(trace, file=sys.stderr)
        message = repr(e)
        result["run_error"] = message
        result["traceback"] = trace
        write_output(output_path, results)

        with open(status_path, "w") as status_file:
            status_file.write("Fail")
        return

    results = {}
    all_complete = True
    for cluster, inputs in clusters.items():
        result = results.setdefault(cluster, {})
        try:
            inputs.update(inputs_all)

            select_counts = np.logical_and(
                inputs["counts1"] >= 1, 
                inputs["counts2"] >= 1, 
                np.logical_not(np.isnan(inputs["overdispersion"]))
            ) 
            inputs["hap1"] = inputs["hap1"][select_counts]
            inputs["hap2"] = inputs["hap2"][select_counts]
            inputs["counts1"] = inputs["counts1"][select_counts]
            inputs["counts2"] = inputs["counts2"][select_counts]
            inputs["counts_total"] = inputs["counts_total"][select_counts]
            inputs["sample_names"] = inputs["sample_names"][select_counts]

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
                result["causal_set_indep"], result["ppas_indep"], result["size_probs_indep"], result["z_phi"], result["z_beta"] , result["phi"], result["beta"] = run_model(
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


#!/usr/bin/env python3

import numpy as np
import os
import sys
import traceback
import pickle
import gc
import subprocess
import glob

if __name__ == '__main__' and __package__ is None:
    __package__ = 'run'
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    sys.path.insert(0, "/agusevlab/awang/plasma")
    
from . import Finemap

def run_plink_ld(gwas_gen_path, marker_ids, num_snps, contig):
    in_path = os.path.join("/tmp", str(np.random.randint(100000000)))
    out_path_base = os.path.join("/tmp", str(np.random.randint(100000000)))
    out_path = out_path_base + ".ld"
    cmd = [
        "/agusevlab/awang/plink/plink", 
        "--bfile", gwas_gen_path + contig, 
        "--r",
        "--extract", in_path, 
        "--out", out_path_base
    ]
    # print(" ".join(cmd)) ####
    
    with open(in_path, "w") as in_file:
        in_file.writelines([i + "\n" for i in marker_ids])
    out = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print(out) ####
    os.remove(in_path)

    ld = np.ones((num_snps, num_snps),)
    marker_map = dict([(val, ind) for ind, val in enumerate(marker_ids)])

    with open(out_path, "r") as out_file:
        next(out_file)
        for line in out_file:
            # print(line) ####
            # if line == "\n":
            #     continue
            data = line.strip().split()
            # print(data) ####
            id1 = data[2]
            id2 = data[5]
            corr = float(data[6])
            idx1 = marker_map[id1]
            idx2 = marker_map[id2]
            ld[idx1, idx2] = corr
            ld[idx2, idx1] = corr

    for path in glob.glob(out_path_base):
        os.remove(path)

    return ld

def restore_informative(shape, values, informative_snps, default):
    vals_all = np.full(shape, default)
    # print(vals_all) ####
    # print(informative_snps) ####
    # print(vals_all[informative_snps]) ####
    # print(values) ####
    # np.put(vals_all, informative_snps, values)
    vals_all[informative_snps] = values
    # print(vals_all) ####
    return vals_all

def run_model(model_cls, inputs, input_updates, informative_snps):
    model_inputs = inputs.copy()
    model_inputs.update(input_updates)

    model = model_cls(**model_inputs)
    model.initialize()
    # print(model.total_exp_var_prior)####

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

    shape_orig = inputs["num_snps_orig"]
    # print(shape_orig) ####
    # print(informative_snps.shape) ####

    causal_set_inf = model.get_causal_set(inputs["confidence"])
    # print(causal_set_inf) ####
    causal_set = restore_informative(shape_orig, causal_set_inf, informative_snps, 1)
    ppas_inf = model.get_ppas()
    # print(ppas_inf) ####
    ppas = restore_informative(shape_orig, ppas_inf, informative_snps, np.nan)
    # print(ppas) ####
    size_probs = model.get_size_probs()    
        
    return causal_set, ppas, size_probs

def write_output(output_path, result):
    with open(output_path, "wb") as output_file:
        pickle.dump(result, output_file)

    gc.collect()

def colocalize(gene_name, data_dir, params_path, filter_path, gwas_dir, gwas_gen_path, boundaries_map_path, status_path):
    with open(status_path, "w") as status_file:
        status_file.write("")

    gene_dir = os.path.join(data_dir, gene_name)
    gene_path = os.path.join(gene_dir, "gene_data.pickle")
    finemap_path = os.path.join(gene_dir, "combined", "plasma_i0.pickle")
    os.makedirs(os.path.join(gene_dir, "coloc"), exist_ok=True)
    
    all_complete = True

    with open(params_path, "rb") as params_file:
        params = pickle.load(params_file)
    
    try:
        with open(gene_path, "rb") as gene_file:
            gene_data = pickle.load(gene_file)
        with open(finemap_path, "rb") as finemap_file:
            finemap_data = pickle.load(finemap_file)
    except Exception as e:
        trace = traceback.format_exc()
        print(trace, file=sys.stderr)
        status_file.write("Complete")

    with open(boundaries_map_path, "rb") as boundaries_map_file:
        boundaries_map = pickle.load(boundaries_map_file)
    if filter_path == "all":
        snp_filter = False
    else:
        with open(filter_path, "rb") as filter_file:
            snp_filter = pickle.load(filter_file)

    contig, start, end = boundaries_map[gene_name]

    studies = os.listdir(gwas_dir)

    for study in studies:
        if study == "gen":
            continue
        # print(study) ####
        try:
            gwas_path = os.path.join(gwas_dir, study)
            gwas_name = study.split(".")[0]
            output_path = os.path.join(gene_dir, "coloc", "{0}.pickle".format(gwas_name))

            with open(gwas_path, "rb") as gwas_file:
                gwas_data = pickle.load(gwas_file)

            inputs = {
                "snp_ids": gene_data["marker_ids"],
                "snp_pos": gene_data["markers"],
                "z_beta": np.array([gwas_data.get(i, np.nan) for i in gene_data["marker_ids"]]),
                "num_snps_orig": len(gene_data["marker_ids"])
            }
            # print(inputs) ####
            inputs.update(params)

            inputs["num_ppl_total_exp"] = gwas_data["_size"]
        
            result = {"z_beta": inputs["z_beta"].copy()}
            informative_snps = np.logical_not(np.isnan(inputs["z_beta"]))
            result["informative_snps"] = informative_snps
            inputs["total_exp_stats"] = inputs["z_beta"][informative_snps]
            # print(inputs["total_exp_stats"]) ####
            inputs["snp_ids"] = np.array(inputs["snp_ids"])[informative_snps]
            inputs["num_snps"] = inputs["total_exp_stats"].size

            inputs["num_causal_prior"] = inputs["num_causal"]

            if inputs["num_snps"] == 0:
                result["data_error"] = "Insufficient Markers"
                write_output(output_path, result)
                return

            inputs["corr_shared"] = run_plink_ld(gwas_gen_path, inputs["snp_ids"], inputs["num_snps"], contig)
            # print(inputs["corr_shared"]) ####

            if inputs["model_flavors_gwas"] == "all":
                model_flavors_gwas = set(["eqtl"])
            else:
                model_flavors_gwas = inputs["model_flavors_gwas"]

            if inputs["model_flavors_qtl"] == "all":
                model_flavors_qtl = set(["full", "indep", "eqtl", "ase", "fmb"])
            else:
                model_flavors_qtl = inputs["model_flavors_qtl"]

            if "eqtl" in model_flavors_gwas:
                updates_eqtl = {"qtl_only": True}
                result["causal_set_eqtl"], result["ppas_eqtl"], result["size_probs_eqtl"] = run_model(
                    Finemap, inputs, updates_eqtl, informative_snps
                )

            # if "fmb" in model_flavors_gwas:
            #     updates_fmb = {"qtl_only": True}
            #     result["causal_set_fmb"], result["ppas_fmb"], result["size_probs_fmb"] = run_model(
            #         FmBenner, inputs, updates_fmb, informative_snps
            #     )
            # print(result["causal_set_eqtl"]) ####

            cluster_results = result.setdefault("clusters", {})
            for cluster, fm_res in finemap_data.items():
                cluster_results.setdefault(cluster, {})
                # print(fm_res.keys()) ####
                for fg in model_flavors_gwas:
                    for fq in model_flavors_qtl:
                        try:
                            # snps_used = result["ppas_{0}".format(fg)] != np.nan
                            snps_used = np.logical_not(np.isnan(result["ppas_{0}".format(fg)]))
                            # scale = np.nansum(fm_res["ppas_{0}".format(fq)][snps_used])
                            scale = 1.
                            # print(scale) ####
                            # print(snps_used) ####
                            # print(result["ppas_{0}".format(fg)]) ####
                            fm_res_scaled = fm_res["ppas_{0}".format(fq)] / scale
                            clpps = fm_res_scaled * result["ppas_{0}".format(fg)]
                            # print(fm_res["ppas_{0}".format(fq)]) ####
                            h4 = np.nansum(clpps)
                            cluster_results[cluster]["clpp_{0}_{1}".format(fq, fg)] = clpps
                            # if study == "BDSCZ_Ruderfer2018.pickle" and cluster == "Ex":
                                # print(cluster, fg, fq) ####
                                # print(fm_res_scaled) ####
                                # print(list(zip(gene_data["marker_ids"], fm_res_scaled, inputs["z_beta"], clpps))) ####
                            cluster_results[cluster]["h4_{0}_{1}".format(fq, fg)] = h4
                        except KeyError:
                            continue

            write_output(output_path, result)

        except Exception as e:
            all_complete = False
            trace = traceback.format_exc()
            print(trace, file=sys.stderr)
            message = repr(e)
            result = {}
            result["run_error"] = message
            result["traceback"] = trace
            write_output(output_path, result)
            return

    with open(status_path, "w") as status_file:
        if all_complete:
            status_file.write("Complete")
        else:
            status_file.write("Fail")

if __name__ == '__main__':
    colocalize(*sys.argv[1:])


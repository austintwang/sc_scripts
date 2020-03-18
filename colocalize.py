#!/usr/bin/env python3

import numpy as np
import os
import sys
import traceback
import pickle
import gc
import subprocess

if __name__ == '__main__' and __package__ is None:
    __package__ = 'run'
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
    sys.path.insert(0, "/agusevlab/awang/plasma")
    
from . import Finemap, FmBenner

def run_plink_ld(gwas_gen_path, marker_ids, contig):
    in_pipe_path = os.path.join("/tmp", str(np.randint(100000000)))
    out_path_base = os.path.join("/scratch", str(np.randint(100000000)))
    out_path = out_path + ".ld"
    cmd = [
        "/agusevlab/awang/plink/plink", 
        "--bfile", gwas_gen_path + "." + contig, 
        "-r",
        "--keep", pipe_path, 
        "--out", out_path_base
    ]
    
    os.mkfifo(in_pipe_path)
    with open(in_pipe_path, "w") as in_pipe:
        in_pipe.writelines([i + "\n" for i in marker_ids])
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    os.remove(in_pipe_path)

    with open(out_path, "r") as out_file:
        for line in out_file:
            print(line) ####
    os.remove(out_path)


def run_model(model_cls, inputs, input_updates, informative_snps):
    model_inputs = inputs.copy()
    model_inputs.update(input_updates)

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
        
    return causal_set, ppas, size_probs

def write_output(output_path, result):
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    output_return = os.path.join(output_path, "output.pickle")
    with open(output_return, "wb") as output_file:
        pickle.dump(result, output_file)

    gc.collect()

def colocalize(gwas_name, gene_name, data_dir, params_path, filter_path, gwas_path, gwas_gen_path, status_path):
    with open(status_path, "w") as status_file:
        status_file.write("")

    gene_dir = os.path.join(data_dir, gene_name)
    gene_path = os.path.join(gene_dir, "gene_data.pickle")
    finemap_path = os.path.join(gene_dir, "plasma.pickle")
    if os.path.isdir(finemap_path):
        finemap_path = os.path.join(finemap_path, "output.pickle")
    output_path = os.path.join(gene_dir, "coloc_{0}.pickle".format(gwas_name))

    try:
        with open(params_path, "rb") as params_file:
            params = pickle.load(params_file)
        with open(gwas_path, "rb") as gwas_file:
            gwas_data = pickle.load(gwas_file)
        with open(gene_path, "rb") as gene_file:
            gene_data = pickle.load(gene_file)
        with open(finemap_path, "rb") as finemap_file:
            finemap_data = pickle.load(finemap_file)
        if filter_path == "all":
            snp_filter = False
        else:
            with open(filter_path, "rb") as filter_file:
                snp_filter = pickle.load(filter_file)

        inputs = {
            "snp_ids": gene_data["marker_ids"],
            "snp_pos": gene_data["markers"],
            "z_beta": np.array([gwas_data.get(i, np.nan) for i in gene_data["marker_ids"]])
        }
        inputs.update(params)
    
        result = {"z_beta": inputs["z_beta"].copy()}
        informative_snps = np.logical_not(np.isnan(inputs["z_beta"]))
        result["informative_snps"] = informative_snps
        inputs["total_exp_stats"] = inputs["z_beta"][informative_snps]
        inputs["snp_ids"] = np.array(inputs["snp_ids"])[informative_snps]
        inputs["num_snps"] = inputs["total_exp_stats"].size

        inputs["num_causal_prior"] = inputs["num_causal"]

        if inputs["num_snps"] == 0:
            result["data_error"] = "Insufficient Markers"
            write_output(output_path, result)
            return

        inputs["corr_shared"] = run_plink_ld(gwas_gen_path, inputs["snp_ids"])

        if inputs["model_flavors_gwas"] == "all":
            model_flavors_gwas = set(["eqtl", "fmb"])
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

        if "fmb" in model_flavors_gwas:
            updates_fmb = {"qtl_only": True}
            result["causal_set_fmb"], result["ppas_fmb"], result["size_probs_fmb"] = run_model(
                FmBenner, inputs, updates_fmb, informative_snps
            )

        cluster_results = results.setdefault("clusters", {})
        for cluster, fm_res in finemap_data.items():
            for fg in model_flavors_gwas:
                for fq in model_flavors_qtl:
                    clpps = fm_res["ppas_{0}".format(fq)] * result["ppas_{0}".format(fg)]
                    h4 = np.nansum(clpps)
                    cluster_results[cluster]["clpp_{0}_{1}".format(fq, fg)] = clpps
                    cluster_results[cluster]["h4_{0}_{1}".format(fq, fg)] = h4

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


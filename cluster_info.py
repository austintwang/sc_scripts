import numpy as np
import scipy.stats
import os
import pickle
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['agg.path.chunksize'] = 10000
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import upsetplot
import tracemalloc ####

def read_data(plasma_data, clusters, gene_name, top_snps=None):
    # print(coloc_data) ####
    # print(gene_name) ####
    data = []
    for c in clusters:
        plasma_clust = plasma_data.get(c, None)
        if plasma_clust is None:
            # print(plasma_data) ####
            continue
        # if top_snps:
        #     print(plasma_clust) ####
        if "causal_set_indep" in plasma_clust:
            # print(plasma_clust["ppas_indep"]) ####
            ppa = True
            if top_snps is None:
                try:
                    top_snp = np.nanargmax(plasma_clust["ppas_indep"])
                    plasma_clust["snp_ids"][top_snp]
                except (ValueError, KeyError):
                    ppa = False
            else:
                # print(gene_name) ####
                try:
                    # top_snp = plasma_clust["snp_ids"].index(top_snps[c])
                    top_snp = np.argwhere(plasma_clust["snp_ids"] == top_snps[c])[0][0]
                except (ValueError, KeyError, IndexError):
                    # print("False") ####
                    ppa = False

            num_snps_total = plasma_clust["num_snps_total"]
            num_snps_informative = plasma_clust.get("num_snps_informative", num_snps_total)
            
            z_phi = plasma_clust["z_phi"][top_snp] if ppa else np.nan
            phi = plasma_clust["phi"][top_snp] if ppa else np.nan
            nlp_phi = -scipy.stats.norm.logsf(np.abs(z_phi)) / np.log(10) - np.log10(2) - np.log10(num_snps_informative)
            
            z_beta = plasma_clust["z_beta"][top_snp] if ppa else np.nan
            beta = plasma_clust["beta"][top_snp] if ppa else np.nan
            nlp_beta = -scipy.stats.norm.logsf(np.abs(z_beta)) / np.log(10) - np.log10(2) - np.log10(num_snps_informative)

            z_comb = (z_phi + z_beta) / np.sqrt(2)
            nlp_comb = -scipy.stats.norm.logsf(np.abs(z_comb)) / np.log(10) - np.log10(2) - np.log10(num_snps_informative)

            data_clust = [
                gene_name, 
                c, 
                plasma_clust["avg_counts_total"],
                plasma_clust["avg_counts_total_scaled"],
                plasma_clust["avg_counts_mapped"],
                np.mean(plasma_clust["overdispersion"]),
                plasma_clust["avg_num_cells"],
                plasma_clust["effective_sample_size"],
                plasma_clust["sample_size"],
                num_snps_informative,
                num_snps_total,
                np.sum(plasma_clust.get("causal_set_indep", np.nan)) - (num_snps_total - num_snps_informative), 
                np.sum(plasma_clust.get("causal_set_ase", np.nan)) - (num_snps_total - num_snps_informative),
                np.sum(plasma_clust.get("causal_set_eqtl", np.nan)) - (num_snps_total - num_snps_informative),
                (np.sum(plasma_clust.get("causal_set_indep", np.nan)) - num_snps_total) / num_snps_informative + 1, 
                (np.sum(plasma_clust.get("causal_set_ase", np.nan)) - num_snps_total) / num_snps_informative + 1,
                (np.sum(plasma_clust.get("causal_set_eqtl", np.nan)) - num_snps_total) / num_snps_informative + 1,
                plasma_clust["ppas_indep"][top_snp] if ppa else np.nan,
                z_comb,
                nlp_comb,
                z_phi,
                phi,
                nlp_phi,
                # np.count_nonzero(-np.log10(scipy.stats.norm.sf(abs(plasma_clust["z_phi"]))*2/plasma_clust["num_snps_informative"]) >= 1.301) if ppa else 0,
                z_beta,
                beta,
                nlp_beta,                
                plasma_clust["snp_ids"][top_snp] if ppa else None,
                plasma_clust.get("split", np.nan),
            ]
            # print(plasma_clust["snp_ids"][top_snp] if ppa else None) ####
            # print(data_clust) ####
            data.append(data_clust)
    return data

def make_df(run_name, split, genes_dir, cluster_map_path, top_snps_dict):
    clusters = load_clusters(cluster_map_path)
    genes = os.listdir(genes_dir)
    # genes = genes[:500] ####
    # tracemalloc.start() ####
    data_lst = []
    for g in genes:
        if (top_snps_dict is not None) and( g not in top_snps_dict):
            continue
        gene_dir = os.path.join(genes_dir, g)
        plasma_path = os.path.join(gene_dir, run_name, "plasma_{0}.pickle")
        try:
            # print(os.path.getsize(plasma_path.format(split))) ####
            with open(plasma_path.format(split), "rb") as plasma_file:
                plasma_data = pickle.load(plasma_file)
        except (FileNotFoundError, pickle.UnpicklingError):
            continue

        data = read_data(plasma_data, clusters, g, top_snps=(top_snps_dict[g] if top_snps_dict is not None else None))
        data_lst.extend(data)

        # print(tracemalloc.get_traced_memory()) ####

    cols = [
        "Gene", 
        "Cluster", 
        "MeanTotalCoverage",
        "MeanTotalCoverageScaled",
        "MeanMappedCoverage",
        "MeanOverdispersion",
        "MeanCellCount",
        "UsableSampleSize",
        "TotalSampleSize",
        "UsableSNPCount",
        "TotalSNPCount",
        "CredibleSetSizeJoint", 
        "CredibleSetSizeAS",
        "CredibleSetSizeQTL",
        "CredibleSetPropJoint", 
        "CredibleSetPropAS",
        "CredibleSetPropQTL",
        "TopSNPPosterior",
        "TopSNPZComb",
        "TopSNPNLPComb",
        "TopSNPZPhi",
        "TopSNPPhi",
        "TopSNPNLPPhi",
        # "NumSigPhi",
        "TopSNPZBeta",
        "TopSNPBeta",
        "TopSNPNLPBeta",
        # "NumSigBeta",
        "TopSNPID",
        "Split",
    ]

    data_df = pd.DataFrame(data_lst, columns=cols)

    return data_df

def load_clusters(cluster_map_path):
    with open(cluster_map_path, "rb") as cluster_map_file:
        cluster_map = pickle.load(cluster_map_file)
    return cluster_map.keys()

def make_violin(
        df,
        var, 
        model_flavors,
        model_names, 
        model_colors,
        title, 
        result_path,
    ):
    sns.set(style="whitegrid", font="Roboto")
    plt.figure(figsize=(4,2))

    palette = [model_colors[m] for m in model_flavors]
    names = [model_names[m] for m in model_flavors]
    chart = sns.violinplot(
        x=var, 
        y="Model", 
        data=df, 
        order=model_flavors, 
        palette=palette,
        cut=0,
        scale="width"
    )
    ax = plt.gca()
    for art in ax.get_children():
        if isinstance(art, matplotlib.collections.PolyCollection):
            art.set_edgecolor((0., 0., 0.))
    plt.xlim(0., 1.)
    chart.set_yticklabels([model_names[m] for m in model_flavors])
    plt.ylabel("")
    plt.title(title)
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()

def make_thresh_barplot(
        df,
        var, 
        model_flavors,
        model_names, 
        threshs,
        title, 
        result_path,
    ):
    threshs = threshs + [np.inf]
    model_flavors = list(model_flavors)
    thresh_data = [[] for _ in range(len(threshs))]
    thresh_data_ret = [[] for _ in range(len(threshs))]
    for m in model_flavors:
        model_data = df.loc[df["Model"] == m, var].to_numpy()
        for i, t in enumerate(threshs):
            thresh_data[i].append(str(np.nanmean((model_data <= t).astype(int))))
            thresh_data_ret[i].append(np.nansum((model_data <= t).astype(int)))

    sns.set(style="whitegrid", font="Roboto")
    plt.figure(figsize=(4,2))
    palette = sns.cubehelix_palette(len(threshs), rot=-.25, light=.7)

    for i, t in enumerate(reversed(threshs)):
        estimator = lambda x: np.nanmean((x <= t).astype(int))
        # print(df[var]) ####
        # print(df[var].dtype) ####
        chart = sns.barplot(
            x=var, 
            y="Model", 
            data=df, 
            label=t, 
            order=model_flavors, 
            color=palette[i], 
            estimator=estimator,
            ci=None
        )
        plt.xlabel("Proportion of Loci")
        chart.set_yticklabels([model_names[m] for m in model_flavors])

    last_marker = [None for _ in range(len(model_flavors))]
    for i, t in enumerate(thresh_data[:-1]):
        # print(threshs[i]) ####
        for j, x in enumerate(t):
            xval = float(x)
            if (last_marker[j] is None and xval >= 0.04) or (last_marker[j] and (xval - last_marker[j]) >= 0.08):
                plt.text(
                    xval,
                    j,
                    threshs[i],
                    size="xx-small",
                    weight="medium",
                    ha="center",
                    va="center",
                    bbox={"boxstyle":"round", "pad":.25, "fc":"white", "ec":"white"}
                )
                last_marker[j] = xval

    plt.ylabel("")
    plt.title(title)
    # plt.legend()
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()

    return thresh_data_ret

def make_upset_plot(df, st, sn, thresh, clusters, result_path):
    df_select = df.loc[df[f"TopSNP{st}{sn}"] >= -np.log10(thresh), :]
    sets_data = {}
    # print(df_select) ####
    for index, row in df_select.iterrows():
        # print(i) ####
        sets_data.setdefault(row["Gene"], set()).add(row["Cluster"])
    counts = {}
    for key, value in sets_data.items():
        idx_tuple = tuple(c in value for c in clusters)
        counts.setdefault(idx_tuple, 0)
        counts[idx_tuple] += 1
    indices, values = zip(*counts.items())
    index = pd.MultiIndex.from_tuples(indices, names=clusters)
    # print(indices) ####
    # print(index) ####
    setcounts = pd.Series(values, index=index).sort_values(ascending=False)
    upsetplot.plot(setcounts)
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()

def calc_nlq(df, sn):
    data = df[f"TopSNPNLP{sn}"]
    count = np.count_nonzero(~np.isnan(data))
    ranks = np.argsort(np.nan_to_num(data))
    nlq = data - np.log10(count) + np.log10(ranks + 1) 
    min_sig = -np.inf
    # print(ranks) ####
    for i in ranks:
        sig = nlq.iloc[i]
        # print(i, sig) ####
        if np.isnan(sig):
            continue
        if sig > min_sig:
            min_sig = sig
        else:
            nlq.iloc[i] = min_sig
        # print(sig, nlq.iloc[i]) ####

    df[f"TopSNPNLQ{sn}"] = nlq

def plot_sets(df, out_dir):
    clusters = {
        "_all": "All Cells",
        "Ex": "Excitatory Neuron",
        "Oligo": "Oligodendrocyte",
        "Astro": "Astroglia",
        "In": "Inhibitory Neuron",
        "Endo": "Endothelial",
        "Microglia": "Microglia",
        "OPC": "Oligodendrocyte Progenitor",
        "Per": "Per"
    }
    model_map_dists = {
        "CredibleSetPropJoint": "PLASMA-J",
        "CredibleSetPropAS": "PLASMA-AS",
        "CredibleSetPropQTL": "QTL-Only"
    }
    model_map_thresh = {
        "CredibleSetSizeJoint": "PLASMA-J",
        "CredibleSetSizeAS": "PLASMA-AS",
        "CredibleSetSizeQTL": "QTL-Only"
    }
    model_map_thresh_prop = {
        "CredibleSetPropJoint": "PLASMA-J",
        "CredibleSetPropAS": "PLASMA-AS",
        "CredibleSetPropQTL": "QTL-Only"
    }
    var_dists = "Credible Set Proportion"
    var_thresh = "Credible Set Size"
    var_thresh_prop = "Credible Set Proportion"
    model_flavors_dists = model_map_dists.keys()
    model_flavors_thresh = model_map_thresh.keys()
    model_flavors_thresh_prop = model_map_thresh_prop.keys()
    pal = sns.color_palette()
    model_colors = {
        "CredibleSetPropJoint": pal[0],
        "CredibleSetPropAS": pal[4],
        "CredibleSetPropQTL": pal[7],
    }
    threshs = [5, 10, 20, 50, 100, 200]
    threshs_prop = [0.01, 0.05, 0.1, 0.2, 0.5]

    dfs_clust = {}

    for cluster in clusters.keys():
        # print(cluster)
        df_clust = df.loc[df["Cluster"] == cluster]
        calc_nlq(df_clust, "Comb")
        calc_nlq(df_clust, "Phi")
        calc_nlq(df_clust, "Beta")

        # print(np.count_nonzero(df_clust["TopSNPNLQComb"] >= -np.log10(0.1)))
        # print(np.count_nonzero(df_clust["TopSNPNLQPhi"] >= -np.log10(0.1)))
        # print(np.count_nonzero(df_clust["TopSNPNLQBeta"] >= -np.log10(0.1)))

        dfs_clust[cluster] = df_clust

        df_dists = pd.melt(
            df.loc[df["Cluster"] == cluster], 
            id_vars=["Gene"], 
            value_vars=model_map_dists.keys(),
            var_name="Model",
            value_name=var_dists
        )
        title = clusters[cluster]
        make_violin(
            df_dists,
            var_dists, 
            model_flavors_dists,
            model_map_dists, 
            model_colors,
            title, 
            os.path.join(out_dir, "sets_{0}.svg".format(cluster)),
        )
        df_thresh = pd.melt(
            df.loc[df["Cluster"] == cluster], 
            id_vars=["Gene"], 
            value_vars=model_map_thresh.keys(),
            var_name="Model",
            value_name=var_thresh
        )
        thresh_data = make_thresh_barplot(
            df_thresh,
            var_thresh, 
            model_flavors_thresh,
            model_map_thresh, 
            threshs,
            title, 
            os.path.join(out_dir, "thresh_{0}.svg".format(cluster)),
        )
        with open(os.path.join(out_dir, "threshdata_{0}.txt".format(cluster)), "w") as tdfile:
            tdfile.write("\t".join(map(str, threshs)) + "\n")
            tdfile.write("\t".join(model_flavors_thresh) + "\n\n")
            # print(thresh_data) ####
            np.savetxt(tdfile, np.array(thresh_data).astype(float))

        df_thresh_prop = pd.melt(
            df.loc[df["Cluster"] == cluster], 
            id_vars=["Gene"], 
            value_vars=model_map_thresh_prop.keys(),
            var_name="Model",
            value_name=var_thresh_prop
        )
        thresh_data_prop = make_thresh_barplot(
            df_thresh_prop,
            var_thresh_prop, 
            model_flavors_thresh_prop,
            model_map_thresh_prop, 
            threshs_prop,
            title, 
            os.path.join(out_dir, "thresh_{0}_prop.svg".format(cluster)),
        )
        with open(os.path.join(out_dir, "threshdata_{0}_prop.txt".format(cluster)), "w") as tdfile:
            tdfile.write("\t".join(map(str, threshs_prop)) + "\n")
            tdfile.write("\t".join(model_flavors_thresh_prop) + "\n\n")
            # print(thresh_data) ####
            np.savetxt(tdfile, np.array(thresh_data_prop).astype(float))

    df_all = pd.concat(dfs_clust.values())
    st = "NLQ"
    thresh = 0.1
    clusters = list(clusters.keys())
    make_upset_plot(df_all, st, "Comb", thresh, clusters, os.path.join(out_dir, f"upset_Comb_{st}_{thresh}.svg"))
    make_upset_plot(df_all, st, "Phi", thresh, clusters, os.path.join(out_dir, f"upset_Phi_{st}_{thresh}.svg"))
    make_upset_plot(df_all, st, "Beta", thresh, clusters, os.path.join(out_dir, f"upset_Beta_{st}_{thresh}.svg"))

    summ_data = []
    for cluster, df_clust in dfs_clust.items():
        data = [
            cluster,
            np.count_nonzero(df_clust["TopSNPNLQComb"] >= -np.log10(0.1)),
            np.count_nonzero(df_clust["TopSNPNLQPhi"] >= -np.log10(0.1)),
            np.count_nonzero(df_clust["TopSNPNLQBeta"] >= -np.log10(0.1))
        ]
        df_merged = pd.merge(
            df_clust, 
            dfs_clust["_all"], 
            on=["Gene"], 
            suffixes=["_clust", "_all"]
        )
        df_merged["TopSNPZ2CombDiff"] = df_merged["TopSNPZComb_clust"]**2 - df_merged["TopSNPZComb_all"]**2
        # threshs = [0.1, 0.2, 0.5]
        # print(df_merged.columns) ####
        df_specific = df_merged.loc[:,("Gene", "TopSNPZ2CombDiff", "TopSNPZComb_clust", "TopSNPZComb_all")].sort_values(by="TopSNPZ2CombDiff", ascending=False)
        # for ind, thresh in enumerate(threshs):
        #     cutoff = int(len(df_merged) * thresh)
        df_specific.to_csv(os.path.join(out_dir, "cell_type_spec", f"{cluster}.csv"), sep="\t", index=False, na_rep="None")   

        data_spec = [
            np.count_nonzero(np.logical_and(
                df_merged["TopSNPNLQComb_clust"] >= -np.log10(0.1),
                df_merged["TopSNPNLQComb_all"] < -np.log10(0.1)
            )),
            np.count_nonzero(np.logical_and(
                df_merged["TopSNPNLQPhi_clust"] >= -np.log10(0.1),
                df_merged["TopSNPNLQPhi_all"] < -np.log10(0.1)
            )),
            np.count_nonzero(np.logical_and(
                df_merged["TopSNPNLQBeta_clust"] >= -np.log10(0.1),
                df_merged["TopSNPNLQBeta_all"] < -np.log10(0.1)
            )),
        ]
        data.extend(data_spec)
        summ_data.append(data)

    cols = ["Cluster", "NumSigComb", "NumSigPhi", "NumSigBeta", "DiffSigComb", "DiffSigPhi", "DiffSigBeta"]
    summ_df = pd.DataFrame(summ_data, columns=cols)
    with open(os.path.join(out_dir, "sumstats.txt"), "w") as txt_file:
        summ_df.to_string(txt_file)

def make_scatter(
        df,
        var_x,
        var_y,
        var_h,
        x_label,
        y_label, 
        h_label,
        x_lim,
        y_lim,
        title, 
        result_path,
    ):
    df_rn = df.rename(columns={var_x: x_label, var_y: y_label, var_h: h_label})
    # print(df_rn) ####
    sns.set(style="whitegrid", font="Roboto")
    f, ax = plt.subplots(figsize=(5, 5))
    # ax.set(xscale="log", yscale="log")

    sns.scatterplot(
        x=x_label, 
        y=y_label, 
        hue=h_label,
        hue_norm=(0, 10),
        data=df_rn, 
    )
    if x_lim is not None:
        plt.xlim(-x_lim, x_lim)
    if y_lim is not None:
        plt.ylim(-y_lim, y_lim)
    plt.title(title)
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()

def plot_xval(df, out_dir):
    clusters = {
        "_all": "All Cells",
        "Ex": "Excitatory Neuron",
        "Oligo": "Oligodendrocyte",
        "Astro": "Astroglia",
        "In": "Inhibitory Neuron",
        "Endo": "Endothelial",
        "Microglia": "Microglia",
        "OPC": "Oligodendrocyte Progenitor",
        "Per": "Per"
    }
    for key, value in clusters.items():
        df_clust = df.loc[df["Cluster"] == key] 
        make_scatter(
            df_clust.loc[df["TopSNPNLPPhi_train"] >= -np.log10(0.05/df["UsableSNPCount_train"])],
            "TopSNPPhi_train",
            "TopSNPPhi_test",
            "TopSNPNLPPhi_test",
            "Train Effect Size",
            "Test Effect Size", 
            "Test -log10 P",
            None,
            None,
            "{0} AS Effect".format(value), 
            os.path.join(out_dir, "xval_phi_{0}.svg".format(key)),
        )
        make_scatter(
            df_clust.loc[df["TopSNPNLPBeta_train"] >= -np.log10(0.05/df["UsableSNPCount_train"])],
            "TopSNPBeta_train",
            "TopSNPBeta_test",
            "TopSNPNLPBeta_test",
            "Train Effect Size",
            "Test Effect Size",
            "Test -log10 P",
            None, 
            None,
            "{0} QTL Effect".format(value), 
            os.path.join(out_dir, "xval_beta_{0}.svg".format(key)),
        )

def make_heatmap(arr, order, title, result_path):
    heat_data = pd.DataFrame(data=arr, index=order, columns=order)
    sns.set(style="whitegrid", font="Roboto")
    f, ax = plt.subplots(figsize=(5, 5))
    sns.heatmap(heat_data, annot=True, fmt=".2g", square=True, cbar=False, annot_kws={"size": 10})
    plt.title(title)
    plt.savefig(result_path, bbox_inches='tight')
    plt.clf()

def plot_xcells(df_train, df_test, out_dir, stat_name_train, stat_name_test, cutoff_train, cutoff_test):
    sn1 = stat_name_train
    sn2 = stat_name_test
    df_tr_sig = df_train.loc[
        df_train[f"TopSNPNLP{sn1}"] >= -np.log10(0.05)
    ]
    # df_ts_sig = df_test.loc[
    #     np.logical_and(
    #         df_train["TopSNPNLPPhi"] >= -np.log10(0.05/df_train["UsableSNPCount"]),
    #         abs(df_train["TopSNPPhi"]) <= 5
    #     )
    # ]
    clusters = {
        "_all": "All Cells",
        "Ex": "Excitatory Neuron",
        "In": "Inhibitory Neuron",
        "Oligo": "Oligodendrocyte",
        "OPC": "Oligodendrocyte Progenitor",
        "Astro": "Astroglia",
        # "Endo": "Endothelial",
        # "Per": "Per"
    }
    cluster_order = list(clusters.keys())
    slopes = np.zeros((len(cluster_order), len(cluster_order),),)
    slope_ses = np.zeros((len(cluster_order), len(cluster_order),),)
    nlp_0s = np.zeros((len(cluster_order), len(cluster_order),),)
    nlp_1s = np.zeros((len(cluster_order), len(cluster_order),),)
    storey_pis = np.zeros((len(cluster_order), len(cluster_order),),)
    for ind_i, i in enumerate(cluster_order):
        for ind_j, j in enumerate(cluster_order):
            df_merged = pd.merge(
                df_tr_sig.loc[df_tr_sig["Cluster"] == i], 
                df_test.loc[df_test["Cluster"] == j], 
                on=["Gene"], 
                suffixes=["_train", "_test"]
            )

            num_sig_train = df_merged.shape[0]
            num_sig_test = np.sum(df_merged[f"TopSNPNLP{sn2}_test"] + np.log10(df_merged["UsableSNPCount_test"]) >= -np.log10(0.05))
            storey_pis[ind_i, ind_j] = num_sig_test / num_sig_train
            
    for ind_i, i in enumerate(cluster_order):
        for ind_j, j in enumerate(cluster_order):
            df_clust = df_train.loc[df_train["Cluster"] == i]
            # calc_nlq(df_clust, "Comb")
            calc_nlq(df_clust, "Phi")
            calc_nlq(df_clust, "Beta")
            df_merged = pd.merge(
                df_clust.loc[df_clust[f"TopSNPNLQ{sn1}"] >= -np.log(0.25)], 
                df_test.loc[df_test["Cluster"] == j], 
                on=["Gene"], 
                suffixes=["_train", "_test"]
            )

            make_scatter(
                df_merged,
                f"TopSNP{sn1}_train",
                f"TopSNP{sn2}_test",
                f"TopSNPNLP{sn2}_test",
                "Train Effect Size",
                "Test Effect Size", 
                "Test -log10 P",
                cutoff_train,
                cutoff_test,
                "{0} to {1}".format(clusters[i], clusters[j]), 
                os.path.join(out_dir, "xcell_{0}_{1}.svg".format(i, j)),
            )
            # print(i, j) ####
            # print(xw) ####
            # print(yw) ####

    title = "Cross-Cell Cross-Validation Storey Pi"
    make_heatmap(storey_pis, cluster_order, title, os.path.join(out_dir, "xcell_stats_pi.svg"))

def plot_xcells_nfold(dfs_train, dfs_test, out_dir, stat_name, cutoff):
    sn = stat_name
    clusters = {
        "_all": "All Cells",
        "Ex": "Excitatory Neuron",
        "In": "Inhibitory Neuron",
        "Oligo": "Oligodendrocyte",
        "OPC": "Oligodendrocyte Progenitor",
        "Astro": "Astroglia",
        # "Endo": "Endothelial",
        # "Per": "Per"
    }
    cluster_order = list(clusters.keys())
    slopes_all = []
    slopes_ses_all = []
    for df_train, df_test in zip(dfs_train, dfs_test):
        slopes = np.zeros((len(cluster_order), len(cluster_order),),)
        slope_ses = np.zeros((len(cluster_order), len(cluster_order),),)
        for ind_i, i in enumerate(cluster_order):
            for ind_j, j in enumerate(cluster_order):
                df_clust = df_train.loc[df_train["Cluster"] == i]
                # calc_nlq(df_clust, "Comb")
                calc_nlq(df_clust, stat_name)
                df_merged = pd.merge(
                    df_clust.loc[df_clust[f"TopSNPNLQ{sn}"] >= -np.log(0.25)], 
                    df_test.loc[df_test["Cluster"] == j], 
                    on=["Gene"], 
                    suffixes=["_train", "_test"]
                )
                x = df_merged[f"TopSNP{sn}_train"]
                y = df_merged[f"TopSNP{sn}_test"]
                se = df_merged[f"TopSNP{sn}_train"] / df_merged[f"TopSNPZ{sn}_train"]
                xw = np.nan_to_num(x / se)
                yw = np.nan_to_num(y / se)
                slope = xw.dot(yw) / xw.dot(xw)
                slopes[ind_i, ind_j] = slope
                res = yw - xw * slope
                slope_se = np.sqrt(res.dot(res) / (res.size * xw.dot(xw)))
                slope_ses[ind_i, ind_j] = slope_se
        slopes_all.append(slopes)
        slopes_ses_all.append(slope_ses)


    slopes_stacked = np.stack(slopes_all)
    slopes_ses_stacked = np.stack(slopes_ses_all)
    weights = 1 / slopes_ses_stacked**2
    norms = np.sum(weights, axis=0)
    slopes_meta = np.sum(slopes_stacked * weights, axis=0) / norms
    slopes_ses_meta = np.sqrt(1 / norms)
    z_0s = slopes_meta / slopes_ses_meta
    z_1s = (1 - slopes_meta) / slopes_ses_meta
    nlp_0s = -scipy.stats.norm.logsf(z_0s) / np.log(10)
    nlp_1s = -scipy.stats.norm.logsf(z_1s) / np.log(10)

    title = "Cross-Cell Cross-Validation Slopes"
    make_heatmap(slopes_meta, cluster_order, title, os.path.join(out_dir, "xcell_stats_slopes.svg"))

    title = "Cross-Cell Cross-Validation Slope Standard Errors"
    make_heatmap(slopes_ses_meta, cluster_order, title, os.path.join(out_dir, "xcell_stats_se.svg"))

    title = "Cross-Cell Cross-Validation Significance from Zero"
    make_heatmap(nlp_0s, cluster_order, title, os.path.join(out_dir, "xcell_stats_nlp_0.svg"))

    title = "Cross-Cell Cross-Validation Significance from One"
    make_heatmap(nlp_1s, cluster_order, title, os.path.join(out_dir, "xcell_stats_nlp_1.svg"))

def get_info(run_name, genes_dir, cluster_map_path, out_dir):
    data_df = make_df(run_name, "i0", genes_dir, cluster_map_path, None)
    data_df.sort_values(by=["TopSNPPosterior"], ascending=False, inplace=True)
    csv_path = os.path.join(out_dir, "cluster_info.csv")
    data_df.to_csv(csv_path, sep="\t", index=False, na_rep="None")
    txt_path = os.path.join(out_dir, "cluster_info.txt")
    with open(txt_path, "w") as txt_file:
        data_df.to_string(txt_file)
    plot_sets(data_df, out_dir)

def get_info_xval(run_name, num_splits, genes_dir, cluster_map_path, out_dir):
    df_train = make_df(run_name, "x0", genes_dir, cluster_map_path, None)
    top_snps_train = {}
    for index, row in df_train.iterrows():
        top_snps_train.setdefault(row["Gene"], {})[row["Cluster"]] = row["TopSNPID"]
    # print(top_snps_train) ####
    df_test = make_df(run_name, "x1", genes_dir, cluster_map_path, top_snps_train)
    df_comb = pd.merge(df_train, df_test, on=["Gene", "Cluster"], suffixes=["_train", "_test"])
    # print(df_train) ####
    # print(df_test) ####
    # print(df_comb) ####
    plot_xval(df_comb, os.path.join(out_dir, "xvals"))
    plot_xcells(df_train, df_test, os.path.join(out_dir, "xcells_phi"), "Phi", "Phi", 5, 5)
    plot_xcells(df_train, df_test, os.path.join(out_dir, "xcells_beta"), "Beta", "Beta", 30, 30)
    plot_xcells(df_train, df_test, os.path.join(out_dir, "xcells_phi_beta"), "Phi", "Beta", 5, 30)
    plot_xcells(df_train, df_test, os.path.join(out_dir, "xcells_beta_phi"), "Beta", "Phi", 30, 5)
    df_train.to_csv(os.path.join(out_dir, "train.csv"), sep="\t", index=False, na_rep="None")
    df_test.to_csv(os.path.join(out_dir, "test.csv"), sep="\t", index=False, na_rep="None")   

def get_info_xval_nfold(run_name, num_splits, genes_dir, cluster_map_path, out_dir):
    dfs_train = []
    dfs_test = []
    for i in range(num_splits):
        df_train = make_df(run_name, "x" + str(i), genes_dir, cluster_map_path, None)
        top_snps_train = {}
        for index, row in df_train.iterrows():
            top_snps_train.setdefault(row["Gene"], {})[row["Cluster"]] = row["TopSNPID"]
        # print(top_snps_train) ####
        df_test = make_df(run_name, "i" + str(i), genes_dir, cluster_map_path, top_snps_train)
        df_train.to_csv(os.path.join(out_dir, "nfold", "train_{0}.csv".format(i)), sep="\t", index=False, na_rep="None")
        df_test.to_csv(os.path.join(out_dir, "nfold" , "test_{0}.csv".format(i)), sep="\t", index=False, na_rep="None")
        dfs_train.append(df_train)
        dfs_test.append(df_test)
    # print(df_train) ####
    # print(df_test) ####
    # print(df_comb) ####
    plot_xcells_nfold(dfs_train, dfs_test, os.path.join(out_dir, "xcells_phi"), "Phi", 5)
    plot_xcells_nfold(dfs_train, dfs_test, os.path.join(out_dir, "xcells_beta"), "Beta", 50)
     

if __name__ == '__main__':
    data_path_kellis = "/agusevlab/awang/sc_kellis"
    # cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map.pickle")
    # genes_dir_kellis = os.path.join(data_path_kellis, "genes")

    # out_path_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis/cluster_info.txt"

    # get_info(genes_dir_kellis, cluster_map_path_kellis, out_path_kellis)

    cluster_map_path_kellis = os.path.join(data_path_kellis, "cluster_map_429.pickle")
    genes_dir_kellis = os.path.join(data_path_kellis, "genes_429")

    out_dir_kellis = "/agusevlab/awang/ase_finemap_results/sc_results/kellis_429"

    get_info("combined", genes_dir_kellis, cluster_map_path_kellis, out_dir_kellis)

    # get_info_xval("split", 2, genes_dir_kellis, cluster_map_path_kellis, out_dir_kellis)

    # get_info_xval_nfold("split5", 5, genes_dir_kellis, cluster_map_path_kellis, out_dir_kellis)

    # out_dir_strict = os.path.join(out_dir_kellis, "strict")

    # get_info_xval("split_strict", 2, genes_dir_kellis, cluster_map_path_kellis, out_dir_strict)


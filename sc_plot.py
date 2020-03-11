import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['agg.path.chunksize'] = 10000
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
import traceback

def plot_clusters(cells, cluster_selection, cluster_data, res_key, labels, barcode_map, proj_map, title, result_path):
	plot_data = []
	for cell in cells:
		cluster = barcode_map[cell]
		if cluster in cluster_selection:
			projection = proj_map[cell]
			val = cluster_data[res_key]
			plot_data.append([projection[0], projection[1], val])

	res_df = pd.DataFrame(dflst, columns=[label, primary_var_name, "Model"])

	sns.set(style="whitegrid", font="Roboto", rc={'figure.figsize':(4,4)})


	plt.title(title)
	plt.savefig(result_path, bbox_inches='tight')
	plt.clf()



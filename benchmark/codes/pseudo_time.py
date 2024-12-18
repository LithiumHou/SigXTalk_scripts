# Generate the pseudo-time 
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from pyslingshot import Slingshot
import sys

num_dims_reduced = 2
ng = int(sys.argv[1])
nb = int(sys.argv[2])
ncb = int(sys.argv[3])

cluster_labels = pd.read_csv('/home/jiawen/myMLnet/benchmark/datasets/sergio_anno.csv',sep=";")
cluster_labels = cluster_labels['cluster']
num_cells = len(cluster_labels)
labels = [item[7:] if item.startswith('Cluster') else item for item in cluster_labels]
labels = np.array([int(element) for element in labels])
start_node = 1

data = pd.read_csv('/home/jiawen/myMLnet/benchmark/datasets/sergio_coord.csv', header = 0, sep = "\t", index_col= 0, escapechar='\\')
data = data.values
from anndata import AnnData

ad = AnnData(np.zeros((nb*ncb, ng)))
ad.obsm["X_umap"] = data
ad.obs["celltype"] = labels
ad

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
custom_xlim = (-12, 12)
custom_ylim = (-12, 12)
# plt.setp(axes, xlim=custom_xlim, ylim=custom_ylim)

slingshot = Slingshot(ad, celltype_key="celltype", obsm_key="X_umap", start_node=start_node, debug_level='verbose')
slingshot.fit(num_epochs=1, debug_axes=axes)
pseudotime = pd.DataFrame(slingshot.unified_pseudotime)
pseudotime.index = ['Cell'+str(i) for i in range(num_cells)]
pseudotime.columns = ['Pseudotime1']
pseudotime.to_csv('/home/jiawen/myMLnet/benchmark/Beeline/inputs/example/bench_rectf/PseudoTime.csv', sep = ",")
pseudotime.to_csv('/home/jiawen/myMLnet/benchmark/Beeline/inputs/example/bench_tftg/PseudoTime.csv', sep = ",")
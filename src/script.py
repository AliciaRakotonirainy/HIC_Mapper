#imports
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import KernelPCA
from sklearn.preprocessing import StandardScaler
import argparse
from cover_complex import *
from scc import *

from gudhi import SimplexTree
from gudhi import CoverComplex

import networkx as nx


parser = argparse.ArgumentParser(description='Required to output a graph')

parser.add_argument('n_sample', type=int, nargs='+', default = 100,
    help='number of samples taken')

parser.add_argument('GAIN', type=float, nargs='+', default = 0.37,
    help='GAIN parameter for the Mapper')

parser.add_argument('RES', type=int, nargs='+', default = 5,
    help='RES parameter for the Mapper')

parser.add_argument('sample_replacement', type=bool, default = False,
    help='replacement in the cell pool while sampling')

parser.add_argument('data_path', type=str, default = "C:\Users\oliha\Documents\GitHub\HIC_Mapper\data",
    help="path to the file with two folders/files, 'processed' and 'data_features.txt")

args = parser.parse_args()

n_samples = args.n_sample
scc_mat = SCCMatrix(hic_matrices_folder= args.data_path +"/processed/",
    n_cells_sampling= 100, # number of cells 
    h=20, # window size for average smoothing
    n_slices_max=100,
    sample_repeat= args.sample_replacement,
    seed = 26 )

X_100 = scc_mat.compute_pairwise_dist()

cell_features = pd.read_table(args.data_path + + "\data_features.txt")
cell_features.index = cell_features["cell_nm"] # set cell names as index
mean_insu = cell_features["mean_insu"]

# cells that where selected :
cells_id = scc_mat.contact_maps_files # retrieve the name of the cells that were used to compute the pairwise distances
for i, cell in enumerate(cells_id): # correct the names (replace "." by "_")
    cells_id[i] = cell.replace(".","_")

mean_insu_cells = mean_insu.loc[cells_id]

# compute PCA of this distance matrix ; it will out filter (= lens)
X_100_scaled = StandardScaler().fit_transform(X_100)
transformer = KernelPCA(n_components=2, kernel='rbf') # we could try other kernels
pca_100 = transformer.fit_transform(X_100_scaled)



RES = args.RES
GAIN = args.GAIN
mapper = MapperComplex(input_type="distance matrix", 
    colors=mean_insu_cells,
    resolutions=np.array([RES, RES]), gains=np.array([GAIN, GAIN]),
    clustering=AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=2.5, affinity="euclidean"))

mapper = mapper.fit(X_100, filters=pca_100) # filters : numpy array of shape (num_points) x (num_filters) Each column of the numpy array defines a scalar function defined on the input points.

# plot graph    
G = mapper.get_networkx()
plt.figure()
nx.draw_networkx(G,
    node_color=[mapper.node_info_[v]["colors"] for v in G.nodes()])
plt.plot()   
plt.close()

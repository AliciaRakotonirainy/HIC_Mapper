import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import KernelPCA
from sklearn.preprocessing import StandardScaler


from src.cover_complex import *
from src.scc import *

from gudhi import SimplexTree
from gudhi import CoverComplex

import networkx as nx

# test the mapper on simple example :

## compute the mapper graph associated to this point cloud : 

np.random.seed(27)
t = np.linspace(0,10,300)
x_circle = np.sin(t)
y_circle = np.cos(t) + np.random.randn(len(t))/10

x_circle_2 = np.sin(t)
y_circle_2 = np.cos(t) -3 + np.random.randn(len(t))/10

x_center = np.random.randn(400)/10
y_center = np.random.randn(400)/3 - 1.5

x = np.concatenate((x_circle, x_circle_2, x_center))
y = np.concatenate((y_circle, y_circle_2, y_center))


X = np.concatenate((x.reshape(-1,1),y.reshape(-1,1)), axis=1)
X.shape

# filter = projection on y axis
X_filter = X[:,1]

mapper = MapperComplex(input_type="point cloud",
                       resolutions=np.array([7]))
mapper.fit(X, filters=X_filter)
G = mapper.get_networkx()

# we obtain the expected graph : 2 linked circles !

fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Test mapper on a simple example')
ax1.scatter(x, y)
ax1.set_title("Data points")
nx.draw_networkx(G, ax = ax2)
ax2.set_title("Mapper Graph")
plt.show()


## Compute the mapper on hi-c dataset

# Compute pairwise distance matrix

scc_mat = SCCMatrix("data/chromosomes.txt",
                    "data/hi-matrices/",
                    100,
                    h=10,
                    n_slices_max=10)
X = scc_mat.compute_pairwise_scc()

# plot the pairwise distance matrix
plt.imshow(X)
plt.colorbar()
plt.show()

# Get colors 
cell_features = pd.read_table("data/features.txt")
cell_features.index = cell_features["cell_nm"]
mean_insu = cell_features["mean_insu"]

# cells that where selected :
cells_id = scc_mat.contact_maps_files
for i, cell in enumerate(cells_id):
    cells_id[i] = cell.replace(".","_")

mean_insu_cells = mean_insu.loc[cells_id]

# compute PCA of this distance matrix
X_scaled = StandardScaler().fit_transform(X)
transformer = KernelPCA(n_components=2, kernel='linear')
pca = transformer.fit_transform(X_scaled)

# compute mapper, using 2 first Principal Components as filter
mapper = MapperComplex(input_type="distance matrix", 
                       colors=mean_insu_cells,
                       C=5, resolutions=np.array([5,5]), gains=np.array([0.6,0.6]))
mapper = mapper.fit(X, filters=pca) # filters : numpy array of shape (num_points) x (num_filters) Each column of the numpy array defines a scalar function defined on the input points.

# plot graph
G = mapper.get_networkx()
plt.figure()
nx.draw_networkx(G,
                node_color=[mapper.node_info_[v]["colors"] for v in G.nodes()])
plt.show()

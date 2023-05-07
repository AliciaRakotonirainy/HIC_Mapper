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
import logging
from tqdm import tqdm

logging.basicConfig()
logger = logging.root
logFormatter = logging.Formatter('{relativeCreated:12.0f}ms {levelname:5s} [{filename}] {message:s}', style='{')
logger.setLevel(logging.INFO)

# get CL arguments
parser = argparse.ArgumentParser(description='Required to output a graph')

parser.add_argument('--n_sample', type=int, default = 100,
    help='number of samples taken')
parser.add_argument('--GAIN', type=float, default = 0.37,
    help='GAIN parameter for the Mapper')
parser.add_argument('--RES', type=int, default = 5,
    help='RES parameter for the Mapper')
parser.add_argument('--sample_replacement', type=bool, default = False,
    help='replacement in the cell pool while sampling')
parser.add_argument('--data_path', type=str, 
    help="path to the data folder. Should end with '/'. It should contain the folder 'hic-matrices', with the folders containing each contact map file. It should also contain the 'features.txt' file.")
parser.add_argument('--output_path', type=str, 
    help="path to folder where the mapper graph will be saved. Should end with '/' .")

args = parser.parse_args()

def main():
    N_SAMPLE = int(args.n_sample)
    DATA_PATH = args.data_path
    OUTPUT_PATH = args.output_path
    SAMPLE_REPLACEMENT = args.sample_replacement
    RES = float(args.RES)
    GAIN = float(args.GAIN)

    logger.info("Loading data...")

    # Initialize the SCCMatrix object (loads the contact maps of N_SAMPLE randomly selected cells)
    scc_mat = SCCMatrix(hic_matrices_folder= DATA_PATH +"hi-matrices/",
        n_cells_sampling= N_SAMPLE, # number of cells 
        h=20, # window size for average smoothing
        n_slices_max=100,
        sample_repeat= SAMPLE_REPLACEMENT,
        seed = 26)
    
    logger.info("Computing pairwise SCC matrix...")
    # Computes the pairwise SCC matrix and the pairwise distance matrix between contact maps
    X = scc_mat.compute_pairwise_dist()

    # load cell cycle features
    cell_features = pd.read_table(DATA_PATH + "features.txt")
    cell_features.index = cell_features["cell_nm"] # set cell names as index
    cell_features_names = ["mean_insu", "f_near_band", "f_mitotic_band", "repli_score"]

    # cells that where selected :
    cells_id = scc_mat.contact_maps_files # retrieve the name of the cells that were used to compute the pairwise distances
    for i, cell in enumerate(cells_id): # correct the names (replace "." by "_")
        cells_id[i] = cell.replace(".","_")

    logger.info("Computing PCA of pairwise SCC matrix...")
    # compute PCA of this distance matrix 
    X_scaled = StandardScaler().fit_transform(X)
    transformer = KernelPCA(n_components=2, kernel='rbf') # we could try other kernels
    pca = transformer.fit_transform(X_scaled)

    logger.info("Computing mapper...")
    # plot mapper graph colored with 4 features that are related to cell cycle
    for f_name in cell_features_names:
        f_name_cells = cell_features[f_name]
        f_name_cells = f_name_cells.loc[cells_id]

        mapper = MapperComplex(input_type="distance matrix", 
            colors=f_name_cells,
            resolutions=np.array([RES, RES]), gains=np.array([GAIN, GAIN]),
            clustering=AgglomerativeClustering(n_clusters=None, linkage="single", distance_threshold=2.5, affinity="euclidean"))

        mapper = mapper.fit(X, filters=pca) # filters : numpy array of shape (num_points) x (num_filters) Each column of the numpy array defines a scalar function defined on the input points.

        # plot mapper graph  
        G = mapper.get_networkx()
        plt.figure()
        nx.draw_networkx(G,
            node_color=[mapper.node_info_[v]["colors"] for v in G.nodes()])
        plt.plot()   
        plt.title(f"Mapper graph - Colored by {f_name} \n N_SAMPLE={N_SAMPLE}, RES={RES}, GAIN={GAIN}")
        plt.savefig(OUTPUT_PATH + f"mapper_{f_name}_{N_SAMPLE}_{RES}_{GAIN}.png")
        plt.close()

    logger.info(f"Mapper saved in {OUTPUT_PATH} folder! \n DONE!")

if __name__ == "__main__":
    main()
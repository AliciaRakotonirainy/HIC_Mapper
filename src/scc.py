import numpy as np
import scipy as sc
import os
import cv2 as cv
import matplotlib.pyplot as plt
import pandas as pd
import random
from tqdm import tqdm

class SCCMatrix():

    def __init__(self,
                 chromosomes_path,
                 hic_matrices_folder,
                 n_cells_sampling,
                 h):
        """ 
        This class computes the pairwise SCC matrix of specified cells.
        """
        
        self.chromosomes_path = chromosomes_path
        self.hic_matrices_folder = hic_matrices_folder
        self.n_cells_sampling = n_cells_sampling
        self.h = h
        
        self.contact_maps_files = None
        self.pairwise_scc_matrix = None

    def load_data(self):
        self.chromosomes = pd.read_table(self.chromosomes_path, header=None)
        self.chromosomes.columns = ["chr", "start", "end"]

        # randomly sample the desired number of cells
        self.contact_maps_files = random.sample(
                                        os.listdir(self.hic_matrices_folder), 
                                        self.n_cells_sampling)
        
    def smooth_average(self, mat, h):
        """ 
        Average smoothing of provided matrix, with window size h
        """
        return cv.blur(mat, (h,h))
    
    
    def compute_slices(self, mat):
        """
        Computes the strata of provided contact matrix
        """
        n = mat.shape[0]
        slices = []
        # we create n slices, with n the number of windows in the data matrix
        for k in range(n-1):
            # the slice k contains all the contacts (i,j) such that abs(j-i) = k
            # ie in slice k, all the contacts are made within [k*b, (k+1)*b] of genomic distance 
            # we take only slice with at least 2 elements
            slices.append(
                np.array([mat[i, i+k] for i in range(n-k)]) 
            )

        return slices
    
    
    def compute_scc(self, mat_1, mat_2, h):
        """ 
        Computes the SCC between 2 contact matrices, average smoothing them with window size h
        """
        # smooth the input matrices
        mat_1_smooth = self.smooth_average(mat_1, h)
        mat_2_smooth = self.smooth_average(mat_2, h)

        # compute slices for each matrix
        slices_1 = self.compute_slices(mat_1_smooth)
        slices_2 = self.compute_slices(mat_2_smooth)
        
        # check we have the same number of slices in both matrices
        assert len(slices_1) == len(slices_2)
        K = len(slices_1)

        SCC = 0
        norm_factor = 0
        for k in range(K):
            Xk = slices_1[k]
            Yk = slices_2[k]

            # we cannot compute pearson cor in this strata when there is no contact
            if sum(Xk) == 0 or sum(Yk) == 0:
                continue

            Nk = len(slices_1[k])
            r2k = np.std(Xk) * np.std(Yk)
            pearson_k = sc.stats.pearsonr(Xk, Yk)[0]
            SCC += Nk*r2k*pearson_k 
            norm_factor += Nk * r2k
        SCC = SCC / norm_factor
        return SCC
    
    def compute_pairwise_scc(self):
        
        self.pairwise_scc_matrix = np.zeros(self.n_cells_sampling,
                                            self.n_cells_sampling)

        for i in tqdm(range(self.n_cells_sampling)):
            for j in range(i, self.n_cells_sampling):
                mat_i = sc.sparse.load_npz(self.contact_maps_files[i] + "cmatrix_500k.npz") # shape 5 234 x 5 234
                mat_j = sc.sparse.load_npz(self.contact_maps_files[j] + "cmatrix_500k.npz")

                mat_i_arr = mat_i.toarray()
                mat_j_arr = mat_j.toarray()

                scc = self.compute_scc(mat_i_arr, mat_j_arr, self.h)

                self.pairwise_scc_matrix[i,j] = scc
                self.pairwise_scc_matrix[j,i] = scc


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
                 h,
                 n_slices_max):
        """ 
        This class computes the pairwise SCC matrix of specified cells.
        """
        
        self.chromosomes_path = chromosomes_path
        self.hic_matrices_folder = hic_matrices_folder
        self.n_cells_sampling = n_cells_sampling
        self.h = h
        self.n_slices_max = n_slices_max
        
        self.contact_maps_files = None
        self.pairwise_scc_matrix = None
        self.pairwise_distance_matrix = None

        self.smooth_matrices = [] # [ [smoothed matrix cell 1 chr1, smoothed matrix cell 1 chr2, ..., smoothed matrix cell 1 chr21], 
                                #     [smoothed matrix cell 2 chr1, smoothed matrix cell 2 chr2, ..., smoothed matrix cell 2 chr21],
                                #     [ ... ],
                                #     [smoothed matrix cell n chr1, smoothed matrix cell n chr2, ..., smoothed matrix cell n chr21]]
        self.slices_matrices = [] # [ [slices cell 1 chr 1, slices cell 1 chr2, ..., slices cell 1 chr21 ],
                                #     [slices cell 2 chr 1, slices cell 2 chr2, ..., slices cell 2 chr21 ],
                                #     [ ... ],
                                #     [slices cell n chr 1, slices cell n chr2, ..., slices cell n chr21 ]

        self.load_data()

    def load_data(self):
        self.chromosomes = pd.read_table(self.chromosomes_path, header=None)
        self.chromosomes.columns = ["chr", "start", "end"]

        # randomly sample the desired number of cells
        self.contact_maps_files = random.sample(
                                        os.listdir(self.hic_matrices_folder), 
                                        self.n_cells_sampling)
        
    def separate_chromosomes(self, mat):
        mat_chrom = []
        for i in range(len(self.chromosomes)):
            mat_chrom.append(mat[
                                self.chromosomes.iloc[i, 1]:self.chromosomes.iloc[i, 2]+1, 
                                self.chromosomes.iloc[i, 1]:self.chromosomes.iloc[i, 2]+1 # +1 because excluded
                                ]
                            )      
        return mat_chrom 
        
    def smooth_average(self, mat, h):
        """ 
        Average smoothing of provided matrix, with window size h
        """
        return cv.blur(mat, (h,h))
    
    
    def compute_slices(self, mat):
        """
        Computes the strata of provided contact matrix
        """
        slices = []
        # we create n slices, with n the number of windows in the data matrix
        n = min(self.n_slices_max, mat.shape[0])
        for k in range(n - 1):
            # the slice k contains all the contacts (i,j) such that abs(j-i) = k
            # ie in slice k, all the contacts are made within [k*b, (k+1)*b] of genomic distance 
            # we take only slice with at least 2 elements
            slices.append(
                np.array([mat[i, i+k] for i in range(n-k)]) 
            )

        return slices
    
    
    def compute_scc(self, slices_1, slices_2, h):
        """ 
        TO BE UPDATED
        """
        # we define the scc as the mean of the intra-chromosomes scc
        scc = []
        for ch in range(len(slices_1)):
            K = len(slices_1[ch])
            num = sum([len(slices_1[ch][k]) * np.cov(
                    np.concatenate(
                        (slices_1[ch][k].reshape(1, -1),
                        slices_2[ch][k].reshape(1, -1))
                        )
                    )[0,1]
            for k in range(K)
            ])
            deno = sum([len(slices_1[ch][k]) * np.std(slices_1[ch][k]) * np.std(slices_2[ch][k])
                    for k in range(K)
            ])
            if deno != 0:
                scc.append(num/deno)

        return np.mean(scc)
    
    def smooth_all(self,h):
        for i in range(self.n_cells_sampling):
            mat_i = sc.sparse.load_npz(self.hic_matrices_folder + self.contact_maps_files[i] + "/cmatrix_500k.npz")
            mat_i_arr = mat_i.toarray()
            mat_i_chrom = self.separate_chromosomes(mat_i_arr)

            mat_i_chrom_smooth = []
            for ch in range(len(mat_i_chrom)):
                mat_i_chrom_smooth.append( self.smooth_average(mat_i_chrom[ch], h) )
            self.smooth_matrices.append(mat_i_chrom_smooth)

    def slices_all(self):
        for i in range(self.n_cells_sampling):
            mat_i_slices = []
            for ch in range(len(self.smooth_matrices[i])):
                mat_i_slices.append( self.compute_slices(self.smooth_matrices[i][ch]) )
            self.slices_matrices.append( 
                                mat_i_slices
                                )

    def compute_pairwise_scc(self):

        self.smooth_all(self.h)
        self.slices_all()
        
        self.pairwise_scc_matrix = np.zeros((self.n_cells_sampling,
                                            self.n_cells_sampling))
        
        self.pairwise_distance_matrix = np.zeros((self.n_cells_sampling,
                                            self.n_cells_sampling))
        
        # first compute diagonal terms
        for i in range(self.n_cells_sampling):
            scc = self.compute_scc(self.slices_matrices[i], self.slices_matrices[i], self.h)
            self.pairwise_scc_matrix[i,i] = scc
            self.pairwise_distance_matrix[i,i] = 0

        for i in tqdm(range(self.n_cells_sampling)):
            for j in range(i+1, self.n_cells_sampling):
                scc = self.compute_scc(self.slices_matrices[i], self.slices_matrices[j], self.h)

                self.pairwise_scc_matrix[i,j] = scc
                self.pairwise_scc_matrix[j,i] = scc

                dist = np.sqrt(self.pairwise_scc_matrix[i,i] + self.pairwise_scc_matrix[j,j] - 2*scc)

                self.pairwise_distance_matrix[i,j] = dist
                self.pairwise_distance_matrix[j,i] = dist

        return self.pairwise_distance_matrix

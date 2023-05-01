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
                 hic_matrices_folder,
                 n_cells_sampling,
                 h,
                 n_slices_max,
                 sample_repeat = False,
                 seed=27):
        """ 
        This class computes the pairwise SCC and distance matrices.

        :hic_matrices_folder: Folder that contains all the folders with cells contact maps. Should end with "/"
        :n_cells_sampling: Number of cells to take in the analysis
        :h: Size of window for average smoothing
        :n_slices_max: maximal distance between 2 DNA bins in order to take the contact into account in the strata
        :seed: Random seed, to reproduce cell sampling
        """
        
        self.hic_matrices_folder = hic_matrices_folder
        self.n_cells_sampling = n_cells_sampling
        self.h = h
        self.n_slices_max = n_slices_max
        self.seed = seed
        self.sample_repeat = sample_repeat
        
        self.contact_maps_files = None
        self.pairwise_scc_matrix = None
        self.pairwise_distance_matrix = None

        self.smooth_matrices = [] # list of smoothed contact matrices
        self.slices_matrices = [] # list of list, containg the slices of each smoothed matrix

        self.load_data()

    def load_data(self):
        # randomly sample the desired number of cells
        np.random.seed(self.seed)
        self.contact_maps_files = []
        if self.sample_repeat == True:
            for x in range(self.n_cells_sampling):
                self.contact_maps_files.extend(random.sample(
                                        os.listdir(self.hic_matrices_folder), 
                                        1))
        else: 
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
        Computes the slices (= strata) of provided contact matrix
        """
        slices = []
        # we take into account only the contacts within self.n_slices_max bins of distance
        for k in range(self.n_slices_max-1):
            # the slice k contains all the contacts (i,j) such that abs(j-i) = k
            # ie in slice k, all the contacts are made within [k*b, (k+1)*b] of genomic distance 
            slices.append(
                np.array([mat[i, i+k] for i in range(self.n_slices_max-k)]) 
            )

        return slices
    
    
    def compute_scc(self, slices_1, slices_2, h):
        """ 
        Computes the SCC between 2 contact matrices, given their respective slices
        """
        
        # check we have the same number of slices in both matrices
        assert len(slices_1) == len(slices_2)
        K = len(slices_1)

        num = sum([len(slices_1[k]) * np.cov(
                    np.concatenate(
                        (slices_1[k].reshape(1, -1),
                        slices_2[k].reshape(1, -1))
                        )
                    )[0,1]
        for k in range(K)
        ])
        deno = sum([len(slices_1[k]) * np.std(slices_1[k]) * np.std(slices_2[k])
                for k in range(K)
        ])

        return num/deno
    
    def smooth_all(self,h):
        """ 
        Smooth all the contact matrices
        """
        for i in range(self.n_cells_sampling):
            mat_i = sc.sparse.load_npz(self.hic_matrices_folder + self.contact_maps_files[i] + "/cmatrix_500k.npz")
            mat_i_arr = mat_i.toarray()
            mat_i_smooth = self.smooth_average(mat_i_arr, h)
            self.smooth_matrices.append(mat_i_smooth)

    def slices_all(self):
        """
        Compute the slices (= strata) of all the contact matrices
        """
        for i in range(self.n_cells_sampling):
            self.slices_matrices.append( 
                                self.compute_slices(self.smooth_matrices[i])
                                )

    def compute_pairwise_dist(self):
        """
        Computes the pairwise SCC matrix and the pairwise distance matrix, for all the cells.
        The distance between 2 matrices X, Y is defined as : np.sqrt(SCC(X,X) + SCC(Y,Y) - 2*SCC(X,Y))
        """

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

        # then compute all the other terms
        for i in tqdm(range(self.n_cells_sampling)):
            for j in range(i+1, self.n_cells_sampling):
                scc = self.compute_scc(self.slices_matrices[i], self.slices_matrices[j], self.h)

                self.pairwise_scc_matrix[i,j] = scc
                self.pairwise_scc_matrix[j,i] = scc

                dist = np.sqrt(self.pairwise_scc_matrix[i,i] + self.pairwise_scc_matrix[j,j] - 2*scc)

                self.pairwise_distance_matrix[i,j] = dist
                self.pairwise_distance_matrix[j,i] = dist

        return self.pairwise_distance_matrix
    

######### Same class, but computing the SCC intra-chromosomally first, then averaging the SCC of each chromosomes : ########

class SCCMatrix_chromwise():

    def __init__(self,
                 chromosomes_path,
                 hic_matrices_folder,
                 n_cells_sampling,
                 h,
                 n_slices_max,
                 seed=27):
        """ 
        This class computes the pairwise SCC and distance matrices. Difference with SCCMatrix() : it 
        computes the SCC values intra-chromosommally. The SCC between 2 cells is defined as the average
        SCC of their corresponding chromosomes.

        :chromosomes_path: Path to file containing the chromosomes info
        :hic_matrices_folder: Folder that contains all the folders with cells contact maps. Should end with "/"
        :n_cells_sampling: Number of cells to take in the analysis
        :h: Size of window for average smoothing
        :n_slices_max: maximal distance between 2 DNA bins in order to take the contact into account in the strata
        :seed: Random seed, to reproduce cell sampling
        """
        
        self.chromosomes_path = chromosomes_path
        self.hic_matrices_folder = hic_matrices_folder
        self.n_cells_sampling = n_cells_sampling
        self.h = h
        self.n_slices_max = n_slices_max
        self.seed = seed
        
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
        np.random.seed(self.seed)
        self.contact_maps_files = random.sample(
                                        os.listdir(self.hic_matrices_folder), 
                                        self.n_cells_sampling)
        
    def separate_chromosomes(self, mat):
        """
        For a given contact matrices, returns the list of the sub matrices corresponding to each chromosome.
        """
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
        Computes the slices (= strata) of provided contact matrix
        """
        slices = []
        # we take into account only the contacts within self.n_slices_max bins of distance
        n = min(self.n_slices_max, mat.shape[0])
        for k in range(n - 1):
            # the slice k contains all the contacts (i,j) such that abs(j-i) = k
            # ie in slice k, all the contacts are made within [k*b, (k+1)*b] of genomic distance 
            slices.append(
                np.array([mat[i, i+k] for i in range(n-k)]) 
            )

        return slices
    
    
    def compute_scc(self, slices_1, slices_2, h):
        """ 
        Computes the SCC between 2 matrices, given their slices.
        It first computes the SCC between their corresponding chromosomes, then it averages 
        the SCC of all chromosomes.
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
        """
        Smooth all the contact matrices
        """
        for i in range(self.n_cells_sampling):
            mat_i = sc.sparse.load_npz(self.hic_matrices_folder + self.contact_maps_files[i] + "/cmatrix_500k.npz")
            mat_i_arr = mat_i.toarray()
            mat_i_chrom = self.separate_chromosomes(mat_i_arr)

            mat_i_chrom_smooth = []
            for ch in range(len(mat_i_chrom)):
                mat_i_chrom_smooth.append( self.smooth_average(mat_i_chrom[ch], h) )
            self.smooth_matrices.append(mat_i_chrom_smooth)

    def slices_all(self):
        """
        Compute the slices (= strata) of all the contact matrices
        """
        for i in range(self.n_cells_sampling):
            mat_i_slices = []
            for ch in range(len(self.smooth_matrices[i])):
                mat_i_slices.append( self.compute_slices(self.smooth_matrices[i][ch]) )
            self.slices_matrices.append( 
                                mat_i_slices
                                )

    def compute_pairwise_dist(self):
        """
        Computes the pairwise SCC matrix and the pairwise distance matrix, for all the cells.
        The distance between 2 matrices X, Y is defined as : np.sqrt(SCC(X,X) + SCC(Y,Y) - 2*SCC(X,Y))
        """

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

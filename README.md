# Code reproducibility

Clone this repository.

Create a new conda environment : 

```bash
conda create --name gmda --file requirements.txt
conda activate gmda
```

Go in the `src/` directory, and launch this command:

```bash
python main.py --n_sample=100 --GAIN=0.37 --RES=5 --sample_replacement=False --data_path="../data/" --output_path="../output/"
```

The mapper graphs (colored with the 3 markers related with cell cycle) are saved in the `output/` folder at the root directory, as `.png` images. The name of the images is of the form: `mapper_${marker_name}_${n_sample}_${RES}_${GAIN}.png`

`--n_sample` : int, Number of cells drawn to do the analysis

`--sample_replacement` : bool, Sample cells with replacement

`--GAIN` : float, Gain parameter (Mapper)

`--RES` : int, Resolution parameter (Mapper)

`--data_path` : str, Path to the data folder. Should end with '/'. It should contain the folder 'hic-matrices', with the folders containing each contact map file. It should also contain the 'features.txt' file.

`--output_path` : str, Path to folder where the mapper graph will be saved. Should end with '/' .

# Context

In this project, we analyzed a data set of single-cell Hi-C contact maps.

Single-cell Hi-C (scHi-c) can quantifies the three-dimensional chromatin organization. These 3D genome features are related with vital genome functionalities. scHi-c experiments can be summarized into so-called "Hi-c contact maps". These are pairwise distance matrices that encode how chromatin is folded in the nucleus of a cell. More precisely, each row and column of the matrix represents a small DNA window, and each entry in the matrix is the spatial distance between these windows in the nucleus.

Yang et al (2017) proposed a new similarity measure called Stratum-Adjusted Correlation Coefficient (SCC). It quantifies the similarity between Hi-C contact maps. It can be used to assess the reproducibility of replicate samples, as well as quantify the distance between Hi-C matrices from different cell types or conditions. This is this latest application that we will consider in this project.

The 3D-structure of DNA varies with the cell cycle. For instance, chromosomes tend to fold up before mitosis, and then unfold in order to enable DNA transcription. The goal of this project was to detect the cell cycle as a big loop in the Mapper Graph, based on contact maps processed using SCC. We first computed the pairwise SCC matrix (that quantifies the distance between contact maps of all pairs of cells). Then, we found a set of parameters for which the cell cycle loop can be detected in the Mapper graph. Then, we proved that this loop indeed represents the cell cycle by coloring the Mapper nodes with various markers correlated with the cell cycle. 

### 0. Install required packages
# Install Seurat
if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages("Seurat")

# Install dplyr
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

# Install spSeudoMap
if (!requireNamespace("dplyr", quietly = TRUE))
  devtools::install_github("spSeudoMap")



### 1. Example of using function "pred_cellf_celldart"
library(SeuratObject)

# Find the directory for active script file (file_path)
file_path <- rstudioapi::getSourceEditorContext()$path
file_path <- strsplit(file_path, split=.Platform$file.sep)
file_path <- paste(file_path[[1]][-length(file_path[[1]])],
                   collapse=.Platform$file.sep)

# Set working directory
setwd(file_path)

# Import the CellDART function
source('spSeudoMap_R.R')

# Make output folder
output_folder_name <- 'CellDART_output'
if (!file.exists(output_folder_name)){
  dir.create(output_folder_name)
}


## Load single-cell and spatial datasets
# Load single-cell dataset (RDS file with Seurat object): GSE115746 (mouse single cell: ALS and VISp)
sc_data <- readRDS('sc_data.rds')

# Load spatial dataset (RDS file with Seurat object): 10X genomics data repository
# V1_Mouse_Brain_Sagittal_Anterior
# V1_Mouse_Brain_Sagittal_Posterior
sp_data <- readRDS('sp_data.rds')


## Check the size of spatial dataset
dim(sp_data)

## Set the number of pseudospots: 5 times the number of spatial spots
npseudo <- 5*dim(sp_data)[2]



## Perform spSeuoMap analysis
# cell type mapping of spatial transcriptomics using unmatched single-cell RNA-seq data

# Seurat object to AnnData conversion algorithm: reference to sceasy
# https://github.com/cellgeni/sceasy

## Input

# adata_sp: spatial data (AnnData object) to be used in predicting cell fraction (default: None)
# -> count matrix should be non-normalied raw data
# adata_sc: single-cell data (AnnData object) to be used in making pseudospots (default: None)
# -> count matrix should be non-normalied raw data

# outdir: the directory to save output files (models and results) (default = '.')
# source.code.dir: directory of python source code (default = 'pred_cellf_spSeudoMap.py')

# sp_subset: whether to subset spatial data and calculate for specific spot cluster (default = FALSE)
# spot.cluster.name: group name of the cluster used for subsetting spatial data (default: 'seurat_clusters')
# spot.cluster.of.interest: name of each spot clusters to be used (default: NULL)
# metadata_celltype: column name for single-cell annotation data in metadata (default: 'celltype')

# env.select: select between using reticulate virtual environment or conda environment ("virtual" or "conda")
# python.install: whether to automatically install python version 3.7.12

# python_path: path for the python 3.7. (default: NULL)
# If NULL, python version 3.7.12 will be installed (valid for Linux) 
# If "current", python interpreter associated with current virtual env (ex: r-reticulate) will be used. (version should be 3.7)
# env.name: name of the virtual or conda environment to use for CellDART analysis (default: 'spSeudoMap')

# gpu: check whether to use gpu (True) or not (False) (default = True)
# metadata_celltype: column name for single-cell annotation data in metadata (default: 'celltype')
# num_markers: number of selected marker genes in each cell-type (default = 20)

# mixture_mode: mode of the pseudospot generation 
# -> 'default': when cell types are similar between single-cell/spatial data (Identical to CellDART)
# -> 'pseudotype': when there is a mismatch of cell types between single-cell/spatial data and cell types exclusively present in spatial data is considered when generating pseudospots

# seed_num: seed to be used in random sampling (default = 0)

# mk_ratio_fix: whether to fix the mk_ratio when selecting the 
# mk_ratio: ratio of number of single-cell markers to virtual pseudotype markers (default = 2)
# pseudo_num_genes: number of the virtual markers genes for pseudotype (default = 40)
# -> the number is used only when the mk_ratio is not fixed and number of pseudotype markers should be manually provided 

# pseudo_frac_m: average of presumed fraction of the pseudotype (cell types exclusively present in spatial data) across all spatial spots (default = 0.5)
# -> determined by cell sorting study or by literature based evidence that enables speculation of average fraction in the tissue
# pseudo_frac_std: standard deviation of the distribution of presumed pseudotype fraction across all spatial spots (default = 0.1)
# num_top_genes: number of top genes having highest log fold change between spatial and single-cell normalized pseudobulk counts (spatial/single-cell)
# -> the top genes are used for generating module score to predict pseudotype fraction in spatial spots

# nmix: the number of cells sampled from single-cell data when making a pseudospot (default = 10)
# npseudo: a total number of pseudospots (default = 20000)

# alpha: loss weights of domain classifier to the source classifier (default = 0.6)
# alpha_lr: learning rate for the domain classifier (alpha_lr*0.001, default = 5)
# emb_dim: output size of dimensions for feature extractor (default = 64)

# batch_size: minibatch size for pseudospots and spatial data during the training (default = 512)
# n_iterations: iteration number for the adversarial learning (default = 3000)
# init_train_epoch: iteration number for the pre-training process (default = 10)

## Output
# sp_data_sub: spatial data (Seurat object) with predicted cell fraction in metadata (@meta.data)

## 1. Example for using virtual environment
sp_data_cellf <- pred_cellf_celldart(sp_data,sc_data,outdir='.',
                                     source.code.dir='.',
                                     sp_subset=TRUE,spot.cluster.name='seurat_clusters',
                                     spot.cluster.of.interest=NULL,
                                     env.select='virtual',python.install=F,
                                     python_path=NULL,env.name='spSeudoMap',
                                     gpu=TRUE,metadata_celltype='celltype',
                                     num_markers=10,mixture_mode='pseudotype',
                                     seed_num=0,
                                     mk_ratio_fix=T, mk_ratio=2,
                                     pseudo_num_genes=40,
                                     pseudo_frac_m=0.5,pseudo_frac_std=0.05, num_top_genes=20,
                                     nmix=10,npseudo=60000,alpha=1,alpha_lr=5,
                                     emb_dim=64,batch_size=512,n_iterations=3000,init_train_epoch=10)

## 2. Example for using conda environment
# For Windows OS, install conda environment first and then call it using pred_cellf_celldartconda.env.name='spatial'
# Conda environment should be already installed via Anaconda
sp_data_cellf <- pred_cellf_celldart(sp_data,sc_data,outdir='.',
                                     source.code.dir='.',
                                     sp_subset=TRUE,spot.cluster.name='seurat_clusters',
                                     spot.cluster.of.interest=NULL,
                                     env.select='conda',python.install=F,
                                     python_path=NULL,env.name='spSeudoMap',
                                     gpu=TRUE,metadata_celltype='celltype',
                                     num_markers=10,mixture_mode='pseudotype',
                                     seed_num=0,
                                     mk_ratio_fix=T, mk_ratio=2,
                                     pseudo_num_genes=40,
                                     pseudo_frac_m=0.5,pseudo_frac_std=0.05, num_top_genes=20,
                                     nmix=10,npseudo=60000,alpha=1,alpha_lr=5,
                                     emb_dim=64,batch_size=512,n_iterations=3000,init_train_epoch=10)

## Save seurat object with cell fraction
saveRDS(sp_data_cellf, file.path(output_folder_name, 'sp_data_cellf.rds'))



### 2. Visualization of spatial cell fraction
# Remove '_cellf' from the column names cell fraction metadata
cellf.data <- sp_data_cellf@meta.data
cellf.data.colname <- sapply(colnames(cellf.data), function(x){
  if (grepl('_cellf',x)){return(strsplit(x, split='_cellf')[[1]][1])}
  else {return(x)}
})
sp_data_cellf.mod <- sp_data_cellf
colnames(sp_data_cellf.mod@meta.data) <- cellf.data.colname

# Visualize the layer-specific excitatory neuons
cell_types <- c("L2.3.IT","L4","L5.IT","L5.PT","L6b","L6.CT","L6.IT")
p <- Seurat::SpatialFeaturePlot(sp_data_cellf.mod, features = cell_types, 
                                ncol = 4, alpha=0.6, combine = FALSE)
patchwork::wrap_plots(p, ncol = 8)
### Function to implement spSeudoMap in python
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
pred_cellf_spSeudoMap <- function(sp_data,sc_data,outdir='.',
                                  source.code.dir='.',
                                  sp_subset=FALSE,spot.cluster.name='seurat_clusters',
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
                                  emb_dim=64,batch_size=512,n_iterations=3000,init_train_epoch=10){
  
  # Suppress warnings
  defaultW <- getOption("warn") 
  options(warn = -1)
  
  if (python.install){
    reticulate::install_python(version = '3.7.12')
  }
  
  python_depend = c("scanpy==1.5.1","pandas==1.3.5","numpy==1.21.6",
                    "h5py==2.10.0", "scipy==1.7.3", "scikit-learn==1.0.2", "jupyter==1.0.0",
                    "keras==2.3.1", "tensorflow==1.14.0", "tensorflow-gpu==1.14.0")
  
  # Select between using reticulate virtual environment or conda environment ("virtual" or "conda")
  if (env.select=="virtual"){
    # Setting virtual environment with reticulate
    if (!(env.name %in% reticulate::virtualenv_list())){
      ## Python dependencies use python version 3.7.12
      if (is.null(python_path)){
        reticulate::virtualenv_create(envname = env.name, version = '3.7.12')
      } else if (python_path=="current") {
        reticulate::virtualenv_create(envname = env.name, python = NULL)
      } else {
        reticulate::virtualenv_create(envname = env.name, python = python_path)
      }

      # Create virtual env and install dependencies
      reticulate::virtualenv_install(env.name, packages = python_depend, ignore_installed=T)
    } 
      reticulate::use_virtualenv(env.name, required = T)
      
  } else if (env.select=="conda"){
    if (!(env.name %in% reticulate::conda_list()[['name']])){
      ## Python dependencies use python version 3.7
      if (is.null(python_path)){
        reticulate::conda_create(envname = env.name, python_version = '3.7.12')
      } else if (python_path=="current") {
        reticulate::conda_create(envname = env.name, python = NULL)
      } else {
        reticulate::conda_create(envname = env.name, python = python_path)
      }
      
      # Create conda env and install dependencies
      reticulate::conda_install(env.name, packages = python_depend, ignore_installed=T,
                                channel = c("bioconda", "conda-forge", "defaults"))
    }
    # Apply conda environment
    reticulate::use_condaenv(env.name, required = T)
  } else {
    stop("'env.select' should be either 'virtual' or 'conda'")
  }
  
  
  ## Import anndata
  scanpy_data <- reticulate::import('anndata', convert = FALSE)
  
  ## Import python function
  reticulate::source_python(file.path(source.code.dir, 'pred_cellf_spSeudoMap.py'))
  
  ## 1. Saving single-cell data in anndata format
  # Define count matrix
  sparse_mtx <- Seurat::GetAssayData(sc_data, slot = "counts", assay = "RNA")
  
  
  # Define obs and var (reference from sceasy library)
  obs <- sc_data@meta.data
  if (!metadata_celltype %in% colnames(obs)){
    stop("Column name for the cell annotation should be provided.")
  } else {
    obs <- obs[metadata_celltype]
    obs[[metadata_celltype]] <- factor(obs[[metadata_celltype]])
  }
  var <- data.frame(matrix(nrow=dim(sc_data)[1],ncol=0,
                           dimnames = list(rownames(sc_data),NULL)))
  var[['name']] <- rownames(var)
  
  adata_sc <- scanpy_data$AnnData(
    X = Matrix::t(sparse_mtx),
    obs = obs,
    var = var
  )
  
  ## 2. Subsetting spatial data and save in anndata format
  if (sp_subset){
    cluster_info <- sp_data[[spot.cluster.name]]
    Seurat::Idents(sp_data) <- spot.cluster.name
  }
  
  if (is.null(spot.cluster.of.interest)){
    sp_data_sub <- sp_data
  } else if (sum(spot.cluster.of.interest%in%levels(cluster_info))==length(spot.cluster.of.interest)){
    sp_data_sub <- subset(sp_data, idents=spot.cluster.of.interest)
  } else {
    stop("'spot.cluster.of.interest' should be among the levels of 'spot.cluster.name' provided")
  }
  
  # Define count matrix
  sparse_mtx <- Seurat::GetAssayData(sp_data_sub, slot = "counts", assay = "Spatial")
  
  # Define obs and var (reference from sceasy library)
  obs <- sp_data_sub@meta.data
  var <- data.frame(matrix(nrow=dim(sp_data_sub)[1],ncol=0,
                           dimnames = list(rownames(sp_data_sub),NULL)))
  var[['name']] <- rownames(var)
  
  adata_sp <- scanpy_data$AnnData(
    X = Matrix::t(sparse_mtx),
    obs = obs,
    var = var
  )
  
  # Assign the output directory for the models generated
  if (!file.exists(outdir)){
    dir.create(file.path(outdir, 'results'))
  }
  out_dir <- file.path(getwd(), outdir, 'results')
  
  # Save original directory and set python source code directory
  orig_dir <- getwd()
  file_path <- file.path(source.code.dir)
  setwd(file_path)
  
  # RUN spSeudoMap
  try({
    df <- pred_cellf_spSeudoMap(adata_sp=adata_sp, adata_sc=adata_sc, count_from_raw = F, 
                                gpu=gpu, celltype=metadata_celltype, num_markers=num_markers,
                                mixture_mode=mixture_mode, seed_num=seed_num, 
                                mk_ratio_fix=mk_ratio_fix, mk_ratio=mk_ratio, 
                                pseudo_num_genes=pseudo_num_genes,
                                pseudo_frac_m=pseudo_frac_m, pseudo_frac_std=pseudo_frac_std,
                                num_top_genes=num_top_genes,
                                nmix=nmix, npseudo=npseudo, alpha=alpha, alpha_lr=alpha_lr, 
                                emb_dim=emb_dim, batch_size=batch_size, 
                                n_iterations=n_iterations, init_train_epoch=init_train_epoch, 
                                outdir=out_dir, 
                                return_format='dataframe')
    
    # Saving cell fraction data into the metadata of spatial Seurat object
    sp_data_sub <- Seurat::AddMetaData(sp_data_sub, df)
  })
  
  # Set to original working directory
  setwd(orig_dir)
  
  options(warn = defaultW)
  
  return(sp_data_sub)
  
}

# spSeudoMap  
spSeudoMap: cell type mapping of spatial transcriptomics using unmatched single-cell RNA-seq data  
![Fig  1](https://user-images.githubusercontent.com/61150422/167383248-417cf084-4ecc-4659-979c-54a64c45ec4f.jpg)

## Optimal parameter choices  
  Number of marker genes per cluster: 40 (>20)  
  m/k ratio = 2 (> 1)  
  pseudo_frac_m = average fraction of negative non-sorted population (literature evidence or cell sorting experiment)  
  pseudo_frac_std = 0.1 (> 0.05)  
  Number of pseudospots = 5 to 10 times the number of real spots (20,000~40,000 per Visium slide)  
  Number of sampled cells in a pseudospot (virtual mixture of single-cell data) = 8 (brain), 10 (breast cancer)  
  Iteration number = 3,000  
  Mini-batch size = 512  
  Loss weights between source and domain classifier (alpha) = 0.6  
  Learning rate = 0.001 * alpha_lr = 0.005  

## Code Example  
  Python example: spSeudoMap_example.ipynb  
  R example: Refer to /vignettes/introduction.Rmd  

## Data Example  
  spatial data: 'V1_Adult_Mouse_Brain_Coronal_Section_1' from 10X Genomics Repository  

## Python function for spSeudoMap (pred_cellf_spSeudoMap)  
### Install conda environment and add kernel (jupyter)  
    conda create -n spSeudoMap python=3.8
    conda activate spSeudoMap  
    pip install git+https://github.com/bsungwoo/spSeudoMap.git  
    python -m ipykernel install --user --name spSeudoMap --display-name spSeudoMap  

### Function and main parameters  
``` Plain Text
from spSeudoMap.pred_cellf_spSeudoMap import pred_cellf_spSeudoMap  
adata_sp = pred_cellf_spSeudoMap(adata_sp=None, adata_sc=None, count_from_raw=False,   
                                 gpu=True, celltype='cluster', num_markers=40,  
                                 mixture_mode='pseudotype', seed_num=0,  
                                 mk_ratio_fix=False, mk_ratio=2, pseudo_num_genes=40,  
                                 pseudo_frac_m=0.5, pseudo_frac_std=0.1, num_top_genes=20,  
                                 nmix=10, npseudo=20000, alpha=0.6, alpha_lr=5, emb_dim=64, 
                                 batch_size=512, n_iterations=3000, init_train_epoch=10, 
                                 outdir='./output', return_format='anndata')  
```
  **(1) adata_sp:** spatial data (AnnData object) with raw count matrix to be used in predicting cell fraction (default: None)    
  **(2) adata_sc:** single-cell data (AnnData object) with raw count matrix to be used in making pseudospots (default: None)  
  **(3) count_from_raw:** whether to extract count matrix frow .raw of AnnData (default: False)  
  -> non-normalized raw count matrix should be contained in the AnnData .raw file  
  -> if False, then utilize the count matrices saved in adata_sp and adata_sc directly  
  **(4) gpu:** check whether to use gpu (True) or not (False) (default = True)  
  **(5) celltype:** column name for single-cell annotation data in .obs (default: 'celltype')  
  **(6) num_markers:** number of selected marker genes in each celltype (default = 40)   
  **(7) mixture_mode:** mode of the pseudospot generation ('default': same as CellDART, 'pseudotype': assuming exclusive cell type in spatial data)  
  **(8) seed_num:** seed to be used in random sampling (default = 0)  
  **(9) mk_ratio_fix:** whether to fix the mk_ratio when selecting the pseudotype markers (default = True)  
  **(10) mk_ratio:** ratio of number of single-cell markers to virtual pseudotype markers (default = 2)  
  **(11) pseudo_frac_m:** average of presumed fraction of the pseudotype (cell types exclusively present in spatial data) across all spatial spots (default = 0.5)  
  -> determined by cell sorting study or by literature based evidence that enables speculation of average negative non-sorted cell fraction in the tissue  
  **(12) pseudo_frac_std:** standard deviation of the distribution of presumed pseudotype fraction across all spatial spots (default = 0.1)  
  **(13) num_top_genes:** number of top genes having highest log fold change between spatial and single-cell normalized pseudobulk counts (spatial/single-cell) (default = 20)  
  **(14) nmix:** sampling number of cells in pseudospot (default = 10)  
  **(15) npseudo:** a total number of pseudospots (default = 20,000); approximately 5~10 times the number of pseudospots  

### Training parameters  
  **(1) alpha:** loss weights of the domain classifier to the source classifier (default = 0.6)  
  **(2) alpha_lr:** learning rate for the domain classifier (alpha_lr*0.001, default = 5)  
  **(3) emb_dim:** output size of dimensions for feature extractor (default = 64)  
  **(4) batch_size:** minibatch size for pseudospots and spatial data during the training (default = 512)  
  **(5) n_iterations:** iteration number for the adversarial learning (default = 3,000)  
  **(6) init_train_epoch:** iteration number for the pre-training process (default = 10)  
  **(7) outdir:** the directory to save output files (models and results)  
  **(8) return_format:** whether to return spatial AnnData file with predicted cell fraction in .obs (default: False)  

## R wrap function for spSeudoMap (spSeudoMap::pred_cellf_spSeudoMap)
    devtools::install_github("bsungwoo/spSeudoMap", build_vignettes = T, force = T)  
    library(spSeudoMap)  
    help(pred_cellf_spSeudoMap) # Explanation for the parameters and short examples  
    browseVignettes("spSeudoMap")  # Browse for the vignettes (/vignettes/introduction.Rmd)  

  ### Installation of virtual or conda environment
  Linux distributions: The environment will be automatically installed by running the function   
  Windows: Install conda environment first and then run the function with env.select = 'conda' and python.install=F  

  ## Function and additional parameters
  ```Plain Text
  Using conda environment (environment will be automatically installed in Linux distributions)
  If using Windows, then install conda environment first and then run the function below with python.install = F
  pseudo_frac_m <- 0.5 (average presumed pseudotype fraction in the tissue)
  npseudo <- dim(sp_data)[2]*5
  sp_data_cellf <- pred_cellf_spSeudoMap(sp_data, sc_data,
                                         outdir=output_folder_name,
                                         sp_subset=F, spot.cluster.name='seurat_clusters',
                                         spot.cluster.of.interest=NULL,
                                         env.select='conda', python.install=T,
                                         python_path=NULL, env.name='spSeudoMap',
                                         gpu=TRUE, metadata_celltype='annotation_1',
                                         num_markers=40, mixture_mode='pseudotype',
                                         seed_num=0,
                                         mk_ratio_fix=T, mk_ratio=4,
                                         pseudo_frac_m=pseudo_frac_m, pseudo_frac_std=0.1,
                                         nmix=8, npseudo=nspeudo, alpha=0.6, alpha_lr=5,
                                         emb_dim=64, batch_size=512, n_iterations=3000, init_train_epoch=10)
 ```
 ``` Plain Text
  Using virtual environment (environment will be automatically installed in Linux distributions)
  Not recommended for Windows
  pseudo_frac_m <- 0.5 (average presumed pseudotype fraction in the tissue)
  npseudo <- dim(sp_data)[2]*5
  sp_data_cellf <- pred_cellf_spSeudoMap(sp_data, sc_data,
                                         outdir=output_folder_name,
                                         sp_subset=F, spot.cluster.name='seurat_clusters',
                                         spot.cluster.of.interest=NULL,
                                         env.select='virtual', python.install=T,
                                         python_path=NULL, env.name='spSeudoMap',
                                         gpu=TRUE, metadata_celltype='annotation_1',
                                         num_markers=40, mixture_mode='pseudotype',
                                         seed_num=0,
                                         mk_ratio_fix=T, mk_ratio=4,
                                         pseudo_frac_m=pseudo_frac_m, pseudo_frac_std=0.1,
                                         nmix=8, npseudo=nspeudo, alpha=0.6, alpha_lr=5,
                                         emb_dim=64, batch_size=512, n_iterations=3000, init_train_epoch=10)
  ```  
  **(1) outdir:** the directory to save output files (models and results) (default = '.')  
  **(2) sp_subset:** whether to subset spatial data and calculate for specific spot cluster (default = FALSE)  
  **(3) spot.cluster.name:** group name of the cluster used for subsetting spatial data (default = 'seurat_clusters')  
  **(4) spot.cluster.of.interest:** name of each spot clusters to be used (default = NULL)  
  **(5) env.select:** select between using reticulate virtual environment or conda environment (default = 'conda')  
  -> either of the selection will search the already installed environment  
  -> if environment is not found, then it will automatically install the new environment  
  **(6) python.install:** whether to automatically install python version 3.7.12 (default = F)  
  -> For Windows, set python.install = F  
  **(7) python_path:** path for the python 3.7.12 (default = NULL)  
  **(8) env.name:** name of the virtual or conda environment to use for the analysis (default = 'spSeudoMap')  
  **(9) metadata_celltype:** column name for single-cell annotation data in metadata (default = 'celltype')  

### Potential errors when installing new conda environment by R wrap function (env.select='conda')
  **GLIBCXX_3.4.26 should be already installed**  
  Check if it is installed: strings /usr/lib/x86_64-linux-gnu/libstdc++.so.6 | grep GLIBCXX  
  reference: https://stackoverflow.com/questions/63190229/glibcxx-3-4-26-not-found-running-cross-complied-program-on-beaglebone  

    sudo apt install wget gcc-8 unzip libssl1.0.0 software-properties-common  
    sudo add-apt-repository ppa:ubuntu-toolchain-r/test  
    sudo apt-get install --only-upgrade libstdc++6  




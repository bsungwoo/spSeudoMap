import scanpy as sc
import numpy as np
import pandas as pd
from keras.utils import to_categorical

    
### Make random mix for pseudospot
# Return subsetted array for single-cell data (Xs), composite gene expression array for pseudospot (Xs_new) and cell fraction (ys_new)

## Input
# adata_sc: original single-cell count matrix (total count normalized, but not log-transformed) that contains all the features
# ys: array of numbers representing the cell type categories that each cell in single-cell data is included
# genes_of_interest: genes of interest to generate pseudospot, the cell mixture from single-cell data
# nmix: the number of cells sampled from single-cell data when making a pseudospot (default = 10)
# n_samples: number of pseudospots to generate (default: 10000)
# seed: set seed number for the random value generation (default: seed=0)

## Output
# Xs: subetted numpy ndarray for total count normalized (but not log-transformed) expression of genes_of_interest in single-cell data
# Xs_new: composite gene expression array across "n_samples" number of pseudospots
# ys_new: fraction of cell types across the "n_samples" number of pseudospots
# example of the shape of the output when the number of pseudospots: 20000, number of genes in 'genes_of_interest': 771, number of cell types: 33
# Xs_new.shape = (20000, 771), ys.shape = (20000, 33)
def random_mix(adata_sc, ys, genes_of_interest, nmix=10, n_samples=10000, seed=0):
    # Define empty lists
    Xs_new, ys_new = [], []
    
    # Convert single-cell anndata to array
    adata_sc_sub = adata_sc[:,genes_of_interest].copy()
    if isinstance(adata_sc.X, np.ndarray):
        Xs = adata_sc_sub.X
    else:
        Xs = adata_sc_sub.X.toarray()
        
    ys_ = to_categorical(ys)
    
    rstate = np.random.RandomState(seed)
    fraction_all = rstate.rand(n_samples, nmix)
    randindex_all = rstate.randint(len(Xs), size=(n_samples,nmix))
    
    for i in range(n_samples):
        # fraction: random fraction across the "nmix" number of sampled cells
        fraction = fraction_all[i]
        fraction = fraction/np.sum(fraction)
        fraction = np.reshape(fraction, (nmix,1))
        
        # Random selection of the single cell data by the index
        randindex = randindex_all[i]
        ymix = ys_[randindex]
        # Calculate the fraction of cell types in the cell mixture
        yy = np.sum(ymix*fraction, axis=0)
        # Calculate weighted gene expression of the cell mixture
        XX = np.asarray(Xs[randindex])*fraction
        XX_ = np.sum(XX, axis=0)
        
        # Add cell type fraction & composite gene expression in the list
        ys_new.append(yy)
        Xs_new.append(XX_)

    Xs_new = np.asarray(Xs_new)
    ys_new = np.asarray(ys_new)

    return Xs, Xs_new, ys_new



### Calculate pseudobulk fold change and extract top gene lists
# Return pandas dataframe for pseuobulk fold change between spatial and single-cell pseudobulk

## Input
# adata_sc: original single-cell AnnData that contains all features (total count normalized, but not log-transformed count matrix)
# adata_sp: original spatial AnnData that contains all features (total count normalized, but not log-transformed count matrix)
# gene_list: genes of interest to calculate pseudobulk fold change
# -> if None, then gene_list is automatically set to the intersecting genes between single-cell and spatial data

# raw_count_layer_name: name of the layer that raw non-normalized count data is saved (default: "raw")
# sp_ratio_thres: cutoff for the pseudobulk count to the total pseudobulk count across all cells in spatial data (default: 0)
# sc_ratio_thres: cutoff for the pseudobulk count to the total pseudobulk count across all spots in single-cell data (default: 0)

## Output
# df_foldchange: pandas dataframe containing pseudobulk count, ratio of pseudobulk count to total count, 
# and pandas dataframe containing log fold change between spatial and single-cell normalized pseudobulk counts (spatial/single-cell)
def pseudobulk_foldchange(adata_sc, adata_sp, gene_list=None, raw_count_layer_name="raw",
                            sp_ratio_thres=0, sc_ratio_thres=0):
    # Find intersecting genes btw single-cell & spatial data & df_HVG_sort, excluding marker genes
    if gene_list is None:
        gene_intersect = [i for i in adata_sp.var.index if (i in adata_sc.var.index)]
    else:
        gene_intersect = [i for i in gene_list if ((i in adata_sp.var.index) and (i in adata_sc.var.index))]

    if isinstance(adata_sp.layers[raw_count_layer_name], np.ndarray):
        Xsp = adata_sp[:,gene_intersect].layers[raw_count_layer_name]
    else:
        Xsp = adata_sp[:,gene_intersect].layers[raw_count_layer_name].toarray()
                
    if isinstance(adata_sc.layers[raw_count_layer_name], np.ndarray):
        Xsc = adata_sc[:,gene_intersect].layers[raw_count_layer_name]
    else:
        Xsc = adata_sc[:,gene_intersect].layers[raw_count_layer_name].toarray()
            
    # Calculate pseudobulk expression
    sp_pseudobulk = np.sum(Xsp, axis=0)
    sc_pseudobulk = np.sum(Xsc, axis=0)
    df_foldchange = pd.DataFrame([sp_pseudobulk, sc_pseudobulk, 
                                    sp_pseudobulk/np.sum(sp_pseudobulk), sc_pseudobulk/np.sum(sc_pseudobulk)]).T
    df_foldchange.columns = ['sp_count','sc_count','sp_ratio','sc_ratio']
    df_foldchange['logfc'] = np.log2(df_foldchange['sp_ratio']) - np.log2(df_foldchange['sc_ratio'])
    df_foldchange.index = gene_intersect

    df_foldchange = df_foldchange[(df_foldchange['sp_ratio']>sp_ratio_thres) & (df_foldchange['sc_ratio']>sc_ratio_thres)]
    df_foldchange.sort_values(by=["logfc"], ascending=False, inplace=True)
    
    return df_foldchange

    
   
### Make random mix for pseudospot by adding pseudotype expression profiles based on pseudobulk transcriptomes
# Return subsetted array for single-cell data (Xsc), composite gene expression array for pseudospots containing pseudotype expression profiles (Xs_new), 
# and cell fraction (ys_new)

## Input: 
# adata_sc: original single cell data that contains all the features
# adata_sp: original spatial data that contains all the features
# ys: array of numbers representing the cell type categories that each cell in single-cell data is included
# genes_of_interest: genes of interest to generate pseudospot mixture
# marker_genes: marker genes for the cell types comprising single-cell data
# -> also the genes that expression in pseudotype is set to 0

# nmix: the number of cells sampled from single-cell data when making a pseudospot (default = 10)
# n_samples: number of pseudospots to generate (default: 10000)

# pseudo_frac_m: average of presumed fraction of the pseudotype (cell types exclusively present in spatial data) across all spatial spots (default = 0.5)
# -> determined by cell sorting study or by literature based evidence that enables speculation of average fraction in the tissue
# pseudo_frac_std: standard deviation of the distribution of presumed pseudotype fraction across all spatial spots (default = 0.1)

# ** Parameters for generating module score for the top genes provided **
# num_top_genes: number of top genes having highest log fold change between spatial and single-cell normalized pseudobulk counts (spatial/single-cell)
# -> the top genes are used for generating module score to predict pseudotype fraction in spatial spots
# n_bin: number of levels to divide gene expression (default: n_bin=25)
# n_ctrl: number of control genes to select when constructing control gene set (default: n_ctrl=50)
# seed: set seed number during the module score generation (default: seed=0)

# df_foldchange: pandas dataframe containing log fold change between spatial and single-cell normalized pseudobulk counts (spatial/single-cell)

## Output
# Xsc: subetted numpy ndarray for total count normalized (but not log-transformed) expression of genes_of_interest in single-cell data
# Xs_new: composite gene expression array containing pseudotype expression profiles across "n_samples" number of pseudospots
# ys_new: fraction of cell types across the "n_samples" number of pseudospots, pseudotype is written as "Others"
# example of the shape of the output when the number of pseudospots: 20000, number of genes in 'genes_of_interest': 771, number of cell types: 33
# Xs_new.shape = (20000, 771), ys.shape = (20000, 34)
def random_mix_pseudobulk(adata_sc, adata_sp, ys, genes_of_interest, marker_genes,
                          nmix=10, n_samples=10000, pseudo_frac_m = 0.5, pseudo_frac_std = 0.1,
                          num_top_genes=20, n_bin=25, n_ctrl=50, seed=0, df_foldchange=None):
    # nclss: number of cell types in the single-cell data
    nclss=len(set(ys))
    Xs_new, ys_new = [], []

    # Check feasibility of df_foldchange
    if not isinstance(df_foldchange, pd.DataFrame) or (df_foldchange is None):
        df_foldchange = pseudobulk_foldchange(adata_sc, adata_sp, raw_count_layer_name="raw")

    # Convert spatial anndata to array
    adata_sp_sub = adata_sp[:,genes_of_interest].copy()
    if isinstance(adata_sp.X, np.ndarray):
        Xsp = adata_sp_sub.X
    else:
        Xsp = adata_sp_sub.X.toarray()
    
    # Convert single-cell anndata to array
    adata_sc_sub = adata_sc[:,genes_of_interest].copy()
    if isinstance(adata_sc.X, np.ndarray):
        Xsc = adata_sc_sub.X
    else:
        Xsc = adata_sc_sub.X.toarray()
        
    ys_ = to_categorical(ys)

    # Assign the index number for the columns
    adata_sc_sub.var['num'] = range(len(genes_of_interest))

    # Finding index for the marker genes and virtual pseudotype markers
    genes_marker_gs = [val for val in genes_of_interest if val in marker_genes]
    genes_pseudo_gs = [val for val in genes_of_interest if val not in marker_genes]
    index_rev = adata_sc_sub[:,genes_marker_gs].var['num']
    
    # Assign random state and extract random number
    rstate = np.random.RandomState(seed)
    fraction_all = rstate.rand(n_samples,nmix)
    randindex_all = rstate.randint(len(Xsc),size=(n_samples,nmix))

    ## Predict the pseudotype fraction from the module score of single-cell mixture marker genes
    pseudo_markers = df_foldchange[df_foldchange['logfc']>0].sort_values(by=["logfc"], ascending=False)
    pseudo_markers = pseudo_markers[pseudo_markers.index.isin(genes_pseudo_gs)]
    pseudo_markers = pseudo_markers[:num_top_genes].index.tolist()
    
    # Calculate module score for single-cell mixture marker marker genes
    sc.pp.log1p(adata_sp)
    sc.tl.score_genes(adata_sp, gene_list=pseudo_markers, ctrl_size=n_ctrl, n_bins=n_bin, 
                        score_name='pseudotype', random_state=seed, use_raw=False)
    pseudo_module = adata_sp.obs['pseudotype']

    print("Genes used for calculation of pseudotype fraction")
    print(pseudo_markers)
    
    # normalize the distribution of module score for single-cell mixture markers
    frac_pseudo = (pseudo_module-np.mean(pseudo_module))/(np.std(pseudo_module))
    # Scale the distribution such that mean value matches with pseudo_frac_m and standard deviation to pseudo_frac_std
    frac_pseudo = pseudo_frac_m + frac_pseudo*pseudo_frac_std
    # Fix the lowest pseudotype fraction value as 0 and highest as 1
    frac_pseudo[frac_pseudo<0] = 0
    frac_pseudo[frac_pseudo>1] = 1

    # Assign random integers for the spot selection (expression & presumed fraction)
    randspatial_sel = rstate.randint(len(Xsp), size=(n_samples,1))

    for i in range(n_samples):
        # fraction: random fraction across the "nmix" number of sampled cells, including pseudotype fraction
        fraction = fraction_all[i]
        fraction_ = fraction*(1-frac_pseudo[randspatial_sel[i][0]])/np.sum(fraction)
        fraction_ = np.append(fraction_, frac_pseudo[randspatial_sel[i][0]])
        fraction_ = np.reshape(fraction_, (nmix+1,1))
        
        # Random selection of the single-cell data by the index
        randindex = randindex_all[i]
        ymix = ys_[randindex]
        
        # Add the pseudotype and its fraction
        ymix_ = np.hstack((ymix, np.zeros((nmix,1))))
        ymix_ = np.vstack((ymix_, np.reshape([0]*nclss+[1],(1,nclss+1))))
        # Calculate the fraction of cell types in the cell mixture
        yy = np.sum((ymix_*fraction_), axis=0)

        # Assign the expression of pseudotype from spatial data
        pseudo_exp = Xsp[randspatial_sel[i]]
        # Put 0 values to the marker gene index (for single-cell) -> only gives weight to virtal pseudotype marker expression
        np.put(pseudo_exp, index_rev, np.zeros((1,len(marker_genes))))
        
        # Add the synthetic gene expression profiles of pseudotype
        Xmix = np.vstack((Xsc[randindex], pseudo_exp))
        # Calculate weighted gene expression of the cell mixture
        XX = Xmix*fraction_
        XX_ = np.sum(XX, axis=0)
        
        # Add fraction of cell type & composite gene expression of modified pseudospot in the list
        ys_new.append(yy)
        Xs_new.append(XX_)
        
    Xs_new = np.asarray(Xs_new)
    ys_new = np.asarray(ys_new)

    return Xsc, Xs_new, ys_new

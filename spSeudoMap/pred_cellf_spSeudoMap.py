### Function to implement spSeudoMap in python
# cell type mapping of spatial transcriptomics using unmatched single-cell RNA-seq data

## Input
# adata_sp: spatial data (AnnData object) to be used in predicting cell fraction (default: None)
# -> count matrix should be non-normalied raw data
# adata_sc: single-cell data (AnnData object) to be used in making pseudospots (default: None)
# -> count matrix should be non-normalied raw data
# count_from_raw: whether to extract count matrix frow .raw of AnnData (default: False)
# -> non-normalized count matrix should be contained in the AnnData .raw file
# -> if False, then utilize the count matrices saved in adata_sp and adata_sc directly

# gpu: check whether to use gpu (True) or not (False) (default = True)
# celltype: column name for single-cell annotation data in metadata(.obs) (default: 'celltype')
# num_markers: number of selected marker genes in each cell-type (default = 40)

# mixture_mode: mode of the pseudospot generation 
# -> 'default': when cell types are similar between single-cell/spatial data (Identical to CellDART)
# -> 'pseudotype': when there is a mismatch of cell types between single-cell/spatial data and cell types exclusively present in spatial data is considered when generating pseudospots

# seed_num: seed to be used in random sampling (default = 0)

# mk_ratio_fix: whether to fix the mk_ratio when selecting the pseudotype markers (default = True)
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

# outdir: the directory to save output files (models and results)
# return_format: which data file to return (default: 'anndata')
# -> if 'dataframe', return dataframe for predicted cell fraction
# -> if 'anndata', return spatial data (AnnData object) with predicted cell fraction in .obs
# -> if others, return None

## Output
# spatial_raw: spatial data (AnnData object) with predicted cell fraction in metadata (.obs)
# df: dataframe for predicted cell fraction across all spatial spots
def pred_cellf_spSeudoMap(adata_sp=None, adata_sc=None, count_from_raw=False, 
                            gpu=True, celltype='celltype', num_markers=40,
                            mixture_mode='pseudotype', seed_num=0, 
                            mk_ratio_fix=True, mk_ratio=2, pseudo_num_genes=40, 
                            pseudo_frac_m=0.5, pseudo_frac_std=0.1, num_top_genes=20,
                            nmix=10, npseudo=20000, alpha=0.6, alpha_lr=5, emb_dim=64, 
                            batch_size=512, n_iterations=3000, init_train_epoch=10, 
                            outdir='./output', return_format='anndata'):

    import os
    if gpu:
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"]= "0" # Use only gpu-0
        print('GPU is available and will be used')
    else:
        os.environ['CUDA_VISIBLE_DEVICES'] = "-1" # Use CPU
        print('CPU will be used')
    
    from warnings import simplefilter
    simplefilter(action='ignore', category=Warning)

    import scanpy as sc
    import pandas as pd
    import numpy as np

    from spSeudoMap import da_cellfraction
    from spSeudoMap import utils_mod
    
    ## Change float variables into integer (during conversion from R to python)
    num_markers, seed_num, pseudo_num_genes, num_top_genes, \
    nmix, npseudo, batch_size, emb_dim, n_iterations, init_train_epoch = \
        int(num_markers), int(seed_num), \
        int(pseudo_num_genes), int(num_top_genes), int(nmix), int(npseudo), int(batch_size), int(emb_dim), \
        int(n_iterations), int(init_train_epoch)
    
    ## Create directory if it does not exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    ## Load and preprocess spatial dataset
    if count_from_raw: spatial_all = adata_sp.raw.to_adata()
    else: spatial_all = adata_sp.copy()
    # Make variable names unique for spatial data
    spatial_all.var_names_make_unique()
    print('Shape of the provided spatial data is',spatial_all.shape)
     
    # Generate spatial data with raw counts
    if mixture_mode == 'pseudotype': 
        spatial_all.layers['raw'] = spatial_all.X.copy()
            
    # Total count normalize
    sc.pp.normalize_total(spatial_all, target_sum=1e4, inplace=True)

    ## Load and preprocess single-cell dataset
    if count_from_raw: single_all = adata_sc.raw.to_adata()
    else: single_all = adata_sc.copy()
    # Make variable names unique for spatial data
    single_all.var_names_make_unique()

    # Check if the column for the cell type is included in metadata(.obs)
    if celltype not in list(single_all.obs):
        raise ValueError('Column for cell type is not found')

    print('Shape of the provided single-cell data is',single_all.shape)
        
    # Generate single-cell data with raw counts
    if mixture_mode == 'pseudotype':
        single_all.layers["raw"] = single_all.X.copy()
            
    # Total normalize
    sc.pp.normalize_total(single_all, target_sum=1e4, inplace=True)
    
    # save the normalized data in raw
    single_all.raw = single_all
    
    # log-transform the count matrix
    sc.pp.log1p(single_all)
    
    # Find marker genes for single cell data
    single_all.obs[celltype] = single_all.obs[celltype].astype('category', copy=False)
    
    # Check if there is any cluster that contains only one cell
    df_test = single_all.obs.groupby(celltype)[celltype].count().sort_values()
    if df_test[0] == 1:
        raise ValueError("The number of cell is 1 in cluster: '"+str(df_test.index[0])+"'")

    sc.tl.rank_genes_groups(single_all, celltype, method='wilcoxon')
    genelists = single_all.uns['rank_genes_groups']['names']
    df_genelists = pd.DataFrame.from_records(genelists)
    
    # Combining top marker genes representing each cell type
    res_genes = []
    for column in df_genelists.head(num_markers): 
        res_genes.extend(df_genelists.head(num_markers)[column].tolist())
    res_genes_ = list(set(res_genes))

    if mixture_mode == 'default':
        # Find intersecting genes
        inter_genes_comb = [val for val in res_genes_ if val in spatial_all.var.index]
        print('Total number of genes for training: '+str(len(inter_genes_comb)))
    
    elif mixture_mode == 'pseudotype':
        # Find intersecting genes
        inter_genes = [val for val in res_genes_ if val in spatial_all.var.index]
        print('Total number of cell type marker genes is '+str(len(inter_genes)))
        
        # Determine the number of pseudogenes when the m/k ratio is fixed
        if mk_ratio_fix:
            pseudo_num_genes = int(len(inter_genes)/mk_ratio)
        else: 
            pseudo_num_genes = pseudo_num_genes
        
        ## Extract virtual pseudotype markers from pseudobulk transcriptomes of single-cell and spatial data
        # Find the fold change list for the pseudobulk analysis
        df_pseudo = utils_mod.pseudobulk_foldchange(single_all, spatial_all, raw_count_layer_name="raw")
        print("Dataframe summarizing counts of pseudobulk transcriptomes and log fold change")
        print(df_pseudo)

        # Extract genes showing high log fold change in spatial > single-cell pseudobulk (descending order)
        # Find the exclusive genes having the positive log fold change between spatial > single-cell pseudobulk
        inter_genes_pseudo = df_pseudo[(df_pseudo["logfc"]>0) & (~df_pseudo.index.isin(inter_genes))]
        # Check if the number of the total positive genes is smaller than the given pseudotype marker numbers
        if len(inter_genes_pseudo) < pseudo_num_genes:
            pseudo_num_genes = len(inter_genes_pseudo)
            print("""Number of pseudogenes exceed the number of genes discovered by pseudobulk: 
            Replacing with maximum possible number""")
        inter_genes_pseudo = inter_genes_pseudo[:pseudo_num_genes].index.tolist()
        inter_genes_comb = inter_genes + inter_genes_pseudo
    
        print('Total number of genes for training: '+str(len(inter_genes_comb)))
    
    else:
        raise ValueError('Input for mixture mode should be either default or pseudotype')
    
    
    # Generation of an array representing cell type number
    df_sc = single_all.obs
    lab_sc_sub = df_sc[celltype]
    sc_sub_dict = dict(zip(range(len(set(lab_sc_sub))), set(lab_sc_sub)))
    sc_sub_dict2 = dict((y,x) for x,y in sc_sub_dict.items())
    lab_sc_num = [sc_sub_dict2[ii] for ii in lab_sc_sub]
    # Make an array of numbers representing the cell type categories that each cell in single-cell data is included
    lab_sc_num = np.asarray(lab_sc_num, dtype='int')
    
    # Call original total normalized count (not log-normalized count)
    adata_final = single_all.raw.to_adata()

    # Generate pseudospot: random mixture of cells
    if mixture_mode == 'default':
        # When the cell types comprising the single-cell and spatial data is similar
        mat_sc, sc_mix, lab_mix = utils_mod.random_mix(adata_final, lab_sc_num, inter_genes_comb, nmix=nmix, 
                                                        n_samples=npseudo, seed=seed_num)

    if mixture_mode == 'pseudotype':
        # When there is a mismatch in cell types comprising the single-cell and spatial data
        # pseudotype which explains cell types exclusively present in the spatial data is defined
        mat_sc, sc_mix, lab_mix = utils_mod.random_mix_pseudobulk(adata_final, spatial_all, lab_sc_num, inter_genes_comb,
                                                                    inter_genes, nmix=nmix, n_samples=npseudo, 
                                                                    pseudo_frac_m=pseudo_frac_m, pseudo_frac_std=pseudo_frac_std,
                                                                    num_top_genes=num_top_genes, seed=seed_num,
                                                                    df_foldchange=df_pseudo)
    
    # Raw file for merged spatial data
    spatial_raw = spatial_all

    # Generate count matrix for spatial data (mat_sp)
    spatial_all = spatial_all[:,inter_genes_comb].copy()
    if isinstance(spatial_all.X, np.ndarray):
        mat_sp = spatial_all.X
    else: 
        mat_sp = spatial_all.X.toarray()
    

    # Log-normalize and scale the data 
    def log_minmaxscale(arr):
        arrd = len(arr)
        arr = np.log1p(arr)
        e = 1e-8 # modified by adding e
        return (arr-np.reshape(np.min(arr,axis=1),(arrd,1)))/np.reshape((np.max(arr,axis=1)-np.min(arr,axis=1))+e,(arrd,1))

    sc_mix_s = log_minmaxscale(sc_mix)
    mat_sp_s = log_minmaxscale(mat_sp)
    mat_sc_s = log_minmaxscale(mat_sc)

    print('Size of spatial, single-cell, pseudospot, and cell fraction data:',
            mat_sp_s.shape, mat_sc_s.shape, sc_mix_s.shape, lab_mix.shape)
        

    # Train the domain adaptation model
    embs, clssmodel = da_cellfraction.train(sc_mix_s, lab_mix, mat_sp_s, 
                                            emb_dim = emb_dim, batch_size = batch_size,
                                            n_iterations = n_iterations,
                                            enable_dann = True,
                                            alpha = alpha, alpha_lr = alpha_lr,
                                            initial_train = True,
                                            initial_train_epochs = init_train_epoch)
    # Prediction of cell fraction in each spot
    pred_sp = pd.DataFrame(clssmodel.predict(mat_sp_s))
    pred_sp.index = spatial_all.obs.index

    # Make directory for the model save
    if not os.path.exists(os.path.join(outdir,'model')):
        os.makedirs(os.path.join(outdir,'model'))

    # Save the cell fraction in metadata(.obs) of spatial_raw file
    if mixture_mode == 'default':
        for visnum in range(len(sc_sub_dict)):
            spatial_raw.obs[str(sc_sub_dict[visnum])+'_cellf'] = pred_sp.iloc[pred_sp.index.isin(spatial_raw.obs.index),visnum]
                
    if mixture_mode == 'pseudotype':
        for visnum in range(len(sc_sub_dict)+1):
            if visnum == len(sc_sub_dict):
                spatial_raw.obs['Others_cellf'] = pred_sp.iloc[pred_sp.index.isin(spatial_raw.obs.index),visnum]
            else:
                spatial_raw.obs[str(sc_sub_dict[visnum])+'_cellf'] = pred_sp.iloc[pred_sp.index.isin(spatial_raw.obs.index),visnum]
        
    # Save cell fraction data
    df = spatial_raw.obs.filter(regex='_cellf', axis=1)
    df.to_csv(os.path.join(outdir,'cellfraction.csv'), header=True, index=True)
    print('Cell fraction data was saved')

    # Save spatial anndata
    spatial_raw.write_h5ad(os.path.join(outdir,'model/sp_data.h5ad'))
    print('Spatial anndata was saved')

    # Save model files
    embs.save_weights(os.path.join(outdir,'model/embedder.h5'))
    clssmodel.save_weights(os.path.join(outdir,'model/classifier.h5'))

    print('Model and python data files were saved')

    if return_format=='anndata': return spatial_raw
    elif return_format=='dataframe': return df
    else: return None

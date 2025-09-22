## load packages
import scanpy as sc
import scvi
from rich import print
from scib_metrics.benchmark import Benchmarker
from scvi.model.utils import mde
import numpy as np
import pymde
import pandas as pd

## set working directory
os.chdir("03_subset-bySample")

adata = sc.read_h5ad('../02_scVI_Integration-bySample/adata-rawCount-meta.h5ad')
#vector_rm = ['_scvi_batch','_scvi_labels','leiden_0.1','leiden_0.2','leiden_0.3','leiden_0.4','leiden_0.5','leiden_0.6','celltype_unknown']
#adata.obs.drop(columns=vector_rm, inplace=True)

filtered_celltypes = ['Myeloid cell', 'T cell', 'B cell', 'Fibroblast', 'Plasma cell', 'Mast cell', 'Endothelial cell', 'Glial cell']
for celltype in filtered_celltypes:
    flaginfo1 = '>>>>>> ' + 'start to process ' + celltype
    print(flaginfo1)
    folder_path = celltype
    if not os.path.exists(celltype):
        os.mkdir(celltype)
    else:
        print("dirctory exist")
    os.chdir(celltype)
    selected_cells = adata[adata.obs['celltype'] == celltype, :]
    ad_sub = sc.AnnData(selected_cells.X, obs=selected_cells.obs, var=selected_cells.var)
    ad_sub.layers["counts"] = ad_sub.X.copy()
    ad_sub.raw = ad_sub
    outfile = "adata-rawCount.h5ad"
    ad_sub.write(outfile)

    ## Normalization
    sc.pp.normalize_total(ad_sub, target_sum=1e4)
    sc.pp.log1p(ad_sub)

    ## identify HVGs
    sc.pp.highly_variable_genes(
        ad_sub,
        #flavor="seurat_v3",
        n_top_genes=2000,
        #layer="counts",
        subset=True
    )

    ## scVI training
    scvi.model.SCVI.setup_anndata(ad_sub, layer="counts", batch_key="Sample")
    vae = scvi.model.SCVI(ad_sub, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    vae_name = 'vae_model-bySample'
    vae.save(vae_name)

    ad_sub.obsm["X_scVI"] = vae.get_latent_representation()
    sc.pp.neighbors(ad_sub, use_rep="X_scVI")
    sc.tl.leiden(ad_sub, resolution = 0.2, key_added='leiden_0.2')
    sc.tl.leiden(ad_sub, resolution = 0.3, key_added='leiden_0.3')
    sc.tl.leiden(ad_sub, resolution = 0.4, key_added='leiden_0.4')
    sc.tl.leiden(ad_sub, resolution = 0.5, key_added='leiden_0.5')
    sc.tl.leiden(ad_sub, resolution = 0.6, key_added='leiden_0.6')
    sc.tl.leiden(ad_sub, resolution = 0.7, key_added='leiden_0.7')
    sc.tl.leiden(ad_sub, resolution = 0.8, key_added='leiden_0.8')
    sc.tl.leiden(ad_sub, resolution = 0.9, key_added='leiden_0.9')

    ad_sub.obsm["X_mde"] = mde(ad_sub.obsm["X_scVI"])
    ## find DEGs
    sc.tl.rank_genes_groups(ad_sub, groupby='leiden_0.2', key_added= 'rank_leiden_0.2' ,method='wilcoxon')
    sc.tl.rank_genes_groups(ad_sub, groupby='leiden_0.3', key_added= 'rank_leiden_0.3' ,method='wilcoxon')
    sc.tl.rank_genes_groups(ad_sub, groupby='leiden_0.4', key_added= 'rank_leiden_0.4' ,method='wilcoxon')
    sc.tl.rank_genes_groups(ad_sub, groupby='leiden_0.5', key_added= 'rank_leiden_0.5' ,method='wilcoxon')
    sc.tl.rank_genes_groups(ad_sub, groupby='leiden_0.6', key_added= 'rank_leiden_0.6' ,method='wilcoxon')
    sc.tl.rank_genes_groups(ad_sub, groupby='leiden_0.7', key_added= 'rank_leiden_0.7' ,method='wilcoxon')
    sc.tl.rank_genes_groups(ad_sub, groupby='leiden_0.8', key_added= 'rank_leiden_0.8' ,method='wilcoxon')
    sc.tl.rank_genes_groups(ad_sub, groupby='leiden_0.9', key_added= 'rank_leiden_0.9' ,method='wilcoxon')

    outfile = "adata-subset-scVI-bySample.h5ad"
    ad_sub.write(outfile)
    flaginfo2 = '>>>>>> ' + 'finished processing ' + celltype
    print(flaginfo2)
    print()
    os.chdir('..')
### For DEG checking
def get_top_genes_for_cluster1(adata, my_leiden='leiden_0.6', cluster ='1', top_n=15):
    my_rank = 'rank_' + my_leiden
    ranked_genes = adata.uns[my_rank]['names'][cluster]
    logfc = adata.uns[my_rank]['logfoldchanges'][cluster]
    pvals = adata.uns[my_rank]['pvals'][cluster]
    
    df = pd.DataFrame({
        'gene': ranked_genes,
        'logfc': logfc,
        'pval': pvals
    })
    df = df[df['pval'] < 0.05]
    df = df[df['logfc'] > 0.25]
    df_sorted = df.sort_values(by=['pval', 'logfc'], ascending=[True, False])
    return df_sorted.head(top_n)

def get_top_genes_for_cluster2(adata, my_leiden='leiden_0.6', cluster ='1', top_n=15):
    my_rank = 'rank_' + my_leiden
    ranked_genes = adata.uns[my_rank]['names'][cluster]
    logfc = adata.uns[my_rank]['logfoldchanges'][cluster]
    pvals = adata.uns[my_rank]['pvals'][cluster]

    df = pd.DataFrame({
        'gene': ranked_genes,
        'logfc': logfc,
        'pval': pvals
    })
    df = df[df['pval'] < 0.05]
    df = df[df['logfc'] > 0.25]
    df_sorted = df.sort_values(by=['logfc','pval'], ascending=[False,True])
    return df_sorted.head(top_n)

def save_top_genes_for_cluster(adata, my_leiden='leiden_0.6'):
    my_rank = 'rank_' + my_leiden
    genes_df = pd.DataFrame()
    leiden_clusters = adata.obs[my_leiden].unique()

    def custom_sort_key(item):
        prefix = ''.join(filter(str.isalpha, item))
        number_part = ''.join(filter(str.isdigit, item))
        return (prefix, int(number_part) if number_part else float('inf'))
    
    # 对列表进行排序
    leiden_clusters = sorted(leiden_clusters, key=custom_sort_key)

    for cluster in leiden_clusters:
        ranked_genes = adata.uns[my_rank]['names'][cluster]
        logfc = adata.uns[my_rank]['logfoldchanges'][cluster]
        pvals = adata.uns[my_rank]['pvals'][cluster]
        df = pd.DataFrame({
            'gene': ranked_genes,
            'logfc': logfc,
            'pval': pvals,
            'cluster': cluster  # 新增一列记录亚群名称
        })
        df_sorted = df.sort_values(by=['pval', 'logfc'], ascending=[True, False])
        genes_df = pd.concat([genes_df, df_sorted], ignore_index=True)  # 使用pd.concat替代append

    output_file1 = my_leiden + '-genes_per_cluster.csv'
    genes_df.to_csv(output_file1, index=False)
    
    genes_df = genes_df[(genes_df['pval'] < 0.05) & (genes_df['logfc'] > 0.25)]
    output_file2 = my_leiden + '-degs_per_cluster.csv'
    genes_df.to_csv(output_file2, index=False)
def plot_dotplot_for_leiden(adata, my_leiden, cmap = 'bwr', vmax=None, vmin=None, dot_max=1, dot_min=0,width=10, high=6):
    my_rank = 'rank_' + my_leiden
    top5_genes = []
    
    leiden_clusters = adata.obs[my_leiden].unique()
    #leiden_clusters = sorted(list(leiden_clusters), key=int)

    # 自定义排序函数
    def custom_sort_key(item):
        # 分离字母和数字
        prefix = ''.join(filter(str.isalpha, item))
        number_part = ''.join(filter(str.isdigit, item))
        
        # 返回元组，首先按字母部分排序，然后按数字部分排序
        return (prefix, int(number_part) if number_part else float('inf'))
    
    # 对列表进行排序
    leiden_clusters = sorted(leiden_clusters, key=custom_sort_key)

    for cluster in leiden_clusters:
        ranked_genes = adata.uns[my_rank]['names'][cluster]
        logfc = adata.uns[my_rank]['logfoldchanges'][cluster]
        pvals = adata.uns[my_rank]['pvals'][cluster]
        
        df = pd.DataFrame({
            'gene': ranked_genes,
            'logfc': logfc,
            'pval': pvals
        })
        df = df[df['pval'] < 0.05]
        df = df[df['logfc'] > 0.25]
        df_sorted = df.sort_values(by=['pval', 'logfc'], ascending=[True, False])
        
        top5_genes.append(df_sorted['gene'].head(4).tolist())
    
    top5_genes_flat = [gene for sublist in top5_genes for gene in sublist]
    def deduplicate(lst):
        return list(OrderedDict.fromkeys(lst))
    top5_genes_flat = deduplicate(top5_genes_flat)
    sc.pl.dotplot(adata, 
                  var_names=top5_genes_flat,  # 使用平面列表中的top5基因
                  groupby=my_leiden,  # 使用提供的leiden键
                  cmap= cmap, 
                  standard_scale='var',  # 对基因表达进行标准化
                  expression_cutoff=1,
                  figsize=(width, high),
                  dot_max=dot_max, dot_min=dot_min,
                  vmax = vmax, vmin = vmin, 
                  use_raw=True,save='subcelltype.png')  # 只绘制表达量大于1的基因

    sc.pl.dotplot(adata, 
                  var_names=top5_genes_flat,  # 使用平面列表中的top5基因
                  groupby=my_leiden,  # 使用提供的leiden键
                  cmap= cmap, #'RdBu_r', # 'viridis', 
                  standard_scale='var',  # 对基因表达进行标准化
                  expression_cutoff=1,
                  figsize=(width, high),
                  dot_max=dot_max, dot_min=dot_min,
                  vmax = vmax, vmin = vmin, 
                  use_raw=True,save='subcelltype.pdf')  # 只绘制表达量大于1的基因
os.chdir('../03_subset-bySample/Myeloid cell/')
adata = sc.read_h5ad("adata-subset-scVI-bySample.h5ad")

sc.pl.embedding(
    adata,
    basis="X_mde",
    color=['leiden_0.2', 'leiden_0.3', 'leiden_0.4', 'leiden_0.5', 'leiden_0.6', 'leiden_0.7', 'leiden_0.8', 'leiden_0.9'],
    use_raw=True,
    size=1,
    legend_loc = 'on data',
    ncols=3,
    frameon=False)
marker_genes_dict = {
    'Mac': ["C1QA", "C1QB", "C1QC", "RNASE1", "DAB2", "LGMN", "PLTP", "MAF", "SLCO2B1"],  # 巨噬细胞
    'mono': ["VCAN", "FCN1", "CD300E", "S100A12", "EREG", "APOBEC3A", "STXBP2", "ASGR1", "CCR2", "NRG1"],  # 单核细胞
    'neutrophils': ["FCGR3B", "CXCR2", "SLC25A37", "G0S2", "CXCR1", "ADGRG3", "PROK2", "STEAP4", "CMTM2"],  # 中性粒细胞
    'pDC': ["GZMB", "SCT", "CLIC3", "LRRC26", "LILRA4", "PACSIN1", "CLEC4C", "MAP1A", "PTCRA", "C12orf75"],  # 浆细胞样树突状细胞
    'DC1': ["CLEC9A", "XCR1", "CLNK", "CADM1", "ENPP1", "SNX22", "NCALD", "DBN1", "HLA-DOB", "PPY"],  # cDC1
    'DC2': ["CD1C", "FCER1A", "CD1E", "CD2", "GPAT3", "CCND2", "ENHO", "PKIB", "CD1B"],  # cDC2
    'DC3': ["HMSD", "ANKRD33B", "LAD1", "CCR7", "LAMP3", "CCL19", "CCL22", "INSM1", "TNNT2", "TUBB2B"],  # 其他或炎症相关DC
    'Mast': ["KIT","TPSAB1","CPA3"],
    'Mono CD14+': ["CD14","S100A9","S100A8"],
    'Mono CD16+':  ["FCGR3A","LST1","LILRB2"]
}
sc.pl.dotplot(
    adata, 
    marker_genes_dict, 
    groupby='leiden_0.8', 
    use_raw=True)
from matplotlib.colors import LinearSegmentedColormap

colors = ["#364398", '#FFECA4',"#AE0927"]
cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=100)

my_leiden = 'leiden_0.8'
plot_dotplot_for_leiden(adata,  my_leiden, cmap = cmap,width=18, high=6)
#def plot_dotplot_for_leiden(adata, my_leiden, cmap = 'bwr', vmax=None, vmin=None, dot_max=1, dot_min=0,width=10, high=6):
annotation = {
    'mono_c1': [1,15],
    'mono_c2': [2],
    'mono_c3': [10],
    'mph_c1': [0],
    'mph_c2': [3],
    'mph_c3': [5,8],
    'cd3': [6],
    'mph_c5': [7],
    'mph_c6':[9],
    'krt':[12],
    'bcell':[13],
    'low':[16],
    'dc_c1': [4],
    'dc_c2': [11],
    'dc_c3': [14]
}
annot_dict = {str(c): ct for ct, clusters in annotation.items() for c in clusters}

adata.obs["celltype_tmp"] = [annot_dict.get(c, "unknown") for c in adata.obs["leiden_0.8"]]
sc.pl.embedding(
    adata,
    basis="X_mde",
    color=['celltype_tmp'],
    frameon=False,
    size=1,
    use_raw =True,
    legend_loc = 'on data',
    #palette=palette  # 使用自定义的调色板
    #save='_Totalcell-Patient.png'
)
celltypes_to_remove = ['cd3', 'krt', 'bcell', 'low','mono-unknown']
adata = adata[~adata.obs['celltype_tmp'].isin(celltypes_to_remove)].copy()

sc.pl.embedding(
    adata,
    basis="X_mde",
    color=['celltype_tmp'],
    frameon=False,
    size=1,
    use_raw =True,
    legend_loc = 'on data',
    #palette=palette  # 使用自定义的调色板
    #save='_Totalcell-Patient.png'
)
## import vae model 
import scvi
vae = scvi.model.SCVI.load('vae_model-bySample', adata)

## scANVI integrate
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="celltype_tmp",
    unlabeled_category="unknown",
)

lvae.train(max_epochs=20, n_samples_per_label=100)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
#adata.obs["C_scANVI"] = lvae.predict(adata)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)
sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)
sc.pl.umap(adata, color=['celltype_tmp', "Class"], wspace=0.45)
sc.pl.umap(adata, color=['celltype_tmp'],legend_loc = 'on data',wspace=0.45)
## identify DEGs
sc.tl.rank_genes_groups(adata, groupby='celltype_tmp', key_added= 'rank_celltype_tmp' ,method='wilcoxon')
from matplotlib.colors import LinearSegmentedColormap

colors = ["#364398", '#FFECA4',"#AE0927"]
cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=100)

my_leiden = 'celltype_tmp'
plot_dotplot_for_leiden(adata,  my_leiden, cmap = cmap,width=18, high=6)
annotation = {
    'Mono_FCN1': ['mono_c1'],
    'Mono_EREG': ['mono_c2'], 
    'Mono_HSPA1A':['mono_c3'],
    'Macro_APOE': ['mph_c1'],
    'Macro_F13A1':['mph_c2'],
    'Macro_SPP1': ['mph_c3'],
    'Macro_MKI67':['mph_c5'],
    'DC_FCGR1A':['dc_c1'],
    'DC_CCL19': ['dc_c2'],
    'DC_CLEC9A':['dc_c3']
}
annot_dict = {str(c): ct for ct, clusters in annotation.items() for c in clusters}

adata.obs["subcelltype"] = [annot_dict.get(c, "unknown") for c in adata.obs["celltype_tmp"]]
sc.tl.rank_genes_groups(adata, groupby='subcelltype', key_added= 'rank_subcelltype' ,method='wilcoxon')
palette=["#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#7E6148FF", "#B09C85FF", "#EFC00099", "#86868699", "#CD534C99", "#7AA6DC99"]

sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['subcelltype'],
    #use_raw=True,
    #frameon=False,
    palette = palette,
    size=5,save='_subcelltype.pdf')
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['subcelltype'],
    #use_raw=True,
    #frameon=False,
    palette = palette,
    size=5,save='_subcelltype.png')

sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['subcelltype'],
    #use_raw=True,
    size=5,
    legend_loc = 'on data',
    #frameon=False,
    palette = palette,
    save='_subcelltype-label.pdf')
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['subcelltype'],
    #use_raw=True,
    size=5,
    legend_loc = 'on data',
    #frameon=False,
    palette = palette,
    save='_subcelltype-label.png')
from matplotlib.colors import LinearSegmentedColormap

colors = ["#364398", '#FFECA4',"#AE0927"]
cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=100)

my_leiden = 'subcelltype'
plot_dotplot_for_leiden(adata,  my_leiden, cmap = cmap, width=12,high=3.5,vmin=0,vmax=1,dot_min=0.1)
#def plot_dotplot_for_leiden(adata, my_leiden, cmap = 'bwr', vmax=None, vmin=None, dot_max=1, dot_min=0,width=10, high=6):

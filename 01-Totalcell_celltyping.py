## load packages
import scanpy as sc
import scvi
from rich import print
from scib_metrics.benchmark import Benchmarker
from scvi.model.utils import mde
import numpy as np
import pymde
import os

## set working directory
#os.mkdir('02_scVI_Integration-bySample') # exist
os.chdir('02_scVI_Integration-bySample') 

## read adata
file_path = "../01_merge_data/adata-filtered.h5ad"
adata = sc.read_h5ad(file_path)

## remove noise
adata.obs.drop(columns=['predicted_doublet','log1p_total_counts_ribo','pct_counts_ribo','total_counts_hb','log1p_total_counts_hb','pct_counts_hb'], inplace=True)
adata.obs.drop(columns=['total_counts_ribo','pct_counts_in_top_500_genes','pct_counts_in_top_50_genes','pct_counts_in_top_100_genes','pct_counts_in_top_200_genes'], inplace=True)
adata.obs.drop(columns=['log1p_n_genes_by_counts','log1p_total_counts','doublet_score'], inplace=True)

## add meta information
meta_file = "../01_merge_data/meta-sample-patient-pair.csv"  # 请替换成实际文件路径
meta_data = pd.read_csv(meta_file, index_col=0)
sample_to_patient = dict(zip(meta_data['Sample'], meta_data['Patient']))
sample_to_match = dict(zip(meta_data['Sample'], meta_data['Match']))

adata.obs['Patient'] = None
adata.obs['Match'] = None
def fill_patient(sample):
    return sample_to_patient.get(sample, None)

def fill_match(sample):
    return sample_to_match.get(sample, None)
adata.obs['Patient'] = adata.obs['Sample'].apply(fill_patient)
adata.obs['Match'] = adata.obs['Sample'].apply(fill_match)

## remove sample with low number cells
sample_cell_counts = adata.obs['Sample'].value_counts()
samples_to_keep = sample_cell_counts[sample_cell_counts >= 200].index
adata = adata[adata.obs['Sample'].isin(samples_to_keep), :]
samples_to_modify = ['KUL31-N', 'C119_T']
adata.obs.loc[adata.obs['Sample'].isin(samples_to_modify), 'Match'] = 'nopair'

## sample raw data in layers
adata.layers["counts"] = adata.X.copy()  
adata.raw = adata
adata.write("adata-rawCount.h5ad")

adata = sc.read_h5ad('adata-rawCount.h5ad')
## Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

## identify HVGs
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    subset=True)

## scVI training
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="Sample")
vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
vae.train()
vae.save('vae_model-bySample')

adata.obsm["X_scVI"] = vae.get_latent_representation()
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.leiden(adata, resolution = 0.1, key_added='leiden_0.1')
sc.tl.leiden(adata, resolution = 0.2, key_added='leiden_0.2')
sc.tl.leiden(adata, resolution = 0.3, key_added='leiden_0.3')
sc.tl.leiden(adata, resolution = 0.4, key_added='leiden_0.4')
sc.tl.leiden(adata, resolution = 0.5, key_added='leiden_0.5')
sc.tl.leiden(adata, resolution = 0.6, key_added='leiden_0.6')
adata.obsm["X_mde"] = mde(adata.obsm["X_scVI"])

## find DEGs
sc.tl.rank_genes_groups(adata, groupby='leiden_0.1', key_added= 'rank_leiden_0.1' ,method='wilcoxon')
sc.tl.rank_genes_groups(adata, groupby='leiden_0.2', key_added= 'rank_leiden_0.2' ,method='wilcoxon')
sc.tl.rank_genes_groups(adata, groupby='leiden_0.3', key_added= 'rank_leiden_0.3' ,method='wilcoxon')
sc.tl.rank_genes_groups(adata, groupby='leiden_0.4', key_added= 'rank_leiden_0.4' ,method='wilcoxon')
sc.tl.rank_genes_groups(adata, groupby='leiden_0.5', key_added= 'rank_leiden_0.5' ,method='wilcoxon')
sc.tl.rank_genes_groups(adata, groupby='leiden_0.6', key_added= 'rank_leiden_0.6' ,method='wilcoxon')
adata.write("adata-scVI.h5ad")

### For DEG checking
from collections import OrderedDict
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

## for DotPlot
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
        
        top5_genes.append(df_sorted['gene'].head(5).tolist())
    
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
os.chdir('/share/home/wuqfLab/xueyu/jupyter/CRC/02_scVI_Integration-bySample')
adata = sc.read_h5ad('adata-scVI.h5ad')
adata.obs

plot_dotplot_for_leiden(adata, 'leiden_0.2', cmap = 'bwr', vmax=None, vmin=None, dot_max=1, dot_min=0,width=15, high=6)
print('leiden_0.2')
sc.pl.dotplot(
    adata, 
    marker_genes_dict, 
    groupby='leiden_0.2', 
    use_raw=True)
sc.pl.embedding(
    adata,
    basis="X_mde",
    color=['leiden_0.2'],
    use_raw=True,
    size=1,
    legend_loc = 'on data')
## annotation
annotation = {
    'T cell': [0],
    'B cell': [5,9],
    'Plasma cell': [3,12],
    'Mast cell':[10],
    'Myeloid cell':[2],
    'Endothelial cell': [8],
    'Fibroblast': [7],
    'Glial cell':[13],
    'Epithelial cell':[1, 4, 6, 11, 14]
}
annot_dict = {str(c): ct for ct, clusters in annotation.items() for c in clusters}
adata.obs["celltype"] = [annot_dict.get(c, "unknown") for c in adata.obs["leiden_0.2"]]
adata.obs["celltype_unknown"] = [
    "known" if ct != "unknown" else ct for ct in adata.obs["celltype"]
]

## plot
desired_order = ['T cell','B cell','Plasma cell','Mast cell','Myeloid cell','Endothelial cell','Fibroblast','Glial cell','Epithelial cell']
adata.obs['celltype'] = pd.Categorical(adata.obs['celltype'], categories=desired_order, ordered=True)
sc.pl.embedding(
    adata,
    basis="X_mde",
    color=['celltype'],
    use_raw=True,
    size=1,
    legend_loc = 'on data')
marker_genes_dict = {
    'T cell':['CD3D','CD3E','CD8A','CD8B','CD2'], # T cell
    'B cell':['CD19','MS4A1','CD22'], # B cell
    'Plasma cell':['CD79A','JCHAIN','MZB1'], # Plasma cell # 'JCHAIN',
    'Myeloid cell':['CD14','C1QA','C1QB','MRC1','LYZ','FCGR3A'], # Myeloid cell
    'Cycling cell':['TOP2A','MKI67'],
    'Endothelial cell':['CLDN5','VWF'], # Endothelial cell
    'Fibroblast':['MMP11','COL10A1','THY1','DCN','ACTA2'], # Fibroblast
    'Epithelial cell':['EPCAM','KRT18','KRT19'], # epithelial cell
    'Mast cell':['TPSAB1','HDC','CTSG','CMA1','KRT1','IL1RAPL1','GATA2'], # Mast
    'NK Cell': ['FGFBP2', 'FCGR3A', 'CX3CR1', 'GNLY', 'NKG7', 'TYROBP', 'PRF1'],  # NK cell
    'pDC': ['CLEC4C', 'SLC32A1', 'LILRA4'],  # pDC
    'Glial cell':['PLP1','SOX10','PTPRZ1']# Glial cell 
}
## plot dotplot for celltype
# 绘制 MDE 图 - 细胞类型分布
marker_genes_dict = {
    'T cell':['CD3D','CD2','CD8A'], # T cell
    'B cell':['MS4A1','CD19','CD22'], # B cell
    'Plasma cell':['CD79A','JCHAIN','MZB1'], # Plasma cell
    'Mast cell':['TPSAB1','CTSG','GATA2'], # Mast
    'Myeloid cell':['CD14','FCGR3A','C1QA'], # Myeloid cell
    'Endothelial cell':['CLDN5','VWF'], # Endothelial cell
    'Fibroblast':['THY1','DCN','ACTA2'], # Fibroblast
    'Glial cell':['PLP1','SOX10','PTPRZ1'], # Glial cell 
    'Epithelial cell':['EPCAM','KRT18','KRT19'] # epithelial cell
}
sc.pl.dotplot(
    adata, 
    marker_genes_dict, 
    groupby='celltype', 
    use_raw=True,vmax=0.9)
adata.write("adata-scVI-anno.h5ad")
sc.pp.neighbors(adata, use_rep="X_scVI")
sc.tl.umap(adata)
sc.pl.umap(adata, color=['celltype', "Class"], wspace=0.45)

adata = sc.read_h5ad('adata-scVI-anno.h5ad')
## import vae model 
import scvi
vae = scvi.model.SCVI.load('vae_model-bySample', adata)

##scANVI integrate
lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=adata,
    labels_key="celltype",
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=20, n_samples_per_label=100)
adata.obsm["X_scANVI"] = lvae.get_latent_representation(adata)

sc.pp.neighbors(adata, use_rep="X_scANVI")
sc.tl.umap(adata)
sc.pl.umap(adata, color=['celltype', "Class"], wspace=0.45)
adata.write("adata-scANVI.h5ad")

import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

## dot plot for classical marker genes in each celltype
marker_genes_dict_refined = {
    'T cell':['CD3D','CD3E','CD8A'], # T cell
    'B cell':['MS4A1','CD19'], # B cell
    'Plasma cell':['CD79A','JCHAIN','MZB1'], # Plasma cell
    'Mast cell':['TPSAB1','GATA2'], # Mast
    'Myeloid cell':['CD14','FCGR3A','C1QA'], # Myeloid cell
    'Endothelial cell':['CLDN5','VWF'], # Endothelial cell
    'Fibroblast':['THY1','DCN','ACTA2'], # Fibroblast
    'Glial cell':['PLP1','SOX10','PTPRZ1'],# Glial cell 
    'Epithelial cell':['EPCAM','KRT18','KRT19'], # epithelial cell
}
cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=100)

sc.pl.dotplot(adata, 
              var_names=marker_genes_dict_refined,  # 使用平面列表中的top5基因
              groupby='celltype', 
              cmap= cmap, #'viridis', 
              standard_scale='var',  # 对基因表达进行标准化
              #categories_order = sort_celltype,
              #expression_cutoff=1,
              mean_only_expressed = True,
              dendrogram=False,use_raw=True,save='TotalCell-celltype.pdf')
sc.pl.dotplot(adata, 
              var_names=marker_genes_dict_refined,  # 使用平面列表中的top5基因
              groupby='celltype', 
              cmap= cmap, #'viridis', 
              standard_scale='var',  # 对基因表达进行标准化
              #categories_order = sort_celltype,
              #expression_cutoff=1,
              mean_only_expressed = True,
              dendrogram=False,use_raw=True,save='TotalCell-celltype.png')

## x_MDE for celltype
desired_order = ['T cell','B cell','Plasma cell','Mast cell','Myeloid cell','Endothelial cell','Fibroblast','Glial cell','Epithelial cell']
adata.obs['celltype'] = pd.Categorical(adata.obs['celltype'], categories=desired_order, ordered=True)
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['celltype'],
    #use_raw=True,
    frameon=False,
    palette = palette,
    size=0.8,save='_Totalcell-celltype.pdf'
)
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['celltype'],
    #use_raw=True,
    frameon=False,
    palette = palette,
    size=0.8,save='_Totalcell-celltype.png'
)

sc.pl.umap(adata, color=['celltype'], wspace=0.45,palette = palette,frameon=False,size=0.6)

## save DEGs for major cell types
result = adata.uns['rank_celltype']  # 获取 rank_genes_groups 结果
groups = result['names'].dtype.names  # 获取所有分组名（即 sub_celltype 中的不同亚簇）

dfs = []
for group in groups:
    df = pd.DataFrame({
        'gene': result['names'][group],  # 基因名称
        'logfc': result['logfoldchanges'][group],  # logfoldchange
        'pval': result['pvals'][group],  # p-value
        'pval_adj': result['pvals_adj'][group],  # 调整后的 p-value
        'scores': result['scores'][group]  # 评分
    })
    df = df.sort_values(by=['pval','logfc'], ascending=[True, False])
    df['group'] = group  # 添加 'group' 列标记当前亚簇
    dfs.append(df)  # 将每个亚簇的数据添加到列表中
final_df = pd.concat(dfs, ignore_index=True)
output_path = "TotalCell-celltype-DEGs.csv"
final_df.to_csv(output_path, index=False)

### Class
print('Class:')
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['Class'],
    #use_raw=True,
    frameon=False,
    palette = ["#00A087FF","#3C5488FF", "#F39B7FFF"],
    size=1,save='_Totalcell-Class.pdf'
)
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['Class'],
    #use_raw=True,
    frameon=False,
    palette = ["#00A087FF","#3C5488FF", "#F39B7FFF"],
    size=1,save='_Totalcell-Class.png'
)

### Batch
print('Batch:')
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['batch'],
    #use_raw=True,
    frameon=False,
    palette = palette,
    size=1,save='_Totalcell-batch.pdf'
)
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['batch'],
    #use_raw=True,
    frameon=False,
    palette = palette,
    size=1,save='_Totalcell-batch.png'
)
print('Patient:')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=102)
palette = [cmap(i) for i in range(cmap.N)]

### Patient
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['Patient'],
    frameon=False,
    size=1,
    palette=palette,  # 使用自定义的调色板
    save='_Totalcell-Patient.pdf'
)
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['Patient'],
    frameon=False,
    size=1,
    palette=palette,  # 使用自定义的调色板
    save='_Totalcell-Patient.png'
)
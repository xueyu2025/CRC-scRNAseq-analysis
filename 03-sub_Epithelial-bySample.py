## load packages
import scanpy as sc
from rich import print
from collections import OrderedDict
import scvi
from scvi.model.utils import mde
import numpy as np
import pandas as pd
import os
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

def plot_dotplot_for_leiden_fc(adata, my_leiden, cmap = 'bwr', vmax=None, vmin=None, dot_max=1, dot_min=0,width=10, high=6):
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
        df_sorted = df.sort_values(by=['logfc','pval'], ascending=[False, True])
        
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
                  use_raw=True,save='subcelltype-fc.png')  # 只绘制表达量大于1的基因

    sc.pl.dotplot(adata, 
                  var_names=top5_genes_flat,  # 使用平面列表中的top5基因
                  groupby=my_leiden,  # 使用提供的leiden键
                  cmap= cmap, #'RdBu_r', # 'viridis', 
                  standard_scale='var',  # 对基因表达进行标准化
                  expression_cutoff=1,
                  figsize=(width, high),
                  dot_max=dot_max, dot_min=dot_min,
                  vmax = vmax, vmin = vmin, 
                  use_raw=True,save='subcelltype-fc.pdf')  # 只绘制表达量大于1的基因


adata = sc.read_h5ad("adata-subset-scVI-sample.h5ad")
sc.pl.embedding(
    adata,
    basis="X_mde",
    color=['leiden_0.1', 'leiden_0.2', 'leiden_0.3', 'Class', 'leiden_0.4', 'leiden_0.5', 'leiden_0.6', 'leiden_0.7'],
    use_raw=True,
    size=1,
    legend_loc = 'on data',
    ncols=3,
    frameon=False)
## import vae model 
import scvi
vae = scvi.model.SCVI.load('vae_model-bySample', adata)
adata.obsm["X_scVI"] = vae.get_latent_representation(adata)
sc.pl.embedding(
    adata,
    basis="X_umap",
    color=['leiden_0.1', 'leiden_0.2', 'leiden_0.3', 'Class', 'leiden_0.4', 'leiden_0.5', 'leiden_0.6', 'leiden_0.7'],
    use_raw=True,
    size=1,
    legend_loc = 'on data',
    ncols=3,
    frameon=False)
cluster_counts = adata.obs['leiden_0.7'].value_counts()
keep_clusters = cluster_counts[cluster_counts >= 100].index.tolist()
keep_mask = adata.obs['leiden_0.7'].isin(keep_clusters)
adata = adata[keep_mask, :].copy()

adata.obs['subcelltype'] = adata.obs['leiden_0.7'].astype(str).copy()
sc.pl.embedding(
    adata,
    basis="X_umap",
    palette = palette,
    color=['subcelltype', 'Class','batch'],
    use_raw=True,
    size=1,
    legend_loc = 'on data',
    ncols=3,
    frameon=False)

adata.obs['subcelltype'] = adata.obs['leiden_0.7'].astype(str).copy()
# 替换并更新分类类别
adata.obs['subcelltype'] = adata.obs['subcelltype'].replace('16', '7')
adata.obs['subcelltype'] = adata.obs['subcelltype'].replace('0', '3')
#adata.obs['subcelltype'] = adata.obs['subcelltype'].replace('0', '1')
adata.obs['subcelltype'].value_counts()

annotation = {
    'Epi01': ['3','4'],
    'Epi02': ['1'], 
    'Epi03':['2'],
    'Epi04': ['5'],
    'Epi05':['6'],
    'Epi06': ['7'],
    'Epi07':['8'],
    'Epi08':['9'],
    'Epi09': ['10'],
    'Epi10':['11']
}
annot_dict = {str(c): ct for ct, clusters in annotation.items() for c in clusters}

adata.obs["subcelltype"] = [annot_dict.get(c, "unknown") for c in adata.obs["subcelltype"]]

print(adata.obs['subcelltype'].value_counts())

sc.pl.embedding(
    adata,
    basis="X_umap",
    palette = palette,
    color=['subcelltype', 'Class','batch'],
    use_raw=True,
    size=1,
    legend_loc = 'on data',
    ncols=3,
    frameon=False)

adata.write('adata-rawCount-Epi.h5ad')
# 先确保 subcelltype 是字符串类型
adata.obs['subcelltype'] = adata.obs['subcelltype'].astype(str)
# 过滤掉 subcelltype 中的细胞
adata_filtered = adata[~adata.obs['subcelltype'].isin(['Epi12', 'Epi13', 'Epi14']), :].copy()
sc.pl.embedding(
    adata_filtered,
    basis="X_umap",
    palette = palette,
    color=['subcelltype', 'Class','batch'],
    use_raw=True,
    size=1,
    legend_loc = 'on data',
    ncols=3,
    frameon=False)


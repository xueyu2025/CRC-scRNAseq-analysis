import os
import scanpy as sc
import pandas as pd
import numpy as np
adata = sc.read_h5ad('../../02_scVI_Integration-bySample/adata-rawCount-meta.h5ad')
epi_meta = pd.read_csv("./adata-epithelialcell-meta.csv", index_col=0)
mye_meta = pd.read_csv("../Myeloid cell/adata_myeloidcell-meta.csv", index_col=0)
adata.obs['subcelltype'] = pd.NA
common_epi = adata.obs_names.intersection(epi_meta.index)
adata.obs.loc[common_epi, 'subcelltype'] = epi_meta.loc[common_epi, 'subcelltype']
common_mye = adata.obs_names.intersection(mye_meta.index)
adata.obs.loc[common_mye, 'subcelltype'] = mye_meta.loc[common_mye, 'subcelltype']
adata_subset = adata[adata.obs['subcelltype'].notna()].copy()
print(adata.obs['subcelltype'].value_counts())
adata_raw = sc.read_h5ad('../../02_scVI_Integration-bySample/adata-rawCount.h5ad')
adata_cells = adata_subset.obs.index
adata_raw_subset = adata_raw[adata_raw.obs.index.isin(adata_cells)]
adata_raw_subset.obs['subcelltype'] = adata_raw_subset.obs.index.map(lambda x: adata_subset.obs['subcelltype'][x])
mat = pd.DataFrame(data=adata_raw_subset.X.todense(), 
                   index=adata_raw_subset.obs_names,  
                   columns=adata_raw_subset.var_names)  
mat.to_hdf("adata_Epi-Mye-mat.h5ad", key="mat")  
adata_raw_subset.obs.to_csv('adata_Epi-Mye-meta.csv')
adata_subset.write('adata-Epi-Mye.h5ad')
immunecell = ['T cell', 'Myeloid cell', 'Fibroblast']
adata_immu = adata[adata.obs['celltype'].isin(immunecell)].copy()
sampled_cells = (
    adata_immu.obs
    .groupby('celltype', group_keys=False)
    .apply(lambda x: x.sample(n=min(500, len(x)), random_state=42))
)
adata_ref = adata_immu[sampled_cells.index].copy()
adata_combined = ad.concat([adata_ref, adata_epi], axis=0, join='outer')
adata_combined.write('TotalCell-sampled-for-inferCNV.h5ad')
df = pd.read_csv('adata-epithelialcell-meta.csv') 
print(df.head()) 
import pandas as pd
np.random.seed(42)
sampled_df = df.groupby('subcelltype', group_keys=False).apply(
    lambda x: x.sample(n=min(500, len(x)), random_state=42)
)
print("抽样后各类别数量：\n", sampled_df['subcelltype'].value_counts())
adata_epi = adata[sampled_df.iloc[:,0]].copy()
cell_names = sampled_df.iloc[:, 0].tolist()
missing_cells = set(cell_names) - set(adata_epi.obs.index)
if missing_cells:
    print(f"警告：以下细胞在 adata_epi 中不存在，将被忽略: {missing_cells}")
subcelltype_series = pd.Series(
    data=sampled_df['subcelltype'].values,
    index=sampled_df.iloc[:, 0],  # 索引为细胞名称
    name='subcelltype'  # 新列名
)
subcelltype_series = subcelltype_series[subcelltype_series.index.isin(adata_epi.obs.index)]
adata_ref.obs['subcelltype'] = adata_ref.obs['celltype']
adata = adata_combined.copy()
adata.layers["counts"] = adata.X
import pandas as pd
raw_counts = pd.DataFrame(
    adata.layers["counts"].toarray(),
    index=adata.obs_names,  
    columns=adata.var_names  
)
raw_counts = raw_counts.astype(int)
raw_counts = raw_counts.T
raw_counts.to_csv("infercnv/raw_counts_matrix.txt", sep="\t")
cell_annotations = adata.obs[["subcelltype"]]  
cell_annotations["group"] = cell_annotations["subcelltype"] # ["Normal" if x in reference_cat else "Tumor" for x in cell_annotations["subcelltype"]]
cell_annotations.to_csv("infercnv/cell_annotations.txt", sep="\t", header=False)
sc.pp.normalize_total(adata, target_sum=None) 
sc.pp.log1p(adata)
print("矩阵最大值:", adata.X.max()) 
gene_pos = pd.read_csv("gene_filter_position_rmdup.txt", sep="\t", names=["gene", "chr", "start", "end"])
gene_pos.set_index("gene", inplace=True)
adata.var["chromosome"] = gene_pos.loc[adata.var_names, "chr"].values
adata.var["start"] = gene_pos.loc[adata.var_names, "start"].values
adata.var["end"] = gene_pos.loc[adata.var_names, "end"].values
adata = adata[:, adata.var.chromosome.notna()]
reference_cat = ["T cell", "Myeloid cell", "Fibroblast"]  
cnv.tl.infercnv(
    adata,
    reference_key="subcelltype",       
    reference_cat=reference_cat,    
    exclude_chromosomes = ['MT'],
    n_jobs = 1,
    window_size=250                
)
adata.write("results_cnv.h5ad")
library(scMetabolism)
library(Seurat)
library(qs)
library(tidyverse)
library(rhdf5)
library(Matrix)
library(ggsci)

sce = qread('epiCell-PCA.qs')
CMSmeta = read.csv('../../02_scVI_Integration-bySample/CMS/TotalCell-Seurat-meta-CMS.csv',row.names = 1)
CMSmeta = CMSmeta %>% dplyr::select(Sample,CMS) %>% unique()
meta = sce@meta.data
meta$cellid = rownames(meta)
meta.1 = merge(meta,CMSmeta,by='Sample')
rownames(meta.1) = meta.1$cellid
meta.1 = meta.1[rownames(meta),]
sce@meta.data = meta.1
sce$Class_subcelltype = paste0(sce$Class, '_',sce$subcelltype)
head(sce@meta.data$Class_subcelltype)
cell_counts <- table(sce$Class_subcelltype)
valid_class <- names(cell_counts[cell_counts > 2])
sce <- subset(sce, Class_subcelltype %in% valid_class)
sce <- sce %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(features = rownames(sce))
sce =  RunPCA(sce, features = VariableFeatures(sce))
#Idents(sce) <- "Class_subcelltype"
countexp.Seurat <- sc.metabolism.Seurat(obj = sce,  #Seuratde单细胞object
                                      method = "AUCell",
                                      imputation = F,
                                      ncores = 10,
                                      metabolism.type = "KEGG")
score <- countexp.Seurat@assays$METABOLISM$score
score_change <- score %>%
  select_all(~str_replace_all(., "\\.", "-"))
print(identical(colnames(score_change) , rownames(countexp.Seurat@meta.data)))
score <- countexp.Seurat@assays$METABOLISM$score
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change) )
qsave(countexp.Seurat, 'EpithelialCell-scMetabolism.qs')
scm = rownames(sce@assays$METABOLISM$score)
meta_data_filtered <- sce@assays$METABOLISM$score #[scm,] #sce.sub@meta.data[, scm]
sorted_pathway = sort(apply(meta_data_filtered,1,mean))[seq(82/2+1,82)]
sorted_pathway
pathway_half = names(sorted_pathway)
for (pathway
     in pathway_half){
   # pathway = "Glycolysis / Gluconeogenesis"
    p <- VlnPlot(sce,features = pathway,group.by='subcelltype',pt.size = 0,cols = colorseq) + NoLegend() + theme(axis.title = element_blank())
    pathwayName = gsub(' / ','_',pathway)
    pathwayName = gsub(' ','_',pathwayName)
    print(pathwayName)
    ggsave(filename = paste0('metabolism-pathway-subcelltype/',pathwayName,'.pdf'),plot = p,width=5,height = 3.5)
    ggsave(filename = paste0('metabolism-pathway-subcelltype/',pathwayName,'.jpeg'),plot = p,width=5,height=3.5)
    }
Idents(sce) <- "subcelltype"
exclude_groups <- c("Epi03","Epi05","Epi07","Epi08")
sce_sub <- subset(sce, idents = exclude_groups, invert = TRUE)
scts = unique(sce$subcelltype)
scts = setdiff(scts, c('Epi03','Epi05','Epi07','Epi08'))
sce.sub = subset(sce, subcelltype %in% scts)
meta = sce.sub@meta.data
library(dplyr)
meta.rank = meta %>% dplyr::select(subcelltype,`Glycolysis / Gluconeogenesis`:`Drug metabolism - other enzymes`)library(dplyr)
meta.rank = meta %>% dplyr::select(subcelltype,`Glycolysis / Gluconeogenesis`:`Drug metabolism - other enzymes`)
meta.rank = data.table::melt(meta.rank)
meta.rank = meta.rank[meta.rank$variable %in% pathway_half,]
meta.rank = meta.rank %>% group_by(subcelltype) %>% summarize(mean(value)) %>% arrange(desc(`mean(value)`))
seq = rev(meta.rank$subcelltype)
## DotPlot
sce.sub$subcelltype = factor(sce.sub$subcelltype,levels = seq) # 指定Class_subcelltype的顺序
p=DotPlot(
    sce.sub,
    features = pathway_half,  
    group.by = "subcelltype",                   
    dot.scale = 6,                            
    cluster.idents =FALSE) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank())+  # 旋转X轴标签
          #plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")+
    scale_color_gradientn(values = seq(0,1,0.2), colours = colorRampPalette(c("#364398", '#FFECA4',"#AE0927"))(n = 500));p
ggsave('EpithelialCell-tumor-subcelltype-DotPlot-metabolism.pdf',width=14,height = 5.5)
ggsave('EpithelialCell-tumor-subcelltype-DotPlot-metabolism.png',width=14,height = 5.5)
## Boxplot
meta.rank <- sce.sub@meta.data %>%
  dplyr::select(subcelltype, `Glycolysis / Gluconeogenesis`:`Drug metabolism - other enzymes`) %>%
  data.table::melt() %>%
  dplyr::filter(variable %in% pathway_half)  


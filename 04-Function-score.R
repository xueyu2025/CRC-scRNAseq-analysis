rm(list=ls())
getwd()
library(fgsea)
library(magrittr)
library(ggplot2)
library(nlme)
library(Seurat)
library(qs)
library(readr)
library(tidyverse)
source("../00_script/data/scores.r")
HMPathways <- gmtPathways("../00_script/data/h.all.v2024.1.Hs.symbols.gmt")
sce = qread('adata_sample-2w.qs')
sce = subset(sce, Class != 'Border')
HMPathwaysscores <- score_cells(seur=sce, names=HMPathways, combine_genes='mean', groups=NULL, group_stat='mean', cells.use=NULL)
write_rds(HMPathwaysscores,path="scores/HMPathwaysScoreBycssmillieallcell.rds.gz",compress = "gz")
#HMPathwaysscores <- readr::read_rds("scores/HMPathwaysScoreBycssmillieallcell.rds.gz")
HMPathwaysscoresMatrix <- as.matrix(HMPathwaysscores)
HMPathwaysscoresMatrix <- t(HMPathwaysscoresMatrix )
HMPathwaysscoresMatrix <- apply(HMPathwaysscoresMatrix,2,function(x)signif(x,digits = 3))
colnames(HMPathwaysscoresMatrix) <- rownames(sce@meta.data)
HMPathwaysscoresMatrixSeurat <- CreateSeuratObject(counts = HMPathwaysscoresMatrix )
HMPathwaysscoresMatrixSeurat@meta.data <- sce@meta.data
sce.hm = HMPathwaysscoresMatrixSeurat
sce.hm@assays$Hallmark <- sce.hm@assays$RNA
glial_cells <- sce@meta.data[sce@meta.data$celltype == "Glial cell", ]
glial_classes <- glial_cells$Class
table(glial_classes) 
celltypes <- unique(sce.hm$celltype)  
checkClass = sort(unique(sce.hm$Class))

all_DEfeatures <- lapply(celltypes, function(ct) {
  sce.sub <- subset(sce.hm, celltype == ct)
  if (idx) {
    Idents(sce.sub) = sce.sub$Class
    DEfeatures = FindAllMarkers(sce.sub, pct = 0.1, logfc.threshold = 0.25, 
                                pseudocount.use = 0.1, only.pos = TRUE)
    
    # 检查是否有差异表达基因
    if (nrow(DEfeatures) == 0) {
      # 如果没有差异表达基因，跳过当前细胞类型
      return(NULL)
    }
    
    DEfeatures$celltype <- ct
    return(DEfeatures)  # 返回 DEfeatures 数据框
  } else {
    # 如果 Class 不匹配，返回 NULL 跳过当前迭代
    return(NULL)  # 或者返回 data.frame() 也可以
  }
})
all_DEfeatures <- all_DEfeatures[!sapply(all_DEfeatures, is.null)]
all_DEfeatures <- do.call(rbind, all_DEfeatures)
write.csv(all_DEfeatures, "scores/HM_score-TvsN.csv", row.names = FALSE)
final_DEfeatures = all_DEfeatures %>% filter(cluster == 'Tumor')
length(unique(final_DEfeatures$gene))
sce.hm = sce.hm %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000)
sce.hm = sce.hm %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 200)%>% 
  ScaleData(features = rownames(sce.hm))
hallmarkalias = read.table('../00_script/data/hallmark-alias.txt',header=TRUE,sep='\t')
hallmarkalias = hallmarkalias[hallmarkalias$hallmarktype !='Others',]
print(dim(hallmarkalias))
head(hallmarkalias)
hallmarkalias$pathway = gsub('_','-',hallmarkalias$pathway)
hallmarkalias$pathway[!hallmarkalias$pathway %in% data$pathway]
data = data %>% arrange(hallmarktype)
seq = unique(data$alias)
data$alias = factor(data$alias, levels = seq)

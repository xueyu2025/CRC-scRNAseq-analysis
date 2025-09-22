library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(qs)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(igraph)
lr_network <- readRDS('../00_database/Nichenet/lr_network_human_21122021.rds')
ligand_target_matrix <- readRDS("../00_database/Nichenet/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("../00_database/Nichenet/weighted_networks_nsga2r_final.rds")
lr_network <- lr_network %>% distinct(from, to)
head(lr_network)
seuratObj <- qread("../03_subset-bySample/Epithelial cell/Epi-Mye-PCA.qs")
seuratObj@meta.data %>% head()
nichenet_output_agnostic <- nichenet_seuratobj_aggregate(
  seurat_obj            = seuratObj,
  sender                = "Epi01",
  receiver              = "Macro_SPP1",
  condition_colname     = "Class",
  condition_oi          = "Tumor",
  condition_reference   = "Normal",
  expression_pct        = 0.05,
  ligand_target_matrix  = ligand_target_matrix,
  lr_network            = lr_network,
  weighted_networks     = weighted_networks
)
nichenet_output_agnostic$ligand_expression_dotplot
nichenet_output_agnostic$ligand_differential_expression_heatmap
nichenet_output_agnostic$ligand_target_heatmap
enrichGO(gene         = nichenet_output_agnostic$top_targets,
         OrgDb        = org.Mm.eg.db,
         keyType      = "SYMBOL",
         ont          = "BP")
edges <- nichenet_output_agnostic$ligand_target_df
g <- graph_from_data_frame(edges, directed = FALSE)
plot(g, vertex.label.cex = 0.8)
library(RColorBrewer)
fts <- nichenet_output_agnostic$top_targets[1:50]
fts <- fts[!is.na(fts)]
fts <- fts[fts %in% rownames(seuratObj)]
cls_levels <- unique(seuratObj@meta.data$Class)
n_cls <- length(cls_levels)
my_cols <- brewer.pal(min(n_cls, 2), "Set3")  
DotPlot(seuratObj %>% subset(idents = "Epi01"),
        features = nichenet_output_agnostic$top_targets[1:50] %>%
          rev(), split.by = "Class") + coord_flip()
nichenet_output_agnostic$ligand_activity_target_heatmap
ligands <- colnames(nichenet_output_agnostic$ligand_receptor_matrix)
Idents(seuratObj) <- "subcelltype"
epi_sub <- subset(
  seuratObj,
  idents = paste0("Epi", sprintf("%02d", 1:11))
)
avg_exp <- AverageExpression(
  epi_sub,
  assays        = "RNA",
  slot          = "data",
  features      = ligands,
  return.seurat = FALSE
)$RNA
mat_z <- t(scale(t(avg_exp)))
mat_z_reordered <- mat_z[
  nrow(mat_z):1,
  paste0("Epi", sprintf("%02d", 11:1))
]
col_name_epi01 <- "Epi01"
epi01_vals <- mat_z_reordered[, col_name_epi01]
row_order <- order(epi01_vals, decreasing = TRUE)
mat_z_reordered <- mat_z_reordered[row_order, ]
pheatmap(
  mat_z_reordered,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  color          = colorRampPalette(c("navy","white","firebrick3"))(50),
  main           = "Predictated Ligand expression on Epithelial cells",
  filename       = "Epi01-Macro_Spp1-LigandExpression.pdf",
  width          = 4,
  height         = 8,
  units          = "in"
)
pheatmap(
  mat_z_reordered,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  color          = colorRampPalette(c("navy","white","firebrick3"))(50),
  main           = "Predictated Ligand expression on Epithelial cells",
  filename       = "Epi01-Macro_Spp1-LigandExpression.jpg",
  width          = 4,
  height         = 8,
  units          = "in",
  dpi            = 300
)
nichenet_output_agnostic$ligand_target_heatmap
pdf("ligand_target_heatmap.pdf", width = 10, height = 8) 
print(nichenet_output_agnostic$ligand_target_heatmap)
dev.off()
jpeg("ligand_target_heatmap.jpeg",
     width = 10, height = 8, units = "in", res = 300)
print(nichenet_output_agnostic$ligand_target_heatmap)
dev.off()
nichenet_output_agnostic$ligand_receptor_matrix
nichenet_output_agnostic$ligand_receptor_heatmap
ligand_order <- rownames(mat_z_reordered)
head(ligand_order)
length(ligand_order)
lr_mat <- nichenet_output_agnostic$ligand_receptor_matrix
common_ligands <- intersect(ligand_order, colnames(lr_mat))
lr_mat2 <- lr_mat[, common_ligands, drop=FALSE]
others_lig <- setdiff(colnames(lr_mat), common_ligands)
lr_mat2 <- cbind(lr_mat2, lr_mat[, others_lig, drop=FALSE])
macro_vals   <- mat_z_rot["Macro_SPP1", ]
receptor_order <- colnames(mat_z_rot)[order(macro_vals, decreasing = TRUE)]
common_receptors <- intersect(receptor_order, rownames(lr_mat2))
lr_mat_final     <- lr_mat2[common_receptors, , drop=FALSE]
others_rec       <- setdiff(rownames(lr_mat2), common_receptors)
lr_mat_final     <- rbind(lr_mat_final, lr_mat2[others_rec, , drop=FALSE])
pheatmap(
  t(lr_mat_final),
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  color          = colorRampPalette(c("white", "#CC228D"))(50),
  border_color   = "grey80",
  main           = "Ligand–Receptor",
  filename       = "Epi01-Macro_Spp1-LigandReceptor.pdf",
  width          = 10,    # 英寸
  height         = 6.7,
  units          = "in"
)
pheatmap(
  t(lr_mat_final),
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  color          = colorRampPalette(c("white", "#CC228D"))(50),
  border_color   = "grey80",
  main           = "Ligand–Receptor",
  filename       = "Epi01-Macro_Spp1-LigandReceptor.jpg",
  width          = 10,
  height         = 6.7,
  units          = "in",
  dpi            = 300
)
receptors <- rownames(nichenet_output_agnostic$ligand_receptor_matrix)
Idents(seuratObj) <- "subcelltype"
clusters <- c("Macro_SPP1",
  "Macro_APOE", "Macro_F13A1", "Macro_MKI67","Mono_EREG", "Mono_FCN1", "Mono_HSPA1A",
  "DC_CCL19", "DC_CLEC9A", "DC_FCGR1A"
)
sub_obj <- subset(seuratObj, idents = clusters)
avg_exp <- AverageExpression(
  sub_obj,
  assays        = "RNA",
  slot          = "data",
  features      = receptors,
  return.seurat = FALSE
)$RNA
mat_z <- t(scale(t(avg_exp)))
new_cols <- c("Macro_SPP1", setdiff(clusters, "Macro_SPP1"))
mat_z <- mat_z[, new_cols]
mat_z <- mat_z[rev(rownames(mat_z)), ]
mat_z_rot <- t(mat_z)
mat_z_rot <- mat_z_rot[, rev(colnames(mat_z_rot))]
order_cols <- order(mat_z_rot["Macro_SPP1", ], decreasing = TRUE)
mat_z_rot <- mat_z_rot[, order_cols]
pheatmap(
  mat_z_rot,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,    # 这里显示 subtypes（Macro_SPP1 在首行）
  show_colnames  = TRUE,    # 这里显示 receptors，已按 Macro_SPP1 活性降序排序
  color          = colorRampPalette(c("navy","white","firebrick3"))(50),
  main           = "Predictated receptors expression on myeloid Cells",
  filename       = "Epi01-Macro_SPP1-ReceptorExpression.pdf",
  width          = 10,
  height         = 3,
  units          = "in"
)
pheatmap(
  mat_z_rot,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  color          = colorRampPalette(c("navy","white","firebrick3"))(50),
  main           = "Predictated receptors expression on myeloid Cells",
  filename       = "Epi01-Macro_SPP1-ReceptorExpression.jpg",
  width          = 10,
  height         = 3,
  units          = "in",
  dpi            = 300
)
nichenet_output <- nichenet_seuratobj_aggregate(
  seurat_obj            = seuratObj,
  sender                = "Macro_SPP1",
  receiver              = "Epi01",
  condition_colname     = "Class",
  condition_oi          = "Tumor",
  condition_reference   = "Normal",
  expression_pct        = 0.05,
  ligand_target_matrix  = ligand_target_matrix,
  lr_network            = lr_network,
  weighted_networks     = weighted_networks
)
top_ligands <- nichenet_output$top_ligands
selected_tf        <- "MXI1"
selected_ligands   <- c("LTB", "TFPI", "LGALS3")
selected_receptors <- c("TNFRSF1A","CD40","LTBR",
                        "SDC4","LRP1","F10","F3",
                        "ANXA2","PTPRK")
targets_provided <- c(
  "ACP1", "ENO1", "GAPDH", "PGK1", "LDHA", "SLC2A1",
  "SCD", "INSIG1", "HMGCS1", "HMGCS2", "PLIN2",
  "ASS1", "ODC1", "PYCR1", "ALDH18A1", "SLC1A5",
  "UQCRC2", "MDH2", "NQO1", "GPX2", "TXNIP", "BNIP3L",
  "FHL2", "GLO1", "GGH"
)
tf_lig_edges <- tibble(
  from = selected_tf,
  to   = selected_ligands
) %>%
  mutate(
    Weight = 1,
    Interaction = "TF->Ligand"
  )
lr_edges <- lr_network %>%
  filter(from %in% selected_ligands,
         to   %in% selected_receptors) %>%
  mutate(
    Weight = 1,
    Interaction = "Ligand->Receptor"
  ) %>%
  select(Source = from, Target = to, Interaction, Weight)
ligand_target_df <- ligand_target_matrix %>%
  as.data.frame() %>%
  t() %>%                                  
  as.data.frame() %>%
  rownames_to_column("ligand") %>%         
  pivot_longer(
    cols      = -ligand,
    names_to  = "target",
    values_to = "weight"
  ) %>%
  filter(
    ligand %in% selected_ligands,
    target %in% targets_provided,
    weight > 0
  )
rec_tgt_edges <- lr_edges %>%
  rename(ligand = Source, receptor = Target) %>%
  inner_join(ligand_target_df, by = "ligand") %>%
  transmute(
    Source      = receptor,
    Target      = target,
    Weight      = weight,
    Interaction = "Receptor->Target"
  )
edge_list <- bind_rows(
  tf_lig_edges  %>% select(Source = from, Target = to, Interaction, Weight),
  lr_edges,
  rec_tgt_edges
)
nichenet_output$ligand_activity_target_heatmap
p <- nichenet_output$ligand_activity_target_heatmap
ggsave(
  filename = "Macro_SPP1-Epi01-ligand_activity_target_heatmap.pdf",
  plot     = p,
  width    = 35,
  height   = 9,
  units    = "in"
)
ggsave(
  filename = "Macro_SPP1-Epi01-ligand_activity_target_heatmap.jpeg",
  plot     = p,
  width    = 35,
  height   = 9,
  units    = "in",
  dpi      = 300,
  device   = "jpeg"
)
top_ligs   <- nichenet_output$top_ligands
top_targs50 <- nichenet_output$top_targets[1:50]
top_targs50
mat50 <- nichenet_output$ligand_target_matrix[
  top_ligs,
  top_targs50,
  drop = FALSE
]
pheatmap(
  mat50,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,    # 这里显示 subtypes（Macro_SPP1 在首行）
  show_colnames  = TRUE,    # 这里显示 receptors，已按 Macro_SPP1 活性降序排序
  color          = colorRampPalette(c("whitesmoke", "purple"))(50),
  main           = "Top 50 Targets for Prioritized Ligands",
  filename       = "Macro_SPP1-Epi01-Top50-Targets.pdf",
  width          = 9.5,
  height         = 6,
  border_color   = "white",
  units          = "in"
)

pheatmap(
  mat50,
  cluster_rows   = FALSE,
  cluster_cols   = FALSE,
  show_rownames  = TRUE,
  show_colnames  = TRUE,
  color          = colorRampPalette(c("whitesmoke", "purple"))(50),
  main           = "Top 50 Targets for Prioritized Ligands",
  filename       = "Macro_SPP1-Epi01-Top50-Targets.jpg",
  width          = 9.5,
  height         = 6,
  units          = "in",
  border_color   = "white",
  dpi            = 300
)
row_order <- rev(rownames(mat100))
heatmap_reordered <- nichenet_output$ligand_receptor_heatmap +
  scale_y_discrete(limits = row_order)
print(heatmap_reordered)
ggsave("Macro_SPP1-Epi01-ligand_receptor_heatmap.pdf", heatmap_reordered,
       width = 10, height = 8, units = "in")
ggsave("Macro_SPP1-Epi01-ligand_receptor_heatmap.jpeg", heatmap_reordered,
       width = 10, height = 8, units = "in", dpi = 300)
DotPlot(seuratObj, features = nichenet_output$top_ligands %>% rev(), cols = "RdYlBu") +
  RotatedAxis()
clusters <- c(
  "Macro_SPP1","Macro_APOE", "Macro_F13A1", "Macro_MKI67",
  "Mono_EREG",  "Mono_FCN1",   "Mono_HSPA1A",
  "DC_CCL19",   "DC_CLEC9A",   "DC_FCGR1A"
)
Idents(seuratObj) <- "subcelltype"
seuratObj@active.ident <- factor(
  x      = seuratObj@active.ident,
  levels = clusters
)
subset_obj <- subset(seuratObj, idents = clusters)
p <- DotPlot(
      object   = subset_obj,
      features = rev(nichenet_output$top_ligands),
      cols     = "RdYlBu"
    ) +
    coord_flip() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    ggtitle("Top Ligands in Myeloid Subcelltypes")
ggsave(
  filename = "Macro_SPP1-Epi01-TopLigands-dotplot.pdf",
  plot     = p,
  width    = 6,
  height   = 8,
  units    = "in"
)
ggsave(
  filename = "Macro_SPP1-Epi01-TopLigands-dotplot.jpeg",
  plot     = p,
  width    = 6,
  height   = 8,
  units    = "in",
  dpi      = 300,
  device   = "jpeg"
)


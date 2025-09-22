library(devtools)
library(xCell)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)
library(dplyr)
library(data.table)
library(tidyverse)
library(ggpubr)
library(survival)
library(survminer)
library(gridExtra)
getwd()
setwd("C:\\Users\\ibmc\\Desktop\\CRC\\01_scVI_Integration_bySample\\KM生存曲线")
clinical_data <- fread(("./TCGA-COAD/data/TCGA-COAD.clinical.tsv.gz"), header=T, sep="\t")
expr_data <- fread(("./TCGA-COAD/data/TCGA-COAD.star_fpkm.tsv.gz"), header=T, sep="\t")
gene_map <- fread("./dataset/gencode.v36.annotation.gtf.gene.probemap", header=T, sep="\t")[, 1:2]
colnames(gene_map) <- c("Ensembl_ID", "Gene_Symbol")
survival_data <- fread('./TCGA-COAD/data/TCGA-COAD.survival.tsv.gz')
expr_annotated <- expr_data %>% 
  inner_join(gene_map, by=c("Ensembl_ID")) %>%
  distinct(Gene_Symbol, .keep_all=T) %>% 
  select(-Ensembl_ID) %>%                 
  column_to_rownames("Gene_Symbol")       
fwrite(expr_annotated, ("./TCGA-COAD/result/TCGA-COAD_fpkm_processed.tsv"), sep="\t", row.names=T)
markers <- c("EPCAM", "CDH1", "KRT18", "KRT8", "KRT19")
marker_expr <- expr_annotated[markers, , drop=FALSE]
log2_expr <- log2(marker_expr + 1)
epi_scores <- colMeans(log2_expr, na.rm = TRUE)
names(epi_scores) <- colnames(expr_annotated)
df_epi <- data.frame(
  sample    = names(epi_scores),
  epi.score = epi_scores,
  stringsAsFactors = FALSE
)
df_surv <- survival_data %>%
  inner_join(df_epi, by = "sample") %>%
  filter(substr(sample, 14, 15) == "01") %>%
  mutate(OS.time = OS.time / 365)
df_surv <- df_surv %>%
  mutate(EpiGroup = ifelse(epi.score >= median(epi.score), "High", "Low"),
         EpiGroup = factor(EpiGroup, levels = c("Low","High")))
fit_epi <- survfit(Surv(OS.time, OS) ~ EpiGroup, data = df_surv)
p_epi <- ggsurvplot(
  fit_epi,
  data       = df_surv,
  pval       = TRUE,
  risk.table = TRUE,
  xlab       = "Overall Survival (Years)",
  legend.title = "EPCAM\nExpression Group",
  legend.labs  = c("Low", "High"),
  palette      = c("#4DBBD5","#E64B35")
)
print(p_epi)
best_thresh <- cutoff::logrank(
  data    = df_surv,
  time    = "OS.time",    
  y       = "OS",         
  x       = "epi.score",  
  cut.num =1,           
  n.per   = 0.2,          
  y.per   = 0.1,          
  p.cut   = 1,            
  round   = 5            
)
best_thresh <- best_thresh[ order(best_thresh$pvalue), ][1, ]
cutoff_value <- best_thresh$cut1
print(best_thresh)  
df_surv <- df_surv %>%
  mutate(
    EpiGroup = ifelse(epi.score >= cutoff_value, "High", "Low"),
    EpiGroup = factor(EpiGroup, levels = c("Low","High"))
  )
fit_epi <- survfit(Surv(OS.time, OS) ~ EpiGroup, data = df_surv)
p_epi <- ggsurvplot(
  fit_epi,
  data         = df_surv,
  pval         = TRUE,
  pval.method  = TRUE,
  pval.size    = 6,
  pval.coord   = c(max(df_surv$OS.time)*0.6, 0.85),
  risk.table   = TRUE,
  surv.median.line = "hv",
  xlab         = "Overall Survival (Years)",
  legend.title = "EPCAM Expression\nGroup",
  legend.labs  = c("Low", "High"),
  palette      = c("#4DBBD5","#E64B35")
)
print(p_epi)
best_threshold_surv <- surv_cutpoint(df_surv,
                                     time = "OS.time",  
                                     event = "OS",      
                                     variables = "epi.score",  
                                     minprop = 0.05,     
                                     progressbar = TRUE)  
summary(best_threshold_surv)
best_threshold_data <- surv_categorize(best_threshold_surv)
surv_obj <- Surv(time = best_threshold_data$OS.time, event = best_threshold_data$OS)

surv_fit <- survfit(surv_obj ~ best_threshold_data$epi.score)

ggsurvplot(surv_fit, data = best_threshold_data, surv.median.line = "hv",
           pval = TRUE,                    
           xlab = "Time (days)", ylab = "Survival Probability",  
           legend.title = "",             
           #break.x.by = 1000,             
           color = "strata",              
           palette = c("#bc5148", "#3090a1"))  
epi_markers <- c("EPCAM","KRT19","KRT18")
my_markers  <- c("FCGR3A","CD14","C1QA")
epi_mat <- expr_annotated[epi_markers, , drop=FALSE]
my_mat  <- expr_annotated[my_markers,  , drop=FALSE]
epi_score <- colMeans(log2(epi_mat + 1), na.rm=TRUE)
my_score  <- colMeans(log2(my_mat  + 1), na.rm=TRUE)
df_scores <- data.frame(
  sample    = names(epi_score),
  epi.score = epi_score,
  my.score  = my_score,
  stringsAsFactors = FALSE
)
df_surv2 <- survival_data %>%
  inner_join(df_scores, by = "sample") %>%
  filter(substr(sample, 14, 15) == "01") %>%   
  mutate(OS.time = OS.time / 365)             
epi_med <- median(df_surv2$epi.score)
my_med  <- median(df_surv2$my.score)
df_surv2 <- df_surv2 %>%
  mutate(
    epi.high = epi.score >= epi_med,
    my.high  = my.score  >= my_med,
    CombinedGroup = ifelse(epi.high & my.high, "BothHigh", "Other"),
    CombinedGroup = factor(CombinedGroup, levels = c("Other","BothHigh"))
  )
fit_comb <- survfit(Surv(OS.time, OS) ~ CombinedGroup, data = df_surv2)
p_comb <- ggsurvplot(
  fit_comb,
  data         = df_surv2,
  pval         = TRUE,
  risk.table   = TRUE,
  xlab         = "Overall Survival (Years)",
  legend.title = "Epi & Myeloid\nBoth High",
  legend.labs  = c("Other","Both High"),
  palette      = c("#7570B3","#E7298A")
)
print(p_comb)
best_epi <- cutoff::logrank(
  data    = df_surv2,
  time    = "OS.time",
  y       = "OS",
  x       = "epi.score",
  cut.num = 1,
  n.per   = 0.2,
  y.per   = 0.1,
  p.cut   = 1,
  round   = 5
) %>% arrange(pvalue) %>% slice(1)
epi_thr <- best_epi$cut1
best_my <- cutoff::logrank(
  data    = df_surv2,
  time    = "OS.time",
  y       = "OS",
  x       = "my.score",
  cut.num = 1,
  n.per   = 0.2,
  y.per   = 0.1,
  p.cut   = 1,
  round   = 5
) %>% arrange(pvalue) %>% slice(1)
my_thr <- best_my$cut1
df_cl1 <- df_surv2 %>%
  mutate(
    epi.high = epi.score >= epi_thr,
    my.high  = my.score  >= my_thr,
    Combined = ifelse(epi.high & my.high, "BothHigh", "Other"),
    Combined = factor(Combined, levels=c("Other","BothHigh"))
  )
fit1 <- survfit(Surv(OS.time, OS) ~ Combined, data = df_cl1)
p1   <- ggsurvplot(
  fit1, data = df_cl1,
  pval       = TRUE,
  risk.table = TRUE,
  xlab       = "Overall Survival (Years)",
  legend.title  = "Epi & Myeloid\nBoth High (logrank)",
  legend.labs   = c("Other","BothHigh"),
  palette       = c("#7570B3","#E7298A")
)
print(p1)
cutp <- surv_cutpoint(
  data      = df_surv2,
  time      = "OS.time",
  event     = "OS",
  variables = c("epi.score","my.score"),
  minprop   = 0.1,
  progressbar = FALSE
)
df_cl2 <- surv_categorize(cutp)

df_cl2 <- as.data.frame(df_cl2)
class(df_cl2) <- "data.frame"
df_cl2 <- df_cl2 %>%
  rename(
    EpiHighGroup = `epi.score`,   
  )
df_cl2 <- df_cl2 %>%
  mutate(
    Combined = ifelse(`EpiHighGroup`=="high" & `MyHighGroup`=="high",
                      "BothHigh","Other"),
    Combined = factor(Combined, levels=c("Other","BothHigh"))
  )
fit2 <- survfit(Surv(OS.time, OS) ~ Combined, data = df_cl2)
p2   <- ggsurvplot(
  fit2, data = df_cl2,
  pval       = TRUE,
  risk.table = TRUE,
  xlab       = "Overall Survival (Years)",
  legend.title  = "Epi & Myeloid\nBoth High (cutpoint)",
  legend.labs   = c("Other","BothHigh"),
  palette       = c("#1B9E77","#D95F02")
)
print(p2)
df_q <- df_surv2 %>%
  mutate(
    epi.Q3 = epi.score >= quantile(epi.score, .75),
    my.Q3  = my.score  >= quantile(my.score,  .75),
    CombinedQ = case_when(
      epi.Q3 & my.Q3 ~ "BothHigh_Q3",
      !epi.Q3 & !my.Q3 ~ "BothLow_Q1",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(CombinedQ)) %>%         
  mutate(CombinedQ = factor(CombinedQ, levels=c("BothLow_Q1","BothHigh_Q3")))
fit_q <- survfit(Surv(OS.time, OS) ~ CombinedQ, data = df_q)
p_q   <- ggsurvplot(
  fit_q, data = df_q,
  pval       = TRUE,
  risk.table = TRUE,
  xlab       = "Overall Survival (Years)",
  legend.title = "Epi+Myeloid Extremes",
  legend.labs  = c("Both Low (Q1)", "Both High (Q3)")
)
print(p_q)
df_surv2 <- df_surv2 %>%
  mutate(risk = epi.score + my.score)
cut_risk <- surv_cutpoint(
  df_surv2,
  time      = "OS.time",
  event     = "OS",
  variables = "risk",
  minprop   = 0.1
)
df_risk <- surv_categorize(cut_risk)   # 会生成 risk Low/High
fit_r  <- survfit(Surv(OS.time, OS) ~ risk, data = df_risk)
p_r    <- ggsurvplot(
  fit_r, data = df_risk,
  pval       = TRUE,
  risk.table = TRUE,
  xlab       = "Overall Survival (Years)",
  legend.title = "TCGA-COAD Infiltration",
  palette      = c("red","blue")
)
print(p_r)
arranged <- arrange_ggsurvplots(
  list(p_r),
  print = FALSE,
  ncol  = 1,     # 一列
  nrow  = 1      # 一行
)
ggsave(
  filename = "TCGA_COAD_Infiltration_KM.pdf",
  plot     = arranged,
  width    = 6, height = 6, units = "in"
)
ggsave(
  filename = "TCGA_COAD_Infiltration_KM.jpg",
  plot     = arranged,
  width    = 6, height = 6, units = "in",
  dpi      = 300
)
library(survival)
library(survminer)
library(cutoff)
library(dplyr)
library(tibble)
epi_markers <- c("EPCAM","KRT19","KRT18")
my_markers  <- c("FCGR3A","CD14","C1QA")
epi_score <- expr_annotated[epi_markers, , drop=FALSE] %>% 
  { log2(. + 1) } %>% 
  colMeans(na.rm=TRUE)
my_score  <- expr_annotated[my_markers, , drop=FALSE] %>%
  { log2(. + 1) } %>% 
  colMeans(na.rm=TRUE)
df_base <- survival_data %>%
  mutate(sample = sample) %>%
  filter(substr(sample, 14, 15) == "01") %>%  # 仅保留 Primary Tumor
  mutate(OS.time = OS.time / 365)
df_epi <- df_base %>%
  inner_join(tibble(sample = names(epi_score), epi.score = epi_score),
             by = "sample")
df_my  <- df_base %>%
  inner_join(tibble(sample = names(my_score),  my.score  = my_score),
             by = "sample")
plot_all_methods <- function(df, score_col, label_prefix){
  df_med <- df %>%
    mutate(Group_med = factor(ifelse(.data[[score_col]] >= median(.data[[score_col]]),
                                     "High", "Low"),
                              levels=c("Low","High")))
  fit_med <- survfit(Surv(OS.time, OS) ~ Group_med, data=df_med)
  p_med   <- ggsurvplot(fit_med, df_med,
                        pval       = TRUE, risk.table=TRUE,
                        xlab       = "Overall Survival (Years)",
                        legend.title = paste0(label_prefix, "\nMedian split"),
                        legend.labs  = c("Low","High"),
                        palette      = c("blue","red"))
  best_lr <- cutoff::logrank(df, time="OS.time", y="OS",
                             x=score_col, cut.num=1,
                             n.per=0.2, y.per=0.1, p.cut=1, round=5) %>%
    arrange(pvalue) %>% slice(1)
  df_lr  <- df %>%
    mutate(Group_lr = factor(ifelse(.data[[score_col]] >= best_lr$cut1,
                                    "High","Low"),
                             levels=c("Low","High")))
  fit_lr <- survfit(Surv(OS.time, OS) ~ Group_lr, data=df_lr)

  p_lr <- ggsurvplot(
    fit_lr, df_lr,
    pval         = TRUE,
    risk.table   = TRUE,        # <-- 这里
    xlab         = "Overall Survival (Years)",
    legend.title = paste0(label_prefix, "\nLogRank cut"),
    legend.labs  = c("Low","High"),
    palette      = c("blue","red")
  )
  cutp  <- surv_cutpoint(df, time="OS.time", event="OS",
                         variables = score_col, minprop=0.1)
  df_sc <- as.data.frame( surv_categorize(cutp) )
  df_sc$Group_sc <- factor(
    ifelse(df_sc[[score_col]] == "high", "High","Low"),
    levels = c("Low","High")
  )
  fit_sc <- survfit(Surv(OS.time, OS) ~ Group_sc, data=df_sc)

  p_sc <- ggsurvplot(
    fit_sc, df_sc,
    pval         = TRUE,
    risk.table   = TRUE,        # <-- 这里
    xlab         = "Overall Survival (Years)",
    legend.title = paste0(label_prefix, "\nSurv_cutpoint"),
    legend.labs  = c("Low","High"),
    palette      = c("blue","red")
  )
  q75 <- quantile(df[[score_col]], 0.75, na.rm=TRUE)
  q25 <- quantile(df[[score_col]], 0.25, na.rm=TRUE)
  df_q <- df %>%
    mutate(Group_q = factor(case_when(
      .data[[score_col]] >= q75 ~ "High_Q3",
      .data[[score_col]] <= q25 ~ "Low_Q1",
      TRUE                      ~ NA_character_
    ), levels=c("Low_Q1","High_Q3"))) %>%
    filter(!is.na(Group_q))
  fit_q <- survfit(Surv(OS.time, OS) ~ Group_q, data=df_q)
  p_q   <- ggsurvplot(fit_q, df_q,
                      pval         = TRUE,
                      risk.table   = TRUE,
                      legend.title = paste0(label_prefix, "\nQ1 vs Q3"),
                      legend.labs  = c("Low (Q1)","High (Q3)"),
                      palette      = c("blue","red"))
  list(median = p_med,
       logrank = p_lr,
       cutpoint = p_sc,
       quartile = p_q)
}
plots_epi <- plot_all_methods(df_epi, "epi.score", "Epi markers")
print(plots_epi$median)
print(plots_epi$logrank)
print(plots_epi$cutpoint)
print(plots_epi$quartile)
library(survminer)
full_lr <- arrange_ggsurvplots(
  list(plots_epi$logrank),
  print = FALSE,
  ncol  = 1,
  nrow  = 1
)
ggsave("Epi_LogRankCut_full.pdf",
       plot   = full_lr,
       width  = 6, height = 8, units = "in")
ggsave("Epi_LogRankCut_full.jpeg",
       plot   = full_lr,
       width  = 6, height = 8, units = "in",
       dpi    = 300)
full_lr <- arrange_ggsurvplots(
  list(plots_epi$cutpoint),
  print = FALSE,
  ncol  = 1,
  nrow  = 1
)
ggsave("Epi_cutpoint_full.pdf",
       plot   = full_lr,
       width  = 6, height = 8, units = "in")
ggsave("Epi_cutpoint_full.jpeg",
       plot   = full_lr,
       width  = 6, height = 8, units = "in",
       dpi    = 300)
plots_my  <- plot_all_methods(df_my,  "my.score",  "Myeloid markers")
print(plots_my$median)
print(plots_my$logrank)
print(plots_my$cutpoint)
print(plots_my$quartile)
full_lr <- arrange_ggsurvplots(
  list(plots_my$median),
  print = FALSE,
  ncol  = 1,
  nrow  = 1
)
ggsave("Mye_median_full.pdf",
       plot   = full_lr,
       width  = 6, height = 8, units = "in")
ggsave("Mye_median_full.jpeg",
       plot   = full_lr,
       width  = 6, height = 8, units = "in",
       dpi    = 300)
full_lr <- arrange_ggsurvplots(
  list(plots_my$logrank),
  print = FALSE,
  ncol  = 1,
  nrow  = 1
)
ggsave("Mye_logrank_full.pdf",
       plot   = full_lr,
       width  = 6, height = 8, units = "in")
ggsave("Mye_logrank_full.jpeg",
       plot   = full_lr,
       width  = 6, height = 8, units = "in",
       dpi    = 300)
all_plots <- list(
  plots_epi$median, plots_epi$logrank,
  plots_epi$cutpoint, plots_epi$quartile,
  plots_my$median,  plots_my$logrank,
  plots_my$cutpoint, plots_my$quartile
)
arranged  <- arrange_ggsurvplots(all_plots, print=FALSE, ncol=2, nrow=4)
ggsave("Epi_vs_Myeloid_KM_all_methods.pdf", plot=arranged,
       width=12, height=16)





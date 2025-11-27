# ================================ /
# LCC Analysis Script /
# Author: Nils Briel, Dr., Dept. of Neurology, USZ, CH /
# Date: 2025-08-13 /
# ================================ /

rm(list = ls())

suppressPackageStartupMessages({
  library(readr)
  library(devtools)
  library(basicPlotteR)
  library(plyr)
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(RColorBrewer)
  library(ggpubr)
  library(ggsci)
  library(pheatmap)
  library(ComplexHeatmap)
  library(ggvenn)
  library(EnhancedVolcano)
  library(enrichplot)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
  library(ggupset)
  library(pathview)
  library(SummarizedExperiment)
  library(WGCNA)
  library(BioNERO)
  library(ggpubr)
  library(GO.db)
  library(tidyverse)
})

ves_col <- rev(RColorBrewer::brewer.pal(11, "RdBu"))
ves_col2 <- rev(RColorBrewer::brewer.pal(3, "RdBu"))[c(1, 3)]
wm_col <- rev(RColorBrewer::brewer.pal(11, "RdGy"))
wm_col2 <- rev(RColorBrewer::brewer.pal(3, "RdGy"))[c(1, 3)]

theme_set(theme_bw())

if(!dir.exists("../results/vessel_dda")){
  dir.create("../results/vessel_dda", recursive = TRUE)
}
if(!dir.exists("../results/wm_dda")){
  dir.create("../results/wm_dda", recursive = TRUE)
}
dir_vessel <- "../results/vessel_dda"
dir_wm <- "../results/wm_dda"

# ====================================================== /
# 1. Data Preparation ####
# ====================================================== /

Vessels <- read.csv("../data/Vessel DDA/Protein_Matrix_Cleanup_NoREV_NoMulti_NoUnmapped.csv",
                    stringsAsFactors = FALSE, row.names = 1, header = TRUE)
names(Vessels)[12:19] <- substring(names(Vessels)[12:19], 2)

WM <- read.csv("../data/White Matter DDA/Protein_Matrix_Cleanup_NoREV_NoMulti_NoUnmapped.csv",
               stringsAsFactors = FALSE, row.names = 1, header = TRUE)
names(WM)[12:19] <- substring(names(WM)[12:19], 2)

Annotation <- read_xlsx("../data/SampleAnnotation.xlsx")

Vessel_Matrix <- Vessels[, 12:19]
Vessel_Matrix <- Vessel_Matrix[, Annotation$ProbeLabel[1:8]]

WM_Matrix <- WM[, 12:19]
WM_Matrix <- WM_Matrix[, Annotation$ProbeLabel[9:16]]

names(Vessel_Matrix) <- Annotation$SampleName[match(names(Vessel_Matrix), Annotation$ProbeLabel)]
Vessel_Matrix <- Vessel_Matrix[-c(1147:1196),]

names(WM_Matrix) <- Annotation$SampleName[match(names(WM_Matrix), Annotation$ProbeLabel)]
WM_Matrix <- WM_Matrix[-c(1989:2068),]

# log2-transform data
Vessel_MatrixNorm <- Vessel_Matrix
Vessel_MatrixNorm[Vessel_MatrixNorm == 0] <- 1
Vessel_MatrixNorm <- log2(Vessel_MatrixNorm)

WM_MatrixNorm <- WM_Matrix
WM_MatrixNorm[WM_MatrixNorm == 0] <- 1
WM_MatrixNorm <- log2(WM_MatrixNorm)


# ====================================================== /
# 2 QC plots ####
# ====================================================== /

pdf("../results/vessel_dda/QC_sample_distribution_Vessels.pdf", height = 3.25, width = 4, pointsize = 8)
boxplot(Vessel_MatrixNorm, outline = FALSE,
        col = c(rep(ves_col2[1], 4), rep(ves_col2[2], 4)),
        ylab = expression("log"[2] * " protein intensity"), las = 2)
title("Vessels")
dev.off()

pdf("../results/wm_dda/QC_sample_distribution_WhiteMatter.pdf", height = 3.25, width = 4, pointsize = 8)
boxplot(WM_MatrixNorm, outline = FALSE,
        col = c(rep(wm_col2[1], 4), rep(wm_col2[2], 4)),
        ylab = expression("log"[2] * " protein intensity"), las = 2)
title("White Matter")
dev.off()

# ====================================================== /
## Venn diagrams — total peptides shared ####
# ====================================================== /

pdf("../results/venn_vessel_wm_total.pdf", height = 3.25, width = 4, pointsize = 8)
ggvenn(
  list(
    Vessel = rownames(Vessel_MatrixNorm),
    WM = rownames(WM_MatrixNorm)
  ),
  fill_color = c("#CD534CFF", "#868686FF"),
  stroke_size = 0.5, set_name_size = 4,
  auto_scale = TRUE
)
dev.off()

# Venn Diagram: by condition
ves_exp <- Vessel_MatrixNorm %>%
  mutate(peptide = rownames(.)) %>%
  tidyr::pivot_longer(
    cols = colnames(.)[-9],
    names_to = "SampleName",
    values_to = "ProtAb"
  ) %>%
  left_join(., Annotation[, c(1, 2)], by = "SampleName") %>%
  group_by(peptide, Condition) %>%
  dplyr::summarise(median_ProtAb = median(ProtAb)) %>%
  .[which(.$median_ProtAb > 0), ]

wm_exp <- WM_MatrixNorm %>%
  mutate(peptide = rownames(.)) %>%
  tidyr::pivot_longer(
    cols = colnames(.)[-9],
    names_to = "SampleName",
    values_to = "ProtAb"
  ) %>%
  left_join(., Annotation[, c(1, 2)], by = "SampleName") %>%
  group_by(peptide, Condition) %>%
  dplyr::summarise(median_ProtAb = median(ProtAb)) %>%
  .[which(.$median_ProtAb > 0), ]

x <- list(
  Vessel_H = subset(ves_exp, Condition == "Healthy")$peptide,
  Vessel_D = subset(ves_exp, Condition == "Diseased")$peptide,
  WM_D = subset(wm_exp, Condition == "Diseased")$peptide,
  WM_H = subset(wm_exp, Condition == "Healthy")$peptide
)

pdf("../results/venn_vessel_wm_subsets.pdf",
    height = 4, width = 5, pointsize = 8
)
ggvenn(
  x,
  fill_color = c(ves_col2[1], ves_col[9], wm_col2[2], wm_col2[1]),
  stroke_size = 0.5, set_name_size = 4
)
dev.off()

# ====================================================== /
# 3. All proteins heatmaps ####
# ====================================================== /

## Vessels: All proteins
ha_column <- HeatmapAnnotation(
  df = data.frame(
    Condition = Annotation$Condition[
      match(names(Vessel_MatrixNorm), Annotation[1:8, ]$SampleName)
    ],
    row.names = names(Vessel_MatrixNorm)
  ),
  col = list(Condition = c("Diseased" = ves_col2[2], "Healthy" = ves_col2[1])),
  simple_anno_size = unit(0.25, "cm"),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "plain"),
    labels_gp = gpar(fontsize = 8)
  )
)

zzz <- t(scale(t(Vessel_MatrixNorm)))
pdf("../results/vessel_dda/Heatmap_Vessels_AllProt.pdf",
    height = 3.25, width = 4
)
draw(
  Heatmap(
    zzz, col = ves_col, name = "Z-Score",
    top_annotation = ha_column,
    cluster_rows = TRUE, show_row_dend = TRUE,
    height = unit(2.25, "in"), width = unit(2, "in"),
    column_dend_height = unit(0.4, "cm"),
    show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 8),
    clustering_method_columns = "complete",
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      legend_width = unit(5, "cm"),
      title_position = "topcenter",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontface = "plain", fontsize = 8)
    )
  ),
  heatmap_legend_side = "left",
  annotation_legend_side = "right"
)
dev.off()

## WM: All proteins
ha_column <- HeatmapAnnotation(
  df = data.frame(
    Condition = Annotation$Condition[
      match(names(WM_MatrixNorm), Annotation[9:16, ]$SampleName)
    ],
    row.names = names(WM_MatrixNorm)
  ),
  col = list(Condition = c("Diseased" = wm_col2[2], "Healthy" = wm_col2[1])),
  simple_anno_size = unit(0.25, "cm"),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "plain"),
    labels_gp = gpar(fontsize = 8)
  )
)

zzz <- t(scale(t(WM_MatrixNorm)))
pdf("../results/wm_dda/Heatmap_WhiteMatter_AllProt.pdf",
    height = 3.25, width = 4
)
draw(
  Heatmap(
    zzz, col = wm_col, name = "Z-Score",
    top_annotation = ha_column,
    cluster_rows = TRUE, show_row_dend = TRUE,
    height = unit(2.25, "in"), width = unit(2, "in"),
    column_dend_height = unit(0.4, "cm"),
    show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 8),
    clustering_method_columns = "complete",
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      legend_width = unit(5, "cm"),
      title_position = "topcenter",
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontface = "plain", fontsize = 8)
    )
  ),
  heatmap_legend_side = "left",
  annotation_legend_side = "right"
)
dev.off()

# ====================================================== /
# 4. Differentially Expressed Proteins ####
# ====================================================== /

## --- Vessels Volcano --- ####
df5_Vessels <- Vessels[-c(1147:1196), ]
df5_Vessels[12:19] <- Vessel_MatrixNorm

for (i in 1:nrow(df5_Vessels)) {
  df5_Vessels$Max.fold.change[i] <- ifelse(
    df5_Vessels$Highest.mean.condition[i] == "Healthy Brain",
    -log2(df5_Vessels$Max.fold.change[i]),
    log2(df5_Vessels$Max.fold.change[i])
  )
}

pdf("../results/vessel_dda/Volcano_DiseasedHealthy_Vessels.pdf", width = 10, height = 5, pointsize = 8)
EnhancedVolcano(
  df5_Vessels[!is.infinite(df5_Vessels$Max.fold.change), ],
  lab = df5_Vessels$ProteinName[!is.infinite(df5_Vessels$Max.fold.change)],
  x = "Max.fold.change",
  y = "q.Value",
  xlim = c(-8, 15),
  ylim = c(0, 4.5),
  pCutoff = 5e-2,
  FCcutoff = 2,
  pointSize = 2.0,
  labSize = 5.0,
  title = "Volcano Plot: Vessels",
  subtitle = "Protein Abundance",
  caption = "Fold Change Cutoff: 2, q-value Cutoff: 5e-2",
  legendPosition = "right",
  legendLabSize = 14,
  col = c("grey30", "grey42", "royalblue", "red2"),
  colAlpha = 0.9,
  drawConnectors = TRUE,
  widthConnectors = 0.5
) + theme_bw()
dev.off()

up_Vessels   <- subset(df5_Vessels, q.Value < 0.05 & Max.fold.change >  1)
down_Vessels <- subset(df5_Vessels, q.Value < 0.05 & Max.fold.change < -1)

# Calculate ranks and order genes by decreasing rank score for GSEA. rank = sign(log2FC) * -log10(p.adj)
df5_Vessels$RankScore <- sign(df5_Vessels$Max.fold.change) * -log10(df5_Vessels$q.Value)
df5_Vessels <- df5_Vessels[order(df5_Vessels$RankScore, decreasing = TRUE), ]

write.csv(df5_Vessels,  "../results/vessel_dda/Volcano_Vessels_log2_rankscore_AllProt.csv")
write.csv(up_Vessels,   "../results/vessel_dda/Volcano_Vessels_log2_upregulated.csv")
write.csv(down_Vessels, "../results/vessel_dda/Volcano_Vessels_log2_downregulated.csv")

## --- Vessels DEP Heatmap --- ####
Vessel_DEP_matrix <- rbind(up_Vessels[, 12:19], down_Vessels[, 12:19])
names(Vessel_DEP_matrix) <- Annotation$SampleName[match(names(Vessel_DEP_matrix), Annotation$ProbeLabel)]

ha_column <- HeatmapAnnotation(
  df  = data.frame(Condition = Annotation$Condition[match(names(Vessel_MatrixNorm), Annotation[1:8, ]$SampleName)],
                   row.names = names(Vessel_DEP_matrix)),
  col = list(Condition = c("Diseased" = ves_col2[2], "Healthy" = ves_col2[1]))
)

zzz <- t(scale(t(Vessel_DEP_matrix)))

pdf("../results/vessel_dda/Heatmap_Vessels_DEPonly.pdf", height = 3.25, width = 4)
draw(
  Heatmap(
    zzz, col = ves_col, name = "Z-Score", top_annotation = ha_column,
    cluster_rows = TRUE, show_row_dend = TRUE,
    height = unit(2.25, "in"), width = unit(2, "in"),
    column_dend_height = unit(0.4, "cm"), show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 8),
    clustering_method_columns = "complete"
  ),
  heatmap_legend_side = "left", annotation_legend_side = "right"
)
dev.off()

# ====================================================== /
## --- WM Volcano --- ####
df5_WM <- WM[-c(1989:2068), ]
df5_WM[12:19] <- WM_MatrixNorm

for (i in 1:nrow(df5_WM)) {
  df5_WM$Max.fold.change[i] <- ifelse(
    df5_WM$Highest.mean.condition[i] == "Healthy Brain",
    -log2(df5_WM$Max.fold.change[i]),
    log2(df5_WM$Max.fold.change[i])
  )
}

pdf("../results/wm_dda/Volcano_DiseasedHealthy_WhiteMatter.pdf", width = 10, height = 5, pointsize = 8)
EnhancedVolcano(
  df5_WM[!is.infinite(df5_WM$Max.fold.change), ],
  lab = df5_WM$ProteinName[!is.infinite(df5_WM$Max.fold.change)],
  x = "Max.fold.change",
  y = "q.Value",
  xlim = c(-10, 22),
  ylim = c(0, 4.5),
  pCutoff = 5e-2,
  FCcutoff = 2,
  pointSize = 2.0,
  labSize = 5.0,
  title = "Volcano Plot: White Matter",
  subtitle = "Protein Abundance",
  caption = "Fold Change Cutoff: 2, q-value Cutoff: 5e-2",
  legendPosition = "right",
  legendLabSize = 14,
  col = c("grey30", "grey42", "royalblue", "red2"),
  colAlpha = 0.9,
  drawConnectors = TRUE,
  widthConnectors = 0.5
) + theme_bw()
dev.off()

up_WM   <- subset(df5_WM, q.Value < 0.05 & Max.fold.change >  1)
down_WM <- subset(df5_WM, q.Value < 0.05 & Max.fold.change < -1)

df5_WM$RankScore <- sign(df5_WM$Max.fold.change) * -log10(df5_WM$q.Value)
df5_WM <- df5_WM[order(df5_WM$RankScore, decreasing = TRUE), ]

write.csv(df5_WM,  "../results/wm_dda/Volcano_WhiteMatter_log2_rankscore_AllProt.csv")
write.csv(up_WM,   "../results/wm_dda/Volcano_WhiteMatter_log2_upregulated.csv")
write.csv(down_WM, "../results/wm_dda/Volcano_WhiteMatter_log2_downregulated.csv")

## --- WM DEP Heatmap --- ####
WM_DEP_matrix <- rbind(up_WM[, 12:19], down_WM[, 12:19])
names(WM_DEP_matrix) <- Annotation$SampleName[match(names(WM_DEP_matrix), Annotation$ProbeLabel)]

ha_column <- HeatmapAnnotation(
  df  = data.frame(Condition = Annotation$Condition[match(names(WM_MatrixNorm), Annotation[9:16, ]$SampleName)],
                   row.names = names(WM_DEP_matrix)),
  col = list(Condition = c("Diseased" = wm_col2[2], "Healthy" = wm_col2[1]))
)

zzz <- t(scale(t(WM_DEP_matrix)))

pdf("../results/wm_dda/Heatmap_WhiteMatter_DEPonly.pdf", height = 3.25, width = 4)
draw(
  Heatmap(
    zzz, col = wm_col, name = "Z-Score", top_annotation = ha_column,
    cluster_rows = TRUE, show_row_dend = TRUE,
    height = unit(2.25, "in"), width = unit(2, "in"),
    column_dend_height = unit(0.4, "cm"), show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 8),
    clustering_method_columns = "complete"
  ),
  heatmap_legend_side = "left", annotation_legend_side = "right"
)
dev.off()

# ====================================================== /
# 5. Pathway Enrichment Analyses ####
# ====================================================== /

## 5.1 Vessels GO Enrichment --- ####
Alldata <- read.csv("../results/vessel_dda/Volcano_Vessels_log2_rankscore_AllProt.csv", stringsAsFactors = FALSE)
names(Alldata)[1] <- "UniprotID"
up <- read.csv("../results/vessel_dda/Volcano_Vessels_log2_upregulated.csv", stringsAsFactors = FALSE, row.names = 1)

FC <- Alldata$Max.fold.change
names(FC) <- Alldata$UniprotID

conflicted::conflicts_prefer(clusterProfiler::simplify)
egoUpreg <- enrichGO(
  gene          = row.names(up),
  universe      = Alldata$UniprotID,
  keyType       = "UNIPROT",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
ego2Upreg <- simplify(egoUpreg)
write.csv(data.frame(ego2Upreg), "../results/vessel_dda/enrichGO_Upreg.csv", row.names = FALSE)

pdf("../results/vessel_dda/enrichGO_Upreg_Barplot.pdf", height = 4, width = 6)
barplot(ego2Upreg, x = "GeneRatio", font.size = 8,
        showCategory = length(ego2Upreg@result$ID)) + ylab("Gene ratio")
dev.off()

## reduce redundant pathways
ego2Upreg <- simplify(egoUpreg)
## save output
cluster_summary <- data.frame(ego2Upreg)
write.csv(cluster_summary,
          file = "../results/vessel_dda/enrichGO_Upreg.csv",
          row.names = F)

## barplot of all 18 significant pathways
pdf("../results/vessel_dda/enrichGO_Upreg_Barplot.pdf",
    height = 4, width = 6)
barplot(ego2Upreg, x = "GeneRatio", font.size = 8,
        showCategory = length(ego2Upreg@result[["ID"]])) +
  ylab(label = "Gene ratio")
dev.off()

## Heatplot
pdf("../results/vessel_dda/EnrichGO_Upreg_Heatplot.pdf",
    height = 5, width = 15)
heatplot(ego2Upreg, foldChange = FC, showCategory = length(ego2Upreg@result[["ID"]]), ) + scale_color_manual(values = "grey4") +
  scale_fill_gradient(name = expression(log[2](FC)), low = ves_col[6], high = ves_col[11]) + theme_test() + coord_equal() +
  theme(axis.title = element_text(size = 6), legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()

## Lollipop plot
p <- ggplot(cluster_summary, aes(x = Count, 
                                 y = tidytext::reorder_within(Description, Count, list()), 
                                 size = -log(p.adjust))) + 
  geom_vline(xintercept = 0, alpha = 0.6) +
  geom_segment(aes(y = tidytext::reorder_within(Description, Count,  list()), 
                   xend=0, xend=Count), size = 0.5, color="grey4") +
  geom_point(stat = "identity", shape = 21, 
             fill = "blue3") + 
  theme_bw() + 
  labs(title = "",
       y = "",
       x = "Gene Count, DEPs") + 
  theme(strip.text.y = element_text(angle = 0))
pdf("../results/vessel_dda/EnrichGO_Upreg_lollipop.pdf", width = 8, height = 4)
plot(p)
dev.off()

## Emapplot
pdf("../results/vessel_dda/enrichGO_Upreg_Emapplot.pdf",
    height = 8, width = 8)
emapplot(pairwise_termsim(x = ego2Upreg))
dev.off()

pdf("../results/vessel_dda/enrichGO_Upreg_Cnetplot.pdf",
    height = 16, width = 8)
p1 <- cnetplot(ego2Upreg, foldChange = FC, node_label="category", layout = "kk",
               showCategory = length(ego2Upreg@result[["ID"]]), 
               color_category =ves_col2[2],
               color_gene = ves_col2[1]) + 
  theme(legend.position = "none")
p2 <- cnetplot(ego2Upreg, foldChange = FC, node_label="gene", layout = "kk",
               showCategory = length(ego2Upreg@result[["ID"]]), 
               color_category =ves_col2[2],
               color_gene = ves_col2[1])
cowplot::plot_grid(p1, p2, nrow=2, labels=LETTERS[1:2], rel_heights = c(1.2,0.8))
dev.off()

## 5.2 Vessels KEGG Enrichment --- ####
if (!dir.exists("../results/kegg")) dir.create("../results/kegg")
sig <- subset(Alldata, q.Value < 0.1)

KEGGUpreg <- enrichKEGG(
  gene          = rownames(up),
  universe      = Alldata$UniprotID,
  keyType       = "uniprot",
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)
KEGGUpreg <- setReadable(KEGGUpreg, OrgDb = org.Hs.eg.db, keyType="UNIPROT")
## save output
cluster_summary <- data.frame(KEGGUpreg)
write.csv(cluster_summary,
          file = "../results/vessel_dda/enrichKEGG_Upreg.csv",
          row.names = F)

# Generdate pathview plots
setwd("../results/kegg/")
ComplCoag <- pathview(gene.data  = FC,
                      pathway.id = KEGGUpreg$ID[1],
                      gene.idtype = "ACCNUM",
                      species    = "hsa",
                      bins = 25,
                      limit      = list(gene=15, cpd=1))

Lupus <- pathview(gene.data  = FC,
                  pathway.id = KEGGUpreg$ID[2],
                  gene.idtype = "ACCNUM",
                  species    = "hsa",
                  bins = 25,
                  limit      = list(gene=15, cpd=1))

skeleton <- pathview(gene.data  = FC,
                     pathway.id = KEGGUpreg$ID[3],
                     gene.idtype = "ACCNUM",
                     species    = "hsa",
                     bins = 25,
                     limit      = list(gene=15, cpd=1))

ferrop <- pathview(gene.data  = FC,
                   pathway.id = KEGGUpreg$ID[4],
                   gene.idtype = "ACCNUM",
                   species    = "hsa",
                   bins = 25,
                   limit      = list(gene=15, cpd=1))

pw_cancer <- pathview(gene.data  = FC,
                      pathway.id = "hsa05200",
                      gene.idtype = "ACCNUM",
                      species    = "hsa",
                      bins = 25,
                      limit      = list(gene=15, cpd=1))

mapk_cancer <- pathview(gene.data  = FC,
                        pathway.id = "hsa04010",
                        gene.idtype = "ACCNUM",
                        species    = "hsa",
                        bins = 25,
                        limit      = list(gene=15, cpd=1))

egr_inh_resis <- pathview(gene.data  = FC,
                          pathway.id = "hsa01521+5156",
                          gene.idtype = "ACCNUM",
                          species    = "hsa",
                          bins = 25,
                          limit      = list(gene=15, cpd=1))

pld_cancer <- pathview(gene.data  = FC,
                       pathway.id = "hsa04072+5154",
                       gene.idtype = "ACCNUM",
                       species    = "hsa",
                       bins = 25,
                       limit      = list(gene=15, cpd=1))

setwd("../../src/")

# ====================================================== /
## 5.3 WM GO Enrichment --- ####
Alldata <- read.csv("../results/wm_dda/Volcano_WhiteMatter_log2_rankscore_AllProt.csv", stringsAsFactors = FALSE)
names(Alldata)[1] <- "UniprotID"
up <- read.csv("../results/wm_dda/Volcano_WhiteMatter_log2_upregulated.csv", stringsAsFactors = FALSE, row.names = 1)

FC <- up$Max.fold.change
names(FC) <- rownames(up)

egoUpreg <- enrichGO(
  gene          = row.names(up),
  universe      = Alldata$UniprotID,
  keyType       = "UNIPROT",
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

## reduce redundant pathways: 10 pathways
ego2Upreg <- simplify(egoUpreg)
## save output
cluster_summary <- data.frame(ego2Upreg)
write.csv(cluster_summary,
          file = "../results/wm_dda/enrichGO_Upreg.csv",
          row.names = F)

## barplot of all 10 significant pathways
pdf("../results/wm_dda/WM_enrichGO_Upreg_Barplot.pdf",
    height = 3.25, width = 6)
barplot(ego2Upreg, x = "GeneRatio", font.size = 8,
        showCategory = length(ego2Upreg@result[["ID"]])) +
  ylab(label = "Gene ratio")
dev.off()

## Heatplot
pdf("../results/wm_dda/WM_EnrichGO_Upreg_Heatplot.pdf",
    height = 5, width = 10)
heatplot(ego2Upreg, foldChange = FC, showCategory = length(ego2Upreg@result[["ID"]]), ) + scale_color_manual(values = "grey4") +
  scale_fill_gradient(name = expression(log[2](FC)), low = wm_col[6], high = wm_col[11]) + theme_test() + coord_equal() +
  theme(axis.title = element_text(size = 6), legend.title = element_text(size = 8),
        legend.text = element_text(size = 8), axis.text.x = element_text(angle = 45, hjust = 1)) 
dev.off()

## Lollipop plot
p <- ggplot(cluster_summary, aes(x = Count, 
                                 y = tidytext::reorder_within(Description, Count, list()), 
                                 size = -log(p.adjust))) + 
  geom_vline(xintercept = 0, alpha = 0.6) +
  geom_segment(aes(y = tidytext::reorder_within(Description, Count,  list()), 
                   xend=0, xend=Count), size = 0.5, color="grey4") +
  geom_point(stat = "identity", shape = 21, 
             fill = "black") + 
  theme_bw() + 
  labs(title = "",
       y = "",
       x = "Gene Count, DEPs") + 
  theme(strip.text.y = element_text(angle = 0))
pdf("../results/wm_dda/EnrichGO_Upreg_lollipop.pdf", width = 8, height = 4)
plot(p)
dev.off()

## Emapplot
pdf("../results/wm_dda/WM_enrichGO_Upreg_Emapplot.pdf",
    height = 6, width = 6)
emapplot(pairwise_termsim(x = ego2Upreg))
dev.off()

pdf("../results/wm_dda/WM_enrichGO_Upreg_Cnetplot.pdf",
    height = 10, width = 8)
p1 <- cnetplot(ego2Upreg, foldChange = FC, node_label="category", layout = "kk",
               showCategory = length(ego2Upreg@result[["ID"]]),
               color_category = "grey42") + 
  scale_color_gradient2(name=expression(log[2](FC)), low = wm_col[6], high = wm_col[11]) + theme(legend.position = "none")
p2 <- cnetplot(ego2Upreg, foldChange = FC, node_label="gene", layout = "kk",
               showCategory = length(ego2Upreg@result[["ID"]]),
               color_category = "grey42") + 
  scale_color_gradient2(name=expression(log[2](FC)), low = wm_col[6], high = wm_col[11])
cowplot::plot_grid(p1, p2, nrow=2, labels=LETTERS[1:2], rel_heights = c(1.2,0.8))
dev.off()

## 5.4 WM KEGG Enrichment --- ####
KEGGUpreg <- enrichKEGG(gene = row.names(up),
                        universe      = Alldata$UniprotID,
                        keyType = "uniprot",
                        organism         = "hsa",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05)
KEGGUpreg <- setReadable(KEGGUpreg, OrgDb = org.Hs.eg.db, keyType="UNIPROT")

# ====================================================== /
# 6 WGCNA ####
# ====================================================== /
## 6.1 Data Preparation ####
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(WGCNA)
  library(BioNERO)
  library(ggpubr)
})
multicoreParam <-BiocParallel::SnowParam(workers = 4)
conflicted::conflict_prefer("filter", "dplyr")
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflicts_prefer(dplyr::mutate)

Metadata <- Annotation %>% as.data.frame() %>%
  `rownames<-`(paste0(.$SampleName,"_",ifelse(.$TissueType == "Vessels", "V", "WM"))) %>%
  mutate(ID = rownames(.),
         col = c(rep(ves_col2[1],4),rep(ves_col2[2],4),c(rep(wm_col2[1],4),rep(wm_col2[2],4)))
  )

AlldataV <- read.csv("../results/vessel_dda/Volcano_Vessels_log2_rankscore_AllProt.csv", stringsAsFactors = F,
                     header = T) %>% .[,c(1:12)]
names(AlldataV)[1] <- "UniprotID"
rowannot <- join(Alldata[1:12],AlldataV, by = "ProteinName") 

data <- join(Vessels[,c(1,12:19)], WM[,c(1,12:19)], by = "ProteinName")  %>%
  `rownames<-`(.$ProteinName) %>% select(-"ProteinName") %>% as.matrix()

#change col and rownames
indices <- match(colnames(data), Metadata$ProbeLabel)
colnames(data) <- Metadata[indices, "ID"]
data <- data[,rownames(Metadata)]

rowannot <- data.frame(row.names = rownames(data), ProteinName = rownames(data))

zma.se <- SummarizedExperiment(assays = SimpleList(counts = as.matrix(data)),
                               colData = DataFrame(Metadata),
                               rowData = DataFrame(rowannot),
                               checkDimnames = T)

if(!file.exists("../results/final_exp.Rds")){
  exp_filt <- BioNERO::replace_na(zma.se)
  exp_filt <- remove_nonexp(exp_filt, method = "median", min_exp = 10)
  exp_filt <- filter_by_variance(exp_filt, n = 2000)
  exp_filt <- ZKfiltering(exp_filt, cor_method = "pearson")
  exp_filt <- PC_correction(exp_filt)
  final_exp <- exp_preprocess(
    zma.se, min_exp = 10, variance_filter = TRUE, n = 2000
  )
  saveRDS(final_exp,"../results/final_exp.Rds")
}else{
  final_exp <- readRDS("../results/final_exp.Rds")
}

## 6.2 Soft Threshold ####
cols <- data.frame(row.names =  rownames(final_exp@colData), col = final_exp$col)
pdf("../results/exp_cor_heatmap.pdf", width = 6, height = 5.5)
plot_heatmap(final_exp, coldata_cols = "Condition", type = "samplecor", show_rownames = F)
dev.off()
pdf("../results/exp_dis_heatmap.pdf", width = 6, height = 7)
plot_heatmap(final_exp[1:50,], type = "expr", coldata_cols = c("Condition"), show_rownames = T, cor_method = "pearson", show_colnames = FALSE)
dev.off()
pdf("../results/exp_pca.pdf", width = 6, height = 5)
plot_PCA(final_exp, metadata_cols = "Condition") + theme_test() + coord_equal(1.5) + scale_color_manual(values = rev(ves_col2))
dev.off()

# Infer a gene coexpression network (GCN) using exp2gcn()
sft <- SFT_fit(final_exp, net_type = "signed hybrid", cor_method = "pearson")
power <- sft$power
pdf("../results/wgcna_sft.pdf", width = 8, height = 3.5)
sft$plot
dev.off()

## 6.3 Dendro and Module Characterization ####
if(!file.exists("../results/net_wgcna.Rds")){
  net <- exp2gcn(
    final_exp,
    net_type = "signed hybrid",
    SFTpower = power,
    cor_method = "pearson"
  )
  saveRDS(net,"../results/net_wgcna.Rds")
}else{
  net <- readRDS("../results/net_wgcna.Rds")
}

pdf("../results/module_dendro.pdf", width = 6, height = 4)
plot(plot_dendro_and_colors(net))
dev.off()
pdf("../results/module_eigengene_network.pdf", width = 6, height = 5)
plot(plot_eigengene_network(net))
dev.off()
pdf("../results/module_ngenes_per_module.pdf", width = 6, height = 4)
plot(plot_ngenes_per_module(net))
dev.off()

module_stability(final_exp, net, nRuns = 5)
ggsave("../results/module_stability.pdf", width = 6, height = 4)
final_exp$Condition_TissueType <- paste0(final_exp$Condition,":",final_exp$TissueType)
MEtrait <- module_trait_cor(exp = final_exp,
                            MEs = net$MEs,
                            cor_method = "spearman",
                            metadata_cols = c("Condition","TissueType","Condition_TissueType"))
pdf("../results/module_trait_correlation_heatmap.pdf", width = 8, height = 6)
plot_module_trait_cor(MEtrait, palette = "RdBu")
dev.off()

p1 <- plot_expression_profile(
  exp = final_exp, 
  metadata_cols = "Condition_TissueType",
  net = net, 
  plot_module = TRUE, 
  modulename = "salmon"
) + scale_fill_manual(values = c("red3","salmon", "grey3", "grey42")) + theme_test() +
  theme(legend.position = "none",  axis.text.x = element_blank())

p2 <- plot_expression_profile(
  exp = final_exp, 
  metadata_cols = "Condition_TissueType",
  net = net, 
  plot_module = TRUE, 
  modulename = "cyan"
) + scale_fill_manual(values = c("cyan4","cyan", "grey3", "grey42")) + theme_test() +
  theme(legend.position = "none", axis.text.x = element_blank())

p3 <- plot_expression_profile(
  exp = final_exp, 
  metadata_cols = "Condition_TissueType",
  net = net, 
  plot_module = TRUE, 
  modulename = "black"
) + scale_fill_manual(values = c("grey3","grey23", "grey", "grey72")) + theme_test()+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

pdf("../results/module_trait_correlation.pdf", width = 6, height = 8)
cowplot::plot_grid(plotlist = list(p1, p2, p3), nrow=3, labels=LETTERS[1:3], rel_heights = c(0.9,0.9,1.2))
dev.off()

## 6.4 Hub genes ####
hubs <- get_hubs_gcn(final_exp, net)

modules <- unique(hubs$Module)

for (module in modules) {
  edges <- get_edge_list(net, module = module)
  edges_filtered <- get_edge_list(
    net, module = module,
    filter = TRUE, method = "pvalue",
    pvalue_cutoff = 0.001,
    nSamples = ncol(final_exp)
  )
  
  pdf(paste0("../results/network_", module, ".pdf"), width = 7, height = 5)
  plot(
    plot_gcn(
      edgelist_gcn = edges_filtered,
      net = net,
      color_by = "module",
      hubs = hubs,
      top_n_hubs = 10,
      show_labels = "allhubs"
    )
  )
  dev.off()
}

## 6.5 Enrichment Analysis ####

go_annotations <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(Alldata$ProteinName),
  columns = c("GO"),
  keytype = "SYMBOL"
)

go_an <- AnnotationDbi::select(
  GO.db,
  keys = go_annotations$GO,
  columns = c("DEFINITION", "GOID", "ONTOLOGY", "TERM"),
  keytype = "GOID"
) %>% drop_na()

if (!file.exists("../results/interpro_enrichment.Rds")) {
  interpro_enrichment <- module_enrichment(
    net = net,
    background_genes = Alldata$ProteinName,
    p = 0.1,
    min_setsize = 5,
    annotation = go_annotations
  )
  saveRDS(interpro_enrichment, "../results/interpro_enrichment.Rds")
} else {
  interpro_enrichment <- readRDS("../results/interpro_enrichment.Rds")
}

# KEGG enrichment per module loop
dat_ls2 <- list()
for (i in unique(hubs$Module)) {
  genes <- subset(net$genes_and_modules, Modules == i)$Genes
  KEGGUpreg <- enrichKEGG(
    gene          = subset(Alldata, ProteinName %in% genes)$UniprotID,
    universe      = Alldata$UniprotID,
    keyType       = "uniprot",
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.1,
    qvalueCutoff  = 1,
    minGSSize     = 5
  )
  
  if (!is.null(KEGGUpreg)) {
    KEGGUpreg <- setReadable(KEGGUpreg, OrgDb = org.Hs.eg.db, keyType = "UNIPROT")
    
    df <- subset(KEGGUpreg@result, pvalue < 0.1)
    if (nrow(df) > 0) {
      df$module <- i
      dat_ls2[[i]] <- df
    }
  }
}

if (length(dat_ls2) > 0) {
  dat2 <- do.call(rbind, dat_ls2) %>%
    group_by(module) %>%
    slice_min(order_by = p.adjust, n = 5, with_ties = FALSE) %>%
    ungroup()
  
  p <- ggplot(dat2, aes(
    x = Count,
    y = tidytext::reorder_within(Description, Count, list()),
    size = -log(p.adjust),
    fill = module
  )) +
    geom_vline(xintercept = 0, alpha = 0.6) +
    geom_segment(aes(
      y = tidytext::reorder_within(Description, Count, list()),
      xend = 0, yend = tidytext::reorder_within(Description, Count, list())
    ), size = 0.5, color = "grey4") +
    geom_point(stat = "identity", shape = 21) +
    facet_grid(module ~ ., drop = TRUE, space = "free", scales = "free") +
    theme_bw() +
    labs(title = "", y = "", x = "Gene Count, WGCNA Modules") +
    scale_fill_manual(values = unique(dat2$module)) +
    theme(strip.text.y = element_text(angle = 0))
  
  pdf("../results/module_kegg.pdf", width = 7, height = 8.5)
  plot(p)
  dev.off()
  
  # Pathview diagrams — example for selected pathways
  setwd("../results/kegg/")
  for (pw in c("hsa03008", "hsa03015", "hsa04370", "hsa01521+7422")) {
    pathview(
      gene.data = FC,
      pathway.id = pw,
      gene.idtype = "ACCNUM",
      species = "hsa",
      bins = 25,
      limit = list(gene = 15, cpd = 1)
    )
  }
}

# ====================================================== /
# 7. RTK Pathway Correlation with WGCNA Modules ####
# ====================================================== /

## 7.1 Define RTK pathway genes (from hsa01521+7422 modified pathway)
rtk_genes <- c(
  # Growth factors and receptors
  "VEGFA", "VEGFB", "VEGFC", "VEGFD",  # VEGF family
  "PGF",  # PDGF, EGF, BDNF represented by placental growth factor
  "FLT1", "KDR", "FLT4",  # VEGFR1/2/3 (RTK receptors)
  "PDGFRA", "PDGFRB",  # PDGFR
  "EGFR",  # EGF receptor
  
  # Downstream signaling molecules
  "SRC", "LCK",  # Src family kinases (SRC, Sck)
  "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG",  # PI3K
  "GRB2",  # GRB2 adaptor
  "RASA1", "RASGRF1", "RASGRF2",  # RasGRP
  "HRAS", "KRAS", "NRAS",  # Ras
  "RAF1", "BRAF", "ARAF",  # Raf
  "MAP2K1", "MAP2K2",  # MEK1/2
  "MAPK1", "MAPK3",  # ERK1/2
  
  # p38 MAPK pathway
  "MAP3K3", "MAP3K6",  # MKK3/6
  "MAPK14",  # p38 MAPK
  "MAPKAPK2",  # MAPK (downstream of p38)
  "HSPB1",  # HSP27
  
  # PKC and DAG pathway
  "PLCG1", "PLCG2",  # PLCγ1
  "PRKCA", "PRKCB", "PRKCD", "PRKCE",  # PKC
  
  # Cytoskeletal and focal adhesion
  "PTK2", "PTK2B",  # PTK2 (FAK)
  "PXN",  # Paxillin
  "MAPT",  # Tau
  "STMN1",  # STMN1
  
  # Ca2+ and NO pathway
  "NOS3",  # NOS (endothelial NOS)
  "PTGS1", "PTGS2",  # PG (prostaglandin synthase/COX)
  
  # GNA12 and actin organization
  "GNA12", "GNA13"  # GNA12
)

## 7.2 Extract gene symbols and create mapping
extract_gene_symbol <- function(desc) { 
  if(is.na(desc)) return(NA)
  if(!grepl("GN=", desc)) return(NA)
  gene <- sub(".*GN=([^ ]+).*", "\\1", desc)
  return(gene)
}

# Create mappings for Vessels & WM ds
vessels_mapping <- data.frame(
  UniprotID = rownames(Vessels),
  ProteinName = Vessels$ProteinName,
  GeneSymbol = sapply(Vessels$Description, extract_gene_symbol),
  Description = Vessels$Description,
  stringsAsFactors = FALSE
)

wm_mapping <- data.frame(
  UniprotID = rownames(WM),
  ProteinName = WM$ProteinName,
  GeneSymbol = sapply(WM$Description, extract_gene_symbol),
  Description = WM$Description,
  stringsAsFactors = FALSE
)

# Find RTK genes in Vessel % wm ds
rtk_in_vessel <- vessels_mapping %>%
  filter(GeneSymbol %in% rtk_genes) %>%
  mutate(Dataset = "Vessel")

rtk_in_wm <- wm_mapping %>%
  filter(GeneSymbol %in% rtk_genes) %>%
  mutate(Dataset = "WM")

rtk_in_data <- bind_rows(
  rtk_in_vessel %>% select(UniprotID, ProteinName, GeneSymbol, Dataset),
  rtk_in_wm %>% select(UniprotID, ProteinName, GeneSymbol, Dataset)
) %>%
  group_by(UniprotID, GeneSymbol, ProteinName) %>%
  summarise(
    Dataset = paste(unique(Dataset), collapse = ", "),
    .groups = "drop"
  )

# Get full data for RTK genes from both dss with metadata
rtk_full_vessel <- Vessels %>%
  rownames_to_column("UniprotID") %>%
  filter(UniprotID %in% rtk_in_data$UniprotID) %>%
  left_join(rtk_in_data %>% select(UniprotID, GeneSymbol), by = "UniprotID") %>%
  mutate(Dataset = "Vessel")

rtk_full_wm <- WM %>%
  rownames_to_column("UniprotID") %>%
  filter(UniprotID %in% rtk_in_data$UniprotID) %>%
  left_join(rtk_in_data %>% select(UniprotID, GeneSymbol), by = "UniprotID") %>%
  mutate(Dataset = "WM")

# Save 
write.csv(rtk_in_data, "../results/rtk_pathway_genes_in_dataset.csv", row.names = FALSE)
write.csv(rtk_full_vessel, "../results/rtk_pathway_genes_vessel_full.csv", row.names = FALSE)
write.csv(rtk_full_wm, "../results/rtk_pathway_genes_wm_full.csv", row.names = FALSE)

## 7.3 Map RTK genes to WGCNA data
wgcna_mapping <- data.frame(
  ProteinName = rownames(final_exp),
  stringsAsFactors = FALSE
) %>%
  left_join(vessels_mapping %>% select(ProteinName, GeneSymbol, UniprotID), 
            by = "ProteinName")
rtk_in_wgcna <- wgcna_mapping %>%
  filter(GeneSymbol %in% rtk_genes) %>%
  filter(!is.na(GeneSymbol))

rtk_expr <- assay(final_exp)[rtk_in_wgcna$ProteinName, , drop = FALSE]
rownames(rtk_expr) <- rtk_in_wgcna$GeneSymbol

## 7.4 Correlations between RTK genes and module eigengenes
module_eigengenes <- net$MEs

if(nrow(rtk_expr) == 0) {
  cat("ERROR: No RTK genes found in WGCNA data. Cannot proceed with correlation analysis.n")
} else {
  cat("Calculating correlations between", nrow(rtk_expr), "RTK genes and", ncol(module_eigengenes), "modules...nn")
  rtk_module_cor <- list()
  rtk_module_pval <- list()
  
  for (gene in rownames(rtk_expr)) {
    gene_expr <- rtk_expr[gene, ]
    
    for (module in colnames(module_eigengenes)) {
      cor_test <- cor.test(gene_expr, module_eigengenes[, module], method = "spearman")
      
      if (is.null(rtk_module_cor[[gene]])) {
        rtk_module_cor[[gene]] <- list()
        rtk_module_pval[[gene]] <- list()
      }
      
      rtk_module_cor[[gene]][[module]] <- cor_test$estimate
      rtk_module_pval[[gene]][[module]] <- cor_test$p.value
    }
  }
  
  cor_matrix <- do.call(rbind, lapply(rtk_module_cor, function(x) unlist(x)))
  pval_matrix <- do.call(rbind, lapply(rtk_module_pval, function(x) unlist(x)))
  
  # Clean up
  colnames(cor_matrix) <- gsub("[.]rho$", "", colnames(cor_matrix))
  colnames(pval_matrix) <- gsub("[.]rho$", "", colnames(pval_matrix))
  
  # Adjust p-v
  pval_adj_matrix <- matrix(
    p.adjust(as.vector(pval_matrix), method = "BH"),
    nrow = nrow(pval_matrix),
    ncol = ncol(pval_matrix),
    dimnames = dimnames(pval_matrix)
  )
  
  # significance annotation
  sig_matrix <- ifelse(pval_adj_matrix < 0.001, "***",
                      ifelse(pval_adj_matrix < 0.01, "**",
                            ifelse(pval_adj_matrix < 0.05, "*", "")))
  
  cat("Correlation matrix dimensions:", dim(cor_matrix), "n")
  cat("Significant correlations (FDR < 0.05):", sum(pval_adj_matrix < 0.05), "n")
  cat("Significant correlations (FDR < 0.01):", sum(pval_adj_matrix < 0.01), "nn")
}

## 7.5 Visualization
if (exists("cor_matrix") && nrow(cor_matrix) > 0) {
  library(ComplexHeatmap)
  library(circlize)
  
  col_fun <- colorRamp2(
    c(-1, -0.5, 0, 0.5, 1),
    c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
  )
  
  # Annotation for modules (columns)
  module_ids <- colnames(cor_matrix)
  module_colors_named <- setNames(gsub("ME", "", module_ids), module_ids)
  
  ha_top <- HeatmapAnnotation(
    Module = factor(module_ids, levels = names(module_colors_named)),
    col = list(Module = module_colors_named),
    annotation_name_side = "left",
    simple_anno_size = unit(0.3, "cm"),
    annotation_name_gp = gpar(fontsize = 8)
  )
  
  ht <- Heatmap(
    cor_matrix,
    name = "Spearman ρ",
    col = col_fun,
    top_annotation = ha_top,
    
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 9),
    show_row_dend = TRUE,
    row_dend_width = unit(1, "cm"),
    
    column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    show_column_dend = TRUE,
    column_dend_height = unit(1, "cm"),
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (sig_matrix[i, j] != "") {
        grid.text(sig_matrix[i, j], x, y, gp = gpar(fontsize = 8, col = "black"))
      }
    },
    
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(4, "cm"),
      direction = "vertical"
    ),
    
    width = unit(0.5 * ncol(cor_matrix), "cm"),
    height = unit(0.5 * nrow(cor_matrix), "cm")
  )
  
  pdf("../results/rtk_module_correlation_heatmap.pdf", width = 10, height = 8)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  cat("Generated: rtk_module_correlation_heatmap.pdf")
}

# ====================================================== /
# Final session info
# ====================================================== /
sessioninfo::session_info()
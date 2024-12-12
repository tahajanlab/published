library(CellChat)
library(patchwork)
library(NMF)
library(circlize)
library(ComplexHeatmap)
library(Seurat)
library(Matrix)
library(magrittr)
library(tidyverse)
library(SingleCellExperiment)
library(remotes)
library(future)
library(ggpubr)
library(ggalluvial)
library(scales)
library(pals)
library(viridis)

# Load Seurat Object, Visualize Cluster Groups
so <- readRDS("analyzed_utricle_mes_se_p4p6.rds")

# Relevel in descending order
clust.desc.order <- Idents(so) %>% table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% select(-Freq)
colnames(clust.desc.order) <- "clust.desc.order"
clust.desc.order
levels(so) <- clust.desc.order %>% as.matrix()
levels(so) # Confirm reorder

ggplot2::update_geom_defaults("point",list(size = 5, stroke = 0.5, alpha = 0.7))

# custom.hex.colors <- c("Mesenchymal Cells"="#84BF96","Type I Hair Cells"="#A6CEE3","Transitional Epithelial Cells"="#F3E587",
#                        "Supporting Cells"="#DE9E83","Type II Hair Cells"="#B15928",
#                        "Glia"="#FBB268","Roof Cells (CNMD+)"="#6DBD57","Pericytes"="#3385BB",
#                        "Schwann Cells"="#9D7BBA","Endothelial Cells"="#F57C7C",
#                        "Macrophages"="#7F9D56","Melanocytes"="#FE8D1A")

#FINAL
# custom.hex.colors <- c("Mesenchymal Cells"="#3366CC","Type I Hair Cells"="#DC3912","Transitional Epithelial Cells"="#FF9900",
#                        "Supporting Cells"="#109618","Type II Hair Cells"="#990099",
#                        "Glia"="#0099C6","Roof Cells (CNMD+)"="#DD4477","Pericytes"="#64C3B7",
#                        "Schwann Cells"="#B82E2E","Endothelial Cells"="#66AA00",
#                        "Macrophages"="#F46FC7","Melanocytes"="#94594D")


# # Cell Chat Colors
# 
# custom.hex.colors <- c("Mesenchymal Cells"="#E41A1C","Type I Hair Cells"="#377EB8","Transitional Epithelial Cells"="#4DAF4A",
#                        "Supporting Cells"="#984EA3","Type II Hair Cells"="#F29403",
#                        "Glia"="#F781BF","Roof Cells (CNMD+)"="#BC9DCC","Pericytes"="#A65628",
#                        "Schwann Cells"="#54B0E4","Endothelial Cells"="#222F75",
#                        "Macrophages"="#1B9E77","Melanocytes"="#B2DF8A")

# Final
custom.hex.colors <- c("Mesenchymal Cells"="#3366CC","Type I Hair Cells"="#DC3912","Transitional Epithelial Cells"="#FF9900",
                       "Supporting Cells"="#109618","Type II Hair Cells"="#990099",
                       "Glia"="#0099C6","Roof Cells (CNMD+)"="#DD4477","Pericytes"="#64C3B7",
                       "Schwann Cells"="#B82E2E","Endothelial Cells"="#66AA00",
                       "Macrophages"="#F46FC7","Melanocytes"="#94594D")


marker.list <- c("Pecam1", "Stard13", "Rgs5", "Mpz", "Iba1", "Dcn", "Coch", "Pmel")

so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Heatmap (adj p-value < 0.01)

sig.markers <- so.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01) %>%
  top_n(n = 25, wt = avg_log2FC)
write_csv(sig.markers, "fig1-heatmap-top25-markers.csv")
p2 <- DoHeatmap(so, features = sig.markers$gene, group.colors = custom.hex.colors) + scale_fill_viridis(na.value = "white")

pdf("heatmap cc.pdf", width = 15, height = 20)
p2
dev.off()



sig.markers <- so.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01)
write_csv(sig.markers, "fig1-heatmap-all-markers.csv")
p1 <- DoHeatmap(so, features = sig.markers$gene, group.colors = custom.hex.colors) + scale_fill_viridis(na.value = "white")

pdf("heatmap-fine cc.pdf", width = 15, height = 20)
p1
dev.off()

DoHeatmap(so, features = sig.markers$gene, group.colors = custom.hex.colors) + scale_fill_viridis(na.value = "white")

#DoHeatmap(so, features = sig.markers$gene, group.colors = custom.hex.colors) + scale_fill_manual(values=parula(), na.value = "white")
# FeaturePlot(so, features = marker.list, order = TRUE) & NoAxes()
# VlnPlot(so, features = marker.list)
# DoHeatmap(so, features = sig.markers$gene) + NoLegend()
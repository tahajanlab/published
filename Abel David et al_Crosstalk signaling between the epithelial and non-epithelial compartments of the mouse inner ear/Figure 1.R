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

# Load Seurat Object, Visualize Cluster Groups
so <- readRDS("analyzed_utricle_mes_se_p4p6.rds")

# Relevel in descending order
clust.desc.order <- Idents(so) %>% table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% select(-Freq)
colnames(clust.desc.order) <- "clust.desc.order"
clust.desc.order
levels(so) <- clust.desc.order %>% as.matrix()
levels(so) # Confirm reorder

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

# Frequency Table

ggplot2::update_geom_defaults("point",list(size = 5, stroke = 0.5, alpha = 0.7))

# custom.hex.colors <- c("Mesenchymal Cells"="#84BF96","Type I Hair Cells"="#A6CEE3","Transitional Epithelial Cells"="#F3E587",
#                        "Supporting Cells"="#DE9E83","Type II Hair Cells"="#B15928",
#                        "Glia"="#FBB268","Roof Cells (CNMD+)"="#6DBD57","Pericytes"="#3385BB",
#                        "Schwann Cells"="#9D7BBA","Endothelial Cells"="#F57C7C",
#                        "Macrophages"="#7F9D56","Melanocytes"="#FE8D1A")


# custom.hex.colors <- c("Mesenchymal Cells"="#3366CC","Type I Hair Cells"="#DC3912","Transitional Epithelial Cells"="#FF9900",
#                        "Supporting Cells"="#109618","Type II Hair Cells"="#990099",
#                        "Glia"="#0099C6","Roof Cells (CNMD+)"="#DD4477","Pericytes"="#66AA00",
#                        "Schwann Cells"="#B82E2E","Endothelial Cells"="#316395",
#                        "Macrophages"="#994499","Melanocytes"="#22AA99")

# Final
custom.hex.colors <- c("Mesenchymal Cells"="#3366CC","Type I Hair Cells"="#DC3912","Transitional Epithelial Cells"="#FF9900",
                       "Supporting Cells"="#109618","Type II Hair Cells"="#990099",
                       "Glia"="#0099C6","Roof Cells (CNMD+)"="#DD4477","Pericytes"="#64C3B7",
                       "Schwann Cells"="#B82E2E","Endothelial Cells"="#66AA00",
                       "Macrophages"="#F46FC7","Melanocytes"="#94594D")

# CellChat

# custom.hex.colors <- c("Mesenchymal Cells"="#E41A1C","Type I Hair Cells"="#377EB8","Transitional Epithelial Cells"="#4DAF4A",
#                        "Supporting Cells"="#984EA3","Type II Hair Cells"="#F29403",
#                        "Glia"="#F781BF","Roof Cells (CNMD+)"="#BC9DCC","Pericytes"="#A65628",
#                        "Schwann Cells"="#54B0E4","Endothelial Cells"="#222F75",
#                        "Macrophages"="#1B9E77","Melanocytes"="#B2DF8A")


pdf("frequency-labels.pdf", width = 8, height = 6)
freq <- Idents(so) %>% table() %>% as.data.frame() %>% arrange(desc(Freq))
colnames(freq) <- c("Cell Type", "Freq")
str(freq)
freq %>% 
  mutate(labels = Freq) %>%
  ggplot(aes(x="", y = Freq, fill=`Cell Type`, colour = `Cell Type`)) + geom_col() + coord_polar("y", start=0) + blank_theme + labs_pubr() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="bottom", axis.text.x=element_blank()) + geom_text(aes(label = labels), colour = "black", position = position_stack(vjust = 0.5), show.legend = FALSE)  + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()

pdf("frequency-nolabels.pdf", width = 8, height = 6)
freq <- Idents(so) %>% table() %>% as.data.frame() %>% arrange(desc(Freq))
colnames(freq) <- c("Cell Type", "Freq")
str(freq)
freq %>% 
  mutate(labels = Freq) %>%
  ggplot(aes(x="", y = Freq, fill=`Cell Type`, colour = `Cell Type`)) + geom_col() + coord_polar("y", start=0) + blank_theme + labs_pubr() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="bottom", axis.text.x=element_blank()) +
  scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
# + geom_text(aes(label = labels), colour = "black", position = position_stack(vjust = 0.5), show.legend = FALSE) 
dev.off()

so$ident <- Idents(so)
# UMAP
pdf("umap-labels.pdf", width = 2, height = 2)
plot <- DimPlot(so, reduction = "umap", group.by = "ident")
LabelClusters(plot = plot, id = "ident") + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()

pdf("umap-nolabels.tiff", width = 2, height = 2)
plot + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()

pdf("umap-labels-bigger.pdf", width = 3, height = 3)
plot <- DimPlot(so, reduction = "umap")
LabelClusters(plot = plot, id = "ident") + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()

pdf("umap-nolabels-bigger.tiff", width = 3, height = 3)
plot + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()

pdf("umap-labels-biggest.pdf", width = 4, height = 4)
plot <- DimPlot(so, reduction = "umap")
LabelClusters(plot = plot, id = "ident") + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()

pdf("umap-nolabels-biggest.pdf", width = 4, height = 4)
plot + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()


# ---
marker.list <- c("Pecam1", "Stard13", "Rgs5", "Mpz", "Iba1", "Dcn", "Coch", "Pmel")

so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Heatmap (adj p-value < 0.01)
sig.markers <- so.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01)
DoHeatmap(so, features = sig.markers$gene) + NoLegend() 


# scale_colour_manual(values = custom.hex.colors)

# DoHeatmap(so, features = marker.list 


FeaturePlot(so, features = marker.list, order = TRUE) & NoAxes()
VlnPlot(so, features = marker.list)
# Mesenchymal Markers from Wilkerson 2021 
features.use <-  c("Wif1", "Hpgd", "Gpc3", "Coch", "Fbxo2", "Car3", "Cldn11", "Crym", "Pam")
FeaturePlot(object = so, features = features.use) & NoAxes()

# Cell Cycle Markers
features.use <-  c("Mki67", "Top2a")
FeaturePlot(object = so, features = features.use) & NoAxes()

# Dark Cell Markers
features.use <-  c("Kcne1", "Lrp2")
FeaturePlot(object = so, features = features.use) & NoAxes()


DoHeatmap(so, features = sig.markers$gene) + NoLegend()


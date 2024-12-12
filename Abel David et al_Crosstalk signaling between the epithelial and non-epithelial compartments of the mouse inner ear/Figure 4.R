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

# Subset
so <- subset(so, idents = c("Mesenchymal Cells"))
so <- FindNeighbors(so, dims = 1:20)
so <- FindClusters(so, resolution = 0.6) # the resolution = 2
Idents(so)
so <- RunUMAP(so, dims = 1:20)
DimPlot(so, reduction = "umap")

# Frequency Table

ggplot2::update_geom_defaults("point",list(size = 5, stroke = 0.5, alpha = 0.7))

so <- RenameIdents(so, "0" = "S0", "1" = "S1", "2" = "S2")

#custom.hex.colors <- c("S0"="#3366CC","S1"="#DC3912","S2"="#FF9900")
custom.hex.colors <- c("S0"="#5DDAFF","S1"="#DC3912","S2"="#FF9900")

pdf("mes-sub frequency-labels.pdf", width = 8, height = 6)
freq <- Idents(so) %>% table() %>% as.data.frame() %>% arrange(desc(Freq))
colnames(freq) <- c("Cell Type", "Freq")
str(freq)
freq %>% 
  mutate(labels = Freq) %>%
  ggplot(aes(x="", y = Freq, fill=`Cell Type`, colour = `Cell Type`)) + geom_col() + coord_polar("y", start=0) + blank_theme + labs_pubr() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="bottom", axis.text.x=element_blank()) + geom_text(aes(label = labels), colour = "black", position = position_stack(vjust = 0.5), show.legend = FALSE) + scale_fill_manual(values = custom.hex.colors, aesthetics = c("colour", "fill"))
dev.off()

pdf("mes-sub frequency-nolabels.pdf", width = 8, height = 6)
freq <- Idents(so) %>% table() %>% as.data.frame() %>% arrange(desc(Freq))
colnames(freq) <- c("Cell Type", "Freq")
str(freq)
freq %>% 
  mutate(labels = Freq) %>%
  ggplot(aes(x="", y = Freq, fill=`Cell Type`, colour = `Cell Type`)) + geom_col() + coord_polar("y", start=0) + blank_theme + labs_pubr() + theme(plot.title = element_text(hjust = 0.5)) + theme(legend.position="bottom", axis.text.x=element_blank()) +
  scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
# + geom_text(aes(label = labels), colour = "black", position = position_stack(vjust = 0.5), show.legend = FALSE) 
dev.off()

so$cluster <- Idents(so)

# UMAP
pdf("mes-sub sub-umap-labels.pdf", width = 2, height = 2)
plot <- DimPlot(so, reduction = "umap", group.by = "cluster")
LabelClusters(plot = plot, id = "cluster") + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()

pdf("mes-sub sub-umap-nolabels.pdf", width = 2, height = 2)
plot + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()

pdf("mes-sub sub-umap-labels-bigger.pdf", width = 3, height = 3)
plot <- DimPlot(so, reduction = "umap", group.by = "cluster")
LabelClusters(plot = plot, id = "cluster") + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()

pdf("mes-sub sub-umap-nolabels-bigger.pdf", width = 3, height = 3)
plot + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()

# Violin Plots

seurat.clusters <- Idents(so)

clust <- "S0"
celltype.bool <- case_when(
  seurat.clusters != clust ~ paste0("Not ", clust),
  seurat.clusters == clust ~ clust
)
so <- AddMetaData(so, celltype.bool, col.name = "celltype.bool")
so$celltype.bool <- so$celltype.bool %>% fct_relevel(clust)
Idents(so) <- so$celltype.bool
celltype.de <- FindMarkers(so, ident.1 = clust, ident.2 = paste0("Not ", clust))
celltype.de.sig <- celltype.de %>% filter(p_val_adj  < 0.01) %>% top_n(n = 5, wt = avg_log2FC)
print(rownames(celltype.de.sig))
pdf(paste0("mes-sub-", clust, ".pdf"), width = 12, height = 2.5)
# VlnPlot(so, features = rownames(celltype.de.sig), cols = c("#3366CC", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
#VlnPlot(so, features = rownames(celltype.de.sig), cols = c("#5DDAFF", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
select.features <- c("Wif1", "Rbp1", "Gpr37", "Abhd2", "Plekhb1")
VlnPlot(so, features = select.features, cols = c("#5DDAFF", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
dev.off()


clust <- "S1"
celltype.bool <- case_when(
  seurat.clusters != clust ~ paste0("Not ", clust),
  seurat.clusters == clust ~ clust
)
so <- AddMetaData(so, celltype.bool, col.name = "celltype.bool")
so$celltype.bool <- so$celltype.bool %>% fct_relevel(clust)
Idents(so) <- so$celltype.bool
celltype.de <- FindMarkers(so, ident.1 = clust, ident.2 = paste0("Not ", clust))
celltype.de.sig <- celltype.de %>% filter(p_val_adj  < 0.01) %>% top_n(n = 5, wt = avg_log2FC)
print(rownames(celltype.de.sig))
pdf(paste0("mes-sub-", clust, ".pdf"), width = 12, height = 2.5)
#VlnPlot(so, features = rownames(celltype.de.sig), cols = c("#DC3912", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
select.features <- c("Col8a1", "Gja1", "Frzb", "Ntrk2", "Rbp4")
VlnPlot(so, features = select.features, cols = c("#DC3912", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
dev.off()

clust <- "S2"
celltype.bool <- case_when(
  seurat.clusters != clust ~ paste0("Not ", clust),
  seurat.clusters == clust ~ clust
)
so <- AddMetaData(so, celltype.bool, col.name = "celltype.bool")
so$celltype.bool <- so$celltype.bool %>% fct_relevel(clust)
Idents(so) <- so$celltype.bool
celltype.de <- FindMarkers(so, ident.1 = clust, ident.2 = paste0("Not ", clust))
celltype.de.sig <- celltype.de %>% filter(p_val_adj  < 0.01) %>% top_n(n = 5, wt = avg_log2FC)
print(rownames(celltype.de.sig))
pdf(paste0("mes-sub-", clust, ".pdf"), width = 12, height = 2.5)
#VlnPlot(so, features = rownames(celltype.de.sig), cols = c("#FF9900", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
select.features <- c("Gpc3", "Dmp1", "Cpm", "Gcg", "Lum")
VlnPlot(so, features = select.features, cols = c("#FF9900", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
dev.off()

# HeatMap
Idents(so) <- seurat.clusters
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


library(viridis)
# Heatmap (adj p-value < 0.01)

sig.markers <- so.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01) %>%
  top_n(n = 25, wt = avg_log2FC)
p2 <- DoHeatmap(so, features = sig.markers$gene, group.colors = custom.hex.colors) + scale_fill_viridis(na.value = "white")

pdf("mes-sub-heatmap.pdf", width = 15, height = 15)
p2
dev.off()

sig.markers <- so.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01)
p1 <- DoHeatmap(so, features = sig.markers$gene, group.colors = custom.hex.colors) + scale_fill_viridis(na.value = "white")

pdf("mes-sub-heatmap-fine.pdf", width = 15, height = 20)
p1
dev.off()

write_csv(sig.markers, "mes-subcluster.de.genes.csv")

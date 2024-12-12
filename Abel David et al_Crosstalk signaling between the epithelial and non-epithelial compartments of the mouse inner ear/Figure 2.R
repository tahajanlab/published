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
#BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(ggpubr)

# Load Seurat Object, Visualize Cluster Groups
so <- readRDS("analyzed_utricle_mes_se_p4p6.rds")

# Relevel in descending order
clust.desc.order <- Idents(so) %>% table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% select(-Freq)
colnames(clust.desc.order) <- "clust.desc.order"
clust.desc.order
levels(so) <- clust.desc.order %>% as.matrix()
levels(so) # Confirm reorder

# Frequency Table

ggplot2::update_geom_defaults("point",list(size = 5, stroke = 0.5, alpha = 0.7))


features.use <-  c("Epcam")
plot <- FeaturePlot(object = so, features = features.use)

pdf("epcam-labels-bigger.pdf", width = 3, height = 3)
plot + ggtitle("") & scale_color_continuous(low = "#FFFDB2", high = "red")
dev.off()

pdf("epcam-nolabels-bigger.pdf", width = 3, height = 3)
plot + NoLegend() + NoAxes() + ggtitle("") & scale_color_continuous(low = "#FFFDB2", high = "red")
dev.off()

#add cell type metadata
cluster.names <- Idents(so)
so <- subset(so, idents = c("Type I Hair Cells", "Type II Hair Cells", "Supporting Cells", "Transitional Epithelial Cells", "Roof Cells (CNMD+)", "Mesenchymal Cells"))
epithelium <- case_when(
  cluster.names %in% c("Type I Hair Cells", "Type II Hair Cells", "Supporting Cells", "Transitional Epithelial Cells", "Roof Cells (CNMD+)") ~ "Epithelial",
  cluster.names %in% c("Mesenchymal Cells") ~ "Non-Epithelial"
)

so <- AddMetaData(so, epithelium, col.name = "epithelium")
Idents(so) <-so$epithelium

epithelium.de <- FindMarkers(so, ident.1 = "Epithelial", ident.2 = "Non-Epithelial")
epithelium.de.sig <- epithelium.de

write.csv(epithelium.de, "rev.epithelium.de.genes.csv")

write.csv(epithelium.de.sig, "rev.epithelium.de.sig.genes.csv")

epithelium.de.labs<- epithelium.de %>% filter(p_val_adj  < 0.01) %>% top_n(n = 200, wt = abs(avg_log2FC))

pdf("rev epithelial volcano.pdf", width = 20, height = 20)
EnhancedVolcano(epithelium.de.sig, 
                rownames(epithelium.de.sig),
                selectLab = rownames(epithelium.de.labs), 
                x ="avg_log2FC", 
                y ="p_val_adj", 
                pCutoff = 0.01,
                FCcutoff = 2.0,
                col=c('black', 'black', 'black', 'red')) + theme_pubr() + NoLegend()
dev.off()

epithelium.de.labs<- epithelium.de %>% filter(p_val_adj  < 0.01) %>% top_n(n = 200, wt = abs(avg_log2FC))

pdf("rev epithelial volcano select.pdf", width = 5, height = 5)
EnhancedVolcano(epithelium.de.sig, 
                rownames(epithelium.de.sig),
                selectLab = c("Epcam", "Krt18", "Dcn"), 
                x ="avg_log2FC", 
                y ="p_val_adj", 
                pCutoff = 0.01,
                FCcutoff = 2.0,
                col=c('black', 'black', 'black', 'red')) + theme_pubr() + NoLegend()
dev.off()


pdf("rev epithelial volcano no labs.pdf", width = 5, height = 5)
EnhancedVolcano(epithelium.de.sig, 
                rownames(epithelium.de.sig),
                labSize = 0,
                x ="avg_log2FC", 
                y ="p_val_adj", 
                pCutoff = 0.01,
                FCcutoff = 2.0,
                col=c('black', 'black', 'black', 'red')) + theme_pubr() + NoLegend()
dev.off()

pdf("rev epithelial volcano no labs orange.pdf", width = 5, height = 5)
EnhancedVolcano(epithelium.de.sig, 
                rownames(epithelium.de.sig),
                labSize = 0,
                x ="avg_log2FC", 
                y ="p_val_adj", 
                pCutoff = 0.01,
                FCcutoff = 2.0,
                col=c('gray', 'gray', 'gray', 'orange')) + theme_pubr() + NoLegend()
dev.off()

pdf("rev epithelial volcano no labs blue.pdf", width = 5, height = 5)
EnhancedVolcano(epithelium.de.sig, 
                rownames(epithelium.de.sig),
                labSize = 0,
                x ="avg_log2FC", 
                y ="p_val_adj", 
                pCutoff = 0.01,
                FCcutoff = 2.0,
                col=c('gray', 'gray', 'gray', 'blue')) + theme_pubr() + NoLegend()
dev.off()

# Final
custom.hex.colors <- c("Mesenchymal Cells"="blue","Type I Hair Cells"="orange","Transitional Epithelial Cells"="orange",
                       "Supporting Cells"="orange","Type II Hair Cells"="orange",
                       "Glia"="blue","Roof Cells (CNMD+)"="orange","Pericytes"="blue",
                       "Schwann Cells"="blue","Endothelial Cells"="blue",
                       "Macrophages"="blue","Melanocytes"="blue")

pdf("div-umap-labels-bigger.pdf", width = 3, height = 3)
plot <- DimPlot(so, reduction = "umap")
LabelClusters(plot = plot, id = "ident") + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()

pdf("div-umap-nolabels-bigger.pdf", width = 3, height = 3)
plot <- DimPlot(so, reduction = "umap")
plot + NoLegend() + NoAxes() + scale_fill_manual(values = custom.hex.colors) + scale_colour_manual(values = custom.hex.colors)
dev.off()
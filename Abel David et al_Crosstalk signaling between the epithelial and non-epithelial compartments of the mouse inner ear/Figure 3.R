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

# Load Seurat Object, Visualize Cluster Groups
so <- readRDS("analyzed_utricle_mes_se_p4p6.rds")

# Relevel in descending order
clust.desc.order <- Idents(so) %>% table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% select(-Freq)
colnames(clust.desc.order) <- "clust.desc.order"
clust.desc.order
levels(so) <- clust.desc.order %>% as.matrix()
levels(so) # Confirm reorder

# List of Differentially Expressed Genes, Top 5-10 for each Cluster

select.clusters <- c("Schwann Cells", "Glia", "Pericytes", "Melanocytes", "Endothelial Cells", "Macrophages")

# so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# sig.markers <- so.markers %>% group_by(cluster) %>% filter(p_val_adj < 0.01) %>% top_n(n = 5, wt = avg_log2FC) %>% filter(cluster %in% select.clusters)

seurat.clusters <- Idents(so)

clust <- "Schwann Cells"
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
pdf(paste0(clust, " Violin Plot.pdf"), width = 12, height = 2.5)
#VlnPlot(so, features = rownames(celltype.de.sig), cols = c("#B82E2E", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
select.features <- c("Mpz", "Plp1", "Mbp", "Pmp22")
VlnPlot(so, features = select.features, cols = c("#B82E2E", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
dev.off()

clust <- "Glia"
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
pdf(paste0(clust, " Violin Plot.pdf"), width = 12, height = 2.5)
#VlnPlot(so, features = rownames(celltype.de.sig), cols = c("#0099C6", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
select.features <- c("Ednrb", "Gas7", "Arpc1b", "Fabp7")
VlnPlot(so, features = select.features, cols = c("#0099C6", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
dev.off()

clust <- "Pericytes"
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
pdf(paste0(clust, " Violin Plot.pdf"), width = 12, height = 2.5)
#VlnPlot(so, features = rownames(celltype.de.sig), cols = c("#64C3B7", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
select.features <- c("Itga1", "Gucy1a3", "Rgs5", "Ednra")
VlnPlot(so, features = select.features, cols = c("#64C3B7", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
dev.off()  

clust <- "Melanocytes"
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
pdf(paste0(clust, " Violin Plot.pdf"), width = 12, height = 2.5)
#VlnPlot(so, features = rownames(celltype.de.sig), cols = c("#94594D", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
select.features <- c("Gpnmb", "Tyrp1", "Pmel", "Kcnj13")
VlnPlot(so, features = select.features, cols = c("#94594D", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
dev.off()  

clust <- "Endothelial Cells"
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
pdf(paste0(clust, " Violin Plot.pdf"), width = 12, height = 2.5)
# VlnPlot(so, features = rownames(celltype.de.sig), cols = c("#66AA00", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
select.features <- c("Cldn5", "Flt1", "Ppp1r16b", "Pecam1")
VlnPlot(so, features = select.features, cols = c("#66AA00", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
dev.off()  

clust <- "Macrophages"
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
pdf(paste0(clust, " Violin Plot.pdf"), width = 12, height = 2.5)
#VlnPlot(so, features = rownames(celltype.de.sig), cols = c("#F46FC7", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
select.features <- c("Csf1r", "Ccl2", "Cd83", "C1qb")
VlnPlot(so, features = select.features, cols = c("#F46FC7", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
dev.off()  

clust <- "Mesenchymal Cells"
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
pdf(paste0(clust, " Violin Plot.pdf"), width = 12, height = 2.5)
#VlnPlot(so, features = rownames(celltype.de.sig), cols = c("#3366CC", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())
select.features <- c("Otor", "Dcn", "Coch", "Car3", "Ptgds")
VlnPlot(so, features = select.features, cols = c("#3366CC", "#BBBBBB"), ncol = 5, pt.size = 0) & theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank())

dev.off()  

ggplot2::update_geom_defaults("point",list(size = 5, stroke = 0.5, alpha = 0.7))

# Custom Colors to Gray Out the Non-Epithelial Cells (#BBBBBB)

# custom.hex.colors <- c("Mesenchymal Cells"="#84BF96","Type I Hair Cells"="#BBBBBB","Transitional Epithelial Cells"="#BBBBBB",
#                        "Supporting Cells"="#BBBBBB","Type II Hair Cells"="#BBBBBB",
#                        "Glia"="#FBB268","Roof Cells (CNMD+)"="#BBBBBB","Pericytes"="#3385BB",
#                        "Schwann Cells"="#9D7BBA","Endothelial Cells"="#F57C7C",
#                        "Macrophages"="#7F9D56","Melanocytes"="#FE8D1A")

custom.hex.colors <- c("Mesenchymal Cells"="#3366CC","Type I Hair Cells"="#BBBBBB","Transitional Epithelial Cells"="#BBBBBB",
                       "Supporting Cells"="#BBBBBB","Type II Hair Cells"="#BBBBBB",
                       "Glia"="#0099C6","Roof Cells (CNMD+)"="#BBBBBB","Pericytes"="#64C3B7",
                       "Schwann Cells"="#B82E2E","Endothelial Cells"="#66AA00",
                       "Macrophages"="#F46FC7","Melanocytes"="#94594D")




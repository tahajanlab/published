library(Seurat)
library(pheatmap)
library(magrittr)
library(dplyr)
library(biomaRt)
library(matrixStats)

colMax <- function(data) sapply(data, max, na.rm = TRUE)

disease_maps_2 <- function(object, genes, defined_threshold) {
  print(head(object$celltype))
  genes <- genes[genes %in% rownames(object)]
  RNA_dat <- as.matrix(object@assays$originalexp@counts)
  RNA_dat <- RNA_dat[rownames(RNA_dat) %in% genes,]
  temp_data <- rowSums(as.data.frame(RNA_dat))
  genes <- names(temp_data[temp_data>10])
  int_dat <- as.matrix(object@assays$RNA@data)
  int_dat <- int_dat[rownames(int_dat) %in% genes,]
  object_data <- as.data.frame(int_dat)
  object_mean <- rowMeans(object_data)
  object_sd <- rowSds(as.matrix(object_data))
   
  new_data <- c()
  for(i in 1:nrow(object_data)) {
    object_data[i, ] <- (object_data[i, ] - object_mean[i]) / object_sd[i]
    temp <- c()
    temp$ID <- as.matrix(object$celltype)
    temp$Value <- object_data[i,]
    temp <- aggregate(unlist(temp$Value) ~ temp$ID, temp, mean)
    temp$`unlist(temp$Value)` <- (temp$`unlist(temp$Value)` - mean(temp$`unlist(temp$Value)`))/ sd(temp$`unlist(temp$Value)`)
    # print(temp)
    new_data <- cbind(new_data, temp$`unlist(temp$Value)`)
    # print(temp$`unlist(temp$Value)`)
    temp <- c()
  }
  
  colnames(new_data) <- rownames(object_data) 
  rownames(new_data) <- levels(factor(as.matrix(object$celltype))) # fixed on 10/30/2023
  max <- colMax(as.data.frame(new_data))
  new_data1 <- new_data[, max >= defined_threshold]
  return(new_data1)
}

# Load Seurat Object
data.seurat <- readRDS("data/analyzed_utricle_mes_se_p4p6.rds")

# Set 'celltype' from Idents()
data.seurat$celltype <- Idents(data.seurat)
data.seurat$celltype %>% table()
levels(data.seurat$celltype)
levels(factor(data.seurat$celltype))

# Troubleshoot expression of Myo3, Myo7a, and Espn
png("UMAP.png", width = 10, height = 10, res= 400, units= "in")
DimPlot(data.seurat, group.by = "celltype", label = TRUE)
dev.off()

FeaturePlot(data.seurat, features = c("Myo7a", "Espn"))

# Cell Types:
# Mesenchymal Cells             
# Supporting Cells
# Endothelial Cells
# Type I Hair Cells
# Transitional Epithelial Cells Type II Hair Cells           
# Pericytes
# Roof Cells (CNMD+)
# Glia
# Macrophages
# Schwann Cells
# Melanocytes 

# Load Gene Lists
hhl <- read.csv("data/hereditary_hearing_loss_genes_mouse.csv", header = T)
hhl_genes <- hhl$gene[!is.na(hhl$gene)]
hhl_genes <- hhl_genes[!duplicated(hhl_genes)]
vestibular <- read.csv("data/human_vestibular_defects_genes_mouse.csv", header = T)
vestibular_genes <- vestibular$Gene[!duplicated(vestibular$Gene)]

summary(hhl_genes %in% vestibular_genes) # 14 shared
shared <- vestibular_genes[vestibular_genes %in% hhl_genes] # shared
vestibular_genes <- vestibular_genes[!(vestibular_genes %in% hhl_genes)] # 41
hhl_genes <- hhl_genes[!(hhl_genes %in% shared)] # 143

# Clean Up Hereditary Hearing Loss Gene List and Categories
dns <- hhl$gene[hhl$type == "autosomal_dom_NSHL"] #dominant non syndromic hearing loss
rns <- hhl$gene[hhl$type == "autosomal_rec_NSHL"] #recessive nonsyndromic hearing loss
xlink_ns <- hhl$gene[hhl$type == "x_linked_NSHL"] #x-linked nonsyndromic hearing loss
d_or_rs <- hhl$gene[hhl$type == "autosomal_dom_rec_SHL"] # autosomal dominant/recessive syndromic hearing loss
ds <- hhl$gene[hhl$type == "autosomal_dom_SHL"] # dominant syndromic hearing loss
rs <- hhl$gene[hhl$type == "autosomal_rec_SHL"] # recessive syndromic hearing loss
xlink_rs <-c("COL4A5", "NDP") #x-linked recessive syndromic hearing loss
hhl_list <- c(dns, rns, xlink_ns, ds, rs, xlink_rs)
hhl_labels <- c(rep("DFNA", length(dns)), rep("DFNB", length(rns)),
                rep("XLR-NS", length(xlink_ns)), rep("DSHL", length(ds)), 
                rep("RSHL", length(rs)), rep("XLR-SHL", length(xlink_rs)))

# Clean Up Vestibular Gene List and Categories
dns <- vestibular$Gene[vestibular$Disease %in% c("DFNA11","DFNA3A", "DFNA15", "DNFA9", "DNFA28")] #dominant nonsyndromic hearing loss
rns <- vestibular$Gene[vestibular$Disease %in% c("DFNB84","DFNB102")] #recessive nonsyndromic hearing loss
usher <- vestibular$Gene[vestibular$Disease %in% c("USH1D","USH1J", "USH1F", "USH3", "USH1C", "USH1G")] # usher syndrome
dnfb <- vestibular$Gene[vestibular$Disease %in% c("DNFB4","DNFB37", "DNFB36")] #recessive nonsyndromic hearing loss
vs <- vestibular$Gene[vestibular$Disease %in% c("vestibular_schwannoma")] # vs
motion_sickness <- vestibular$Gene[vestibular$Disease %in% c("motion_sickness")] # motion sickness
superior_canal_dehiscence_syndrome <- vestibular$Gene[vestibular$Disease %in% c("superior_canal_dehiscence_syndrome")] # menieres
menieres <- vestibular$Gene[vestibular$Disease %in% c("familial_menieres","Menieres")] # menieres
vestibular_migrane <- vestibular$Gene[vestibular$Disease %in% c("vestibular_migrane")] # vestibular migrane
vestibular_list <- c(dns, rns, usher, dnfb, vs, motion_sickness, superior_canal_dehiscence_syndrome, menieres, vestibular_migrane)
vestibular_labels <- c(rep("DFNA", length(dns)), rep("DFNB", length(rns)), rep("DNFB", length(dnfb)),
                rep("USH", length(usher)), rep("Vestibular Migrane", length(vestibular_migrane)), 
                rep("Menieres", length(menieres)), rep("Vestibular Schwannoma", length(vs)), 
                rep("Motion Sickness", length(motion_sickness)), rep("SCDS", length(superior_canal_dehiscence_syndrome)))

# Create Hearing Loss Heatmap Data
hhl_heatmap_genes <- as.data.frame(cbind(hhl_list, hhl_labels))
hhl_heatmap_genes_merged <- aggregate(hhl_labels ~ hhl_list, data = hhl_heatmap_genes, paste, collapse = "/")
rownames(hhl_heatmap_genes_merged ) <- hhl_heatmap_genes_merged$hhl_list
hhl_all_genes <- as.character(hhl_heatmap_genes_merged$hhl_list)
hhl_disease_data <- disease_maps_2(data.seurat, hhl_all_genes, 1)

# Create Vestibular Heatmap Data
vest_heatmap_genes <- as.data.frame(cbind(vestibular_list, vestibular_labels))
vest_heatmap_genes_merged <- aggregate(vestibular_labels ~ vestibular_list, data = vest_heatmap_genes, paste, collapse = "/")
rownames(vest_heatmap_genes_merged ) <- vest_heatmap_genes_merged$vestibular_list
vest_all_genes <- as.character(vest_heatmap_genes_merged$vestibular_list)
vest_disease_data <- disease_maps_2(data.seurat, vest_all_genes, 1)

# Produce Heatmaps
hhl_heatmap <- pheatmap(hhl_disease_data, cluster_rows = F, cluster_cols = T,border_color=NA, treeheight_row = 0, treeheight_col = 0)
vest_heatmap <- pheatmap(vest_disease_data, cluster_rows = F, cluster_cols = T,border_color=NA, treeheight_row = 0, treeheight_col = 0)

png("Hearing Loss Heatmap.png", width = 30, height = 3, res= 400, units= "in")
hhl_heatmap
dev.off()

png("Vestibulopathy Heatmap.png", width = 10, height = 3, res= 400, units= "in")
vest_heatmap
dev.off()

pdf("Hearing Loss Heatmap.pdf", width = 30, height = 3)
hhl_heatmap
dev.off()

pdf("Vestibulopathy Heatmap.pdf", width = 10, height = 3)
vest_heatmap
dev.off()


pdf("Hearing Loss Heatmap - No Clustering.pdf", width = 30, height = 3)
hhl_heatmap_noclust <- pheatmap(hhl_disease_data, cluster_rows = F, cluster_cols = F,border_color=NA, treeheight_row = 0, treeheight_col = 0)
print(hhl_heatmap_noclust)
dev.off()

pdf("Vestibulopathy Heatmap - No Clustering.pdf", width = 10, height = 3)
vest_heatmap_noclust <- pheatmap(vest_disease_data, cluster_rows = F, cluster_cols = F,border_color=NA, treeheight_row = 0, treeheight_col = 0)
print(vest_heatmap_noclust)
dev.off()



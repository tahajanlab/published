library(Seurat)
library(Matrix)
library(DoubletFinder)
library(magrittr)
library(dplyr)
library(tidyverse)
library(scater)
library(SingleCellExperiment)
library(swne)
library(Rcpp)
library(M3Drop)
library(future)
library(remotes)
library(renv)
setwd("~/Desktop/SE and MES-P4 P6-Utricle")
source("pipeline.R")
library(writexl)
library(singleCellHaystack)
library(nichenetr)
library(cowplot)
library(gprofiler2)
library(clustifyr)
library(vioplot)

# Load Seurat Object
so <- readRDS("~/Desktop/SE and MES-P4 P6-Utricle/undmg.p4p6.RDS")
so.raw <- readRDS("~/Desktop/SE and MES-P4 P6-Utricle/undmg.p4p6.RDS")
dim(so) #52628  1155

# Enable future parallelization
plan("multiprocess", workers = 8)

# Undamaged cells, treated only with saline (sl)
so <- subset(so, subset = state == "undamaged")

# Quality Control

# Make SCE Object from Seurat Object
sce <- as.SingleCellExperiment(so)

# Rename Plate
table(sce$plate_name) #p4-dt-wt-mes-plt1 173; p4-dt-wt-mes-plt2 0
sce[,sce$plate_code == "180604-30"]$plate_name<- "p4-dt-wt-mes-plt2"
table(sce$plate_name) #p4-dt-wt-mes-plt1 87; p4-dt-wt-mes-plt2 86

# Removing plate "p6-dt-wt-mes-plt2" had greater than 25% (nearly 50%) outliers
sce <- sce[,sce$plate_name != "p6-dt-wt-mes-plt2"]

# Calculate QC Metrics
sce <- addPerCellQC(sce, 
                    subsets = list(
                      Mito = grep("^mt-", rownames(sce)),
                      ERCC = grep("^ERCC-", rownames(sce))
                    )
)

# Filter library size and genes detected
qc.lib <- sce$sum < 1e5
qc.nexprs <- sce$detected < 1000

# Spike In and Mito
qc.spike <- sce$subsets_ERCC_percent  > 10 
summary(qc.spike) # None removed
qc.mito <- sce$subsets_Mito_percent > 10
summary(qc.mito) # 35 Removed

# Adaptive Thresholds
qc.lib2 <- isOutlier(sce$sum, nmads = 3, log = TRUE, type = "both")
summary(qc.lib2) # Outliers = 109
qc.nexprs2 <- isOutlier(sce$detected, nmads = 3, log = TRUE, type = "both")
summary(qc.nexprs2) # Outliers = 117

# Adaptive threshold cutoffs
attr(qc.lib2, "thresholds") # lower 106440.9, higher 2671637.6 
attr(qc.nexprs2, "thresholds") # lower 907.1321, higher 8679.7008 

# Aggregate Poor Quality Cells

discard2 <- qc.lib2 | qc.nexprs2 | qc.spike | qc.mito
summary(discard2) # Outliers (aggregate) = 148
sce$outliers <- discard2

# QC Plots

# QC violin plots highlighting outliers
pal <- c("black", "#ff2a00")
ggplot2::update_geom_defaults("point",list(size = 2, stroke = 0.5, alpha = 1))

# Increasing number of total genes with increasing total count
plotColData(sce, x = "sum", y="detected", colour_by="cell_type") 

# Library Size
plot.y <- "sum"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce, y=plot.y, colour_by=color.plot) + ggtitle(title) + scale_color_manual(values= pal) +  scale_y_log10() 


# QC Violin highlighting metadata

# Total Counts ('sum'): sequence pool, age, genotype, cell_type, plate_name
plot.x <- "seq_pool"
plot.y <- "sum"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + ggtitle("Total Count")+ scale_color_manual(values= pal)

plot.x <- "Age"
plot.y <- "sum"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + ggtitle("Total Count")+ scale_color_manual(values= pal)

plot.x <- "Genotype"
plot.y <- "sum"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + ggtitle("Total Count")+ scale_color_manual(values= pal)

plot.x <- "cell_type"
plot.y <- "sum"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + ggtitle("Total Count")+ scale_color_manual(values= pal)

plot.x <- "plate_name"
plot.y <- "sum"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + ggtitle("Total Count")+ scale_color_manual(values= pal) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Plate position plot
# plotting function requires the cell position data (alpha-numeric combination) to be stored in an field called 'plate_position'
sce$plate_position <- sce$Cell_position
plotPlatePosition(sce, size_by = "sum", colour_by = NULL, point_alpha = 1)

# Outliers Sum versus nExpres
plotColData(sce, x = "sum", y="detected", colour_by="outliers") + scale_color_manual(values= pal)

plotColData(sce, x = "sum", y="subsets_Mito_percent", colour_by = "outliers", other_fields="Cell_type") + scale_color_manual(values= pal) + facet_wrap(~Cell_type)

plotColData(sce, x = "sum", y = "subsets_Mito_percent", colour_by = "outliers", other_fields="seq_pool") + scale_color_manual(values= pal)  + facet_wrap(~seq_pool)

# Remove Outlier Cells
sce <- sce[,!sce$outliers]
dim(sce) #52628  1007

# Dump MT and ERCC Genes
# Remove MTs genes
is.control <- grepl("^mt-", rownames(sce)) # find all rownames beginning with "ERCC-" or "mt-"
length(which(is.control)) # check that there are 37 control genes 
sce <- sce[!is.control,]  # remove 37 MTs
dim(sce) #52591  1007

# Remove ERCCs genes
is.control <- grepl("^ERCC-", rownames(sce)) # find all rownames beginning with "ERCC-" or "mt-"
length(which(is.control)) # check that there are 92 control genes 
sce <- sce[!is.control,]  # remove 92 ERCCs 
dim(sce) #52499  1007

# Removing genes not expressed in at least 3 cells
keep_feature <- rowSums(counts(sce)) > 3 
summary(keep_feature) # 22996 genes to remove
sce <- sce[keep_feature, ] # new sce w/ only relevant genes
dim(sce) #25815  1007

# Convert SCE Object to Seurat Object
so <- as.Seurat(sce)
dim(so) #25815  1007

# Finding Variable Features

# Normalizing the Data
so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)

# Isolate Out mesenchymal cells only -- Later pipeline analysis
# so <- subset(so, subset = Cell_type == "mes")

# FINDING VARIABLE FEATURES
selection.method <- "vst" #vst, mvp, disp, or m3drop

if(selection.method != "m3drop") {
  so <- FindVariableFeatures(so, selection.method = selection.method, nfeatures = 2000)
} else {
  # [X] M3 Drop -- For Finding Highly Variable Features, Michaelis Menten Modelling
  # Find highly variable features, requires normalized data
  # Convert Seurat object to SCE
  sce2 <- as.SingleCellExperiment(so)
  # exp_mat <- M3DropConvertData(sce2@assays[["normcounts"]], is.log = FALSE, is.counts = FALSE)
  exp_mat <- M3DropConvertData(sce2, is.log = TRUE, is.counts = FALSE)
  # [1] "Removing  0 undetected genes."
  dim(exp_mat) #25695   955
  fits <- M3DropDropoutModels(exp_mat)
  
  # Sum absolute residuals
  data.frame(MM=fits$MMFit$SAr, Logistic=fits$LogiFit$SAr, DoubleExpo=fits$ExpoFit$SAr)
  # Sum squared residuals
  data.frame(MM=fits$MMFit$SSr, Logistic=fits$LogiFit$SSr, DoubleExpo=fits$ExpoFit$SSr)
  DE_genes <- M3DropFeatureSelection(exp_mat, mt_method="fdr", mt_threshold = 0.01)
  dim(DE_genes) # 2945    4
  heat_out <- M3DropExpressionHeatmap(DE_genes$Gene, exp_mat)
  cell_populations <- M3DropGetHeatmapClusters(heat_out, k=4, type = "cell")
  marker_genes <- M3DropGetMarkers(exp_mat, cell_populations)
  dim(marker_genes) #25695     3
  head(marker_genes)
  # Convert SCE2 to SO
  so <- as.Seurat(sce2)
  # END M3Drop
  
  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000)
  VariableFeatures(so) <- rownames(DE_genes) # If wanting to use M3Drop DE_genes
}

top10 <- head(VariableFeatures(so), 10)
plot1 <- VariableFeaturePlot(so) + NoLegend()
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE) + NoLegend()
plot1 + plot2

# Regressing out batch effects using the ScaleData Function, adding an argument vars.to.regress = c(sequence_pool, Age, Sort.Date, Plate Position, AGe, Genotype, Well Volume, Enzyme, Buffer Concentration, Treatment)
# Scaling the Data

all.genes <- rownames(so)

# Cell Cycle Regression

# Cell Cycle Scoring

# Cell cycle markers from Tirosh et al 2015 is loaded with Seurat
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Perform Cell Cycle Scoring
so <- CellCycleScoring(so, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Visualize the distribution of cell cycle markers
RidgePlot(so, features = c("Pcna","Top2a", "Mcm6", "Mki67"), ncol = 2)

# Running PCA Prior to Cell Cycle Regression
so <- ScaleData(so, features = all.genes)
so <- RunPCA(so, features = VariableFeatures(so))
DimPlot(so, reduction = "pca")

# Regress out cell cycle scores during data scaling

vars.to.regress <- c("seq_pool", "plate_position", "Genotype", "Treatment","sort_date", "sac_time", "S.Score", "G2M.Score")
# vars.to.regress <- c("seq_pool", "plate_position", "Genotype", "Treatment","sort_date", "sac_time", "CC.Difference")

so <- ScaleData(so, features = all.genes, vars.to.regress = vars.to.regress)
so <- RunPCA(so, features = VariableFeatures(so))
print(so[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(so, dims = 1:2, reduction = "pca")
DimPlot(so, reduction = "pca")
DimPlot(so, reduction = "pca", group.by = "Cell_type")
DimPlot(so, reduction = "pca", split.by = "Cell_type", group.by = "Cell_type")

DimHeatmap(so, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(so, dims = 1:15, cells = 500, balanced = TRUE)

# Determine Dimensionality of the Dataset
so <- JackStraw(so, num.replicate = 100)
so <- ScoreJackStraw(so, dims = 1:20)
JackStrawPlot(so, dims = 1:20)

ElbowPlot(so)

# ChooseR

# To determine the resolution, resolution  0.8, selection.method = vst

# # Begin ChooseR Workflow
npcs <- 20
resolutions <- c(0.2,0.4,0.6,0.8, 1, 1.6,2,4,6,8,12,16)
assay <- "origexp"
reduction <- "pca"
results_path <- paste0("chooseR-results/", selection.method, "/")

obj <- so

# Run pipeline
for (res in resolutions) {
   message(paste0("Clustering ", res, "..."))
   message("\tFinding ground truth...")

   # "Truths" will be stored at glue::glue("{reduction}.{assay}_res.{res}")
   obj <- find_clusters(
      obj,
      reduction = reduction,
      assay = assay,
      npcs = npcs,   ###change made from original github repository
      resolution = res
   )
   clusters <- obj[[glue::glue("{reduction}.{assay}_res.{res}")]]

   # Now perform iterative, sub-sampled clusters
   results <- multiple_cluster(
      obj,
      n = 100,
      size = 0.8,
      npcs = npcs,
      res = res,
      reduction = reduction,
      assay = assay
   )

   # Now calculate the co-clustering frequencies
   message(paste0("Tallying ", res, "..."))
   # This is the more time efficient vectorisation
   # However, it exhausts vector memory for (nearly) all datasets
   # matches <- purrr::map(columns, find_matches, df = results)
   # matches <- purrr::reduce(matches, `+`)
   columns <- colnames(dplyr::select(results, -cell))
   mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
   i <- 1 # Counter
   for (col in columns) {
      message(paste0("\tRound ", i, "..."))
      mtchs <- Reduce("+", list(
         mtchs,
         find_matches(col, df = results)
      ))
      i <- i + 1
   }

   message(paste0("Scoring ", res, "..."))
   mtchs <- dplyr::mutate_all(
      dplyr::as_tibble(mtchs),
      function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
   )

   # Now calculate silhouette scores
   message(paste0("Silhouette ", res, "..."))
   sil <- cluster::silhouette(
      x = as.numeric(as.character(unlist(clusters))),
      dmatrix = (1 - as.matrix(mtchs))
   )
   saveRDS(sil, paste0(results_path, "silhouette_", res, ".rds"))

   # Finally, calculate grouped metrics
   message(paste0("Grouping ", res, "..."))
   grp <- group_scores(mtchs, unlist(clusters))
   saveRDS(grp, paste0(results_path, "frequency_grouped_", res, ".rds"))
   sil <- group_sil(sil, res)
   saveRDS(sil, paste0(results_path, "silhouette_grouped_", res, ".rds"))
}

saveRDS(obj, paste0(results_path, "clustered_data.rds"))

# Create silhouette plot
# Read in scores and calculate CIs
scores <- purrr::map(
   paste0(results_path, "silhouette_grouped_", resolutions, ".rds"),
   readRDS
)
scores <- dplyr::bind_rows(scores) %>%
   dplyr::group_by(res) %>%
   dplyr::mutate("n_clusters" = dplyr::n()) %>%
   dplyr::ungroup()
meds <- scores %>%
   dplyr::group_by(res) %>%
   dplyr::summarise(
      "boot" = list(boot_median(avg_sil)),
      "n_clusters" = mean(n_clusters)
   ) %>%
   tidyr::unnest_wider(boot)

writexl::write_xlsx(meds, paste0(results_path, "median_ci.xlsx"))

# Find thresholds
threshold <- max(meds$low_med)
choice <- as.character(
   meds %>%
      dplyr::filter(med >= threshold) %>%
      dplyr::arrange(n_clusters) %>%
      tail(n = 1) %>%
      dplyr::pull(res)
)

# Plot Silhouette Scores and Optimal Resolution
ggplot(meds, aes(factor(res), med)) +
   geom_crossbar(
      aes(ymin = low_med, ymax = high_med),
      fill = "grey",
      size = 0.25
   ) +
   geom_hline(aes(yintercept = threshold), colour = "blue") +
   geom_vline(aes(xintercept = choice), colour = "red") +
   geom_jitter(
      data = scores,
      aes(factor(res), avg_sil),
      size = 0.35,
      width = 0.15
   ) +
   scale_x_discrete("Resolution") +
   scale_y_continuous(
      "Silhouette Score",
      expand = c(0, 0),
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.25),
      oob = scales::squish
   ) +
   cowplot::theme_minimal_hgrid() +
   theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
   )

ggsave(
   filename = paste0(results_path, "silhouette_distribution_plot.png"),
   dpi = 300,
   height = 3.5,
   width = 3.5,
   units = "in"
) # Resolution 2.0

# Finally, a dot plot of silhouette scores to help identify less robust clusters
# The initial pipe is to order the clusters by silhouette score
scores %>%
   dplyr::filter(res == choice) %>%
   dplyr::arrange(dplyr::desc(avg_sil)) %>%
   dplyr::mutate_at("cluster", ordered, levels = .$cluster) %>%
   ggplot(aes(factor(cluster), avg_sil)) +
   geom_point() +
   scale_x_discrete("Cluster") +
   scale_y_continuous(
      "Silhouette Score",
      expand = c(0, 0),
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.25),
      oob = scales::squish
   ) +
   cowplot::theme_minimal_grid() +
   theme(
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      axis.line.x = element_line(colour = "black"),
      axis.line.y = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black"),
   )

ggsave(
   filename = paste0(results_path, "silhouette_point_plot_", choice, ".png"),
   dpi = 300,
   height = 3.5,
   width = 3.5,
   units = "in"
)

reduction <- "pca"
assay <- "origexp"
choice <- 2
results_path <- "chooseR-results/"

# Load in the object containing the clustered results
obj <- readRDS(paste0(results_path, "clustered_data.rds"))

# First is a cluster average co-clustering heatmap
# Read the data
grp <- readRDS(paste0(results_path, "frequency_grouped_", choice, ".rds"))

# As the data is symmetrical, we do not need the upper triangle
grp <- grp %>%
   pivot_wider(names_from = "cell_2", values_from = "avg_percent") %>%
   select(str_sort(colnames(.), numeric = T)) %>%
   column_to_rownames("cell_1")
grp[lower.tri(grp)] <- NA
grp <- grp %>%
   as_tibble(rownames = "cell_1") %>%
   pivot_longer(-cell_1, names_to = "cell_2", values_to = "avg_percent") %>%
   mutate_at("cell_2", ordered, levels = unique(.$cell_1)) %>%
   mutate_at("cell_1", ordered, levels = unique(.$cell_1))

# And plot!
plot <- ggplot(grp, aes(factor(cell_1), cell_2, fill = avg_percent)) +
   geom_tile() +
   scale_x_discrete("Cluster", expand = c(0, 0)) +
   scale_y_discrete(
      "Cluster",
      limits = rev(levels(grp$cell_2)),
      expand = c(0, 0)
   ) +
   scale_fill_distiller(
      " ",
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      palette = "RdYlBu",
      na.value = "white"
   ) +
   coord_fixed() +
   theme(
      axis.ticks = element_line(colour = "black"),
      axis.text = element_text(size = 6),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.position = c(0.9, 0.9)
   ) +
   guides(fill = guide_colorbar(barheight = 3, barwidth = 1))

ggsave(
   plot = plot,
   filename = paste0(results_path, "coclustering_heatmap_", choice, ".png"),
   dpi = 300,
   height = 3.5,
   width = 3.5,
   units = "in"
)

# Let's add the silhouette scores to the Seurat object
choice <- 0.8
sil_scores <- readRDS(paste0(results_path, "silhouette_", choice, ".rds"))
sil_scores <- as.data.frame(sil_scores[, 3], row.names = Seurat::Cells(so))
colnames(sil_scores) <- c("sil_score")
so <- AddMetaData(so, metadata = sil_scores)

# UMAP Dimension Reduction

# Cluster the cells
so <- FindNeighbors(so, dims = 1:20)
so <- FindClusters(so, resolution = 0.8) # the resolution = 2

head(Idents(so),5)

# Run non-linear dimensional reduction (UMAP/tSNE)
so <- RunUMAP(so, dims = 1:20)

DimPlot(so, reduction = "umap")
DimPlot(so, reduction = "umap", split.by = "Cell_type")
DimPlot(so, reduction = "umap", group.by = "seq_pool")
DimPlot(so, reduction = "umap", group.by = "Age")
DimPlot(so, reduction = "umap", group.by = "Genotype")
DimPlot(so, reduction = "umap", group.by = "plate_name")
DimPlot(so, reduction = "umap", group.by = "plate_code")
DimPlot(so, reduction = "umap", group.by = "plate_position")

# ChooseR UMAPS

gg_color <- function(n) {
   hues <- seq(15, 375, length = n + 1)
   colours <- hcl(h = hues, c = 100, l = 65)[1:n]
   return(colours)
}

plot <- DimPlot(
   so,
   reduction = "umap",
   group.by = glue::glue("{reduction}.{assay}_res.{choice}"),
   pt.size = 0.5,
   # cols = gg_color(6) # Only necessary if you have ordered your clusters
)

ggsave(
   plot = plot,
   filename = paste0(results_path, choice, "_cluster_umap.png"),
   dpi = 300,
   height = 5,
   width = 5,
   units = "in"
)

# We also find it useful to visualise the silhouette scores on the UMAP!
FeaturePlot(
  so,
  "sil_score",
  reduction = "umap",
  pt.size = 5,
  min.cutoff = -1,
  max.cutoff = 1
) +
  scale_colour_distiller(
    palette = "RdYlBu",
    labels = c(-1, 0, 1),
    breaks = c(-1, 0, 1),
    limits = c(-1, 1)
  )

# Feature Plots with Cell Population Markers for Cell Type Annotation

# Feature Plot

features.use <-  c("Ocm", "Myo7a", "Atoh1", "Spp1", "Anxa4", "Sparcl1", "Dclk1", "Sox2", "Islr", "Pbx18", "PectB", "Tectb", "Pou4f3", "Tbx18")
FeaturePlot(object = so, features = features.use) & NoAxes()

features.use <- c("Mlana", "Cndp2", "Emcn", "Bsnd", "Kcnq1", "Cdh19", "Cav1", "Igf1", "Serpine2", "Ednrb", "Dct", "Tryp1", "Pmel")
FeaturePlot(object = so, features = features.use) & NoAxes()

# Erythrocyte Markers
features.use <- c("Hba1", "Hba2", "Hbb")
FeaturePlot(object = so, features = features.use) & NoAxes()

# Lymphocytic Markers
features.use <- c("Cldn5", "Lyve1", "Prox1")
FeaturePlot(object = so, features = features.use) & NoAxes()

# Endothelial Markers
features.use <- c("Sele", "Cldn5", "Vwf", "Cdh5")
FeaturePlot(object = so, features = features.use) & NoAxes()

# T Cell Markers
features.use <- c("Cd3d", "Cd3g", "Cd3e", "Lck")
FeaturePlot(object = so, features = features.use) & NoAxes()

# Macrophages 
features.use <- c("Lyz", "Aif1", "Hla-dra", "Cd68", "Itgax")
FeaturePlot(object = so, features = features.use) & NoAxes()

# Fibroblast Mesenchymal
features.use <- c("Lum", "Dcn", "Vim", "Pdgfra", "Col1a2")
FeaturePlot(object = so, features = features.use) & NoAxes()

# SWNE Dimension Reduction

# Start SWNE
so@assays[['RNA']] <- so@assays$originalexp
norm.counts <- ExtractNormCounts(so, obj.type = "seurat", rescale.method = "log")
dim(norm.counts) # 25695   955
var.genes <- VariableFeatures(so)
cell.clusters <- Idents(so)
k.range <- seq(2,22,2)
k.err <- FindNumFactors(norm.counts[var.genes,], k.range = k.range, n.cores = 8, do.plot = T)
genes.embed <- c("Ocm", "Myo7a", "Atoh1", "Spp1", "Anxa4", "Sparcl1", "Dclk1")
k <- 20
nmf.res <- RunNMF(norm.counts[var.genes,], k = k) # Stored NMF Results 

snn <- as(so@graphs$originalexp_snn, "dgCMatrix")
knn <- as(so@graphs$originalexp_nn, "dgCMatrix") ## Extract kNN matrix

alpha.exp <- 1.5
snn.exp <- 0.75
n_pull <- 3

swne.embedding <- EmbedSWNE(nmf.res$H, SNN = snn, alpha.exp = alpha.exp, snn.exp = snn.exp, n_pull = n_pull)

swne.embedding$H.coords$name <- ""

nmf.res$W <- ProjectFeatures(norm.counts, nmf.res$H, n.cores = 8)

swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, genes.embed, n_pull = n_pull)

color.seed <- 42
PlotSWNE(swne.embedding, alpha.plot = 0.7, sample.groups = cell.clusters, do.label = F, label.size = 5, pt.size = 3, show.legend = T, seed = color.seed)
# Feature Plots with SWNE

gene.use <- features.use[1]
gene.expr <- norm.counts[gene.use,]
FeaturePlotSWNE(swne.embedding, gene.expr, gene.use, alpha.plot = 0.7, label.size = 5, pt.size = 3)

for(gene.use in features.use)
{
  print(gene.use)
}

# Seurat Clustered Heat Map and DEGs

# Differential Gene Expression - Find All Markers, report only positive ones
so.markers <- FindAllMarkers(so, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top2.byclust <- so.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)
# Idents(so)
FeaturePlot(so, top2.byclust$gene, order = TRUE) & NoAxes()

so.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

# Heatmap (adj p-value < 0.01)
sig.markers <- so.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01)
DoHeatmap(so, features = sig.markers$gene) + NoLegend() + NoAxes()

# Top 10
sig.markers <- so.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(so, features = sig.markers$gene) + NoLegend()

# gProfileR

cluster0.DEG <- so.markers %>% group_by(cluster) %>% filter(cluster == 0) %>% .$gene
gostres.cluster0 <- gost(query = cluster0.DEG, organism = "mmusculus", ordered_query = TRUE)
gostplot(gostres.cluster0, capped = FALSE, interactive = TRUE)

write_lines(cluster0.DEG, "cluster0-DEG.txt")

cluster1.DEG <- so.markers %>% group_by(cluster) %>% filter(cluster == 1) %>% .$gene
gostres.cluster1 <- gost(query = cluster1.DEG, organism = "mmusculus", ordered_query = TRUE)
gostplot(gostres.cluster1, capped = FALSE, interactive = TRUE)

write_lines(cluster1.DEG, "cluster1-DEG.txt")

cluster2.DEG <- so.markers %>% group_by(cluster) %>% filter(cluster == 2) %>% .$gene
gostres.cluster2 <- gost(query = cluster2.DEG, organism = "mmusculus", ordered_query = TRUE)
gostplot(gostres.cluster2, capped = FALSE, interactive = TRUE)

write_lines(cluster2.DEG, "cluster2-DEG.txt")

cluster3.DEG <- so.markers %>% group_by(cluster) %>% filter(cluster == 3) %>% .$gene
gostres.cluster3 <- gost(query = cluster3.DEG, organism = "mmusculus", ordered_query = TRUE)
gostplot(gostres.cluster3, capped = FALSE, interactive = TRUE)

write_lines(cluster3.DEG, "cluster3-DEG.txt")

cluster4.DEG <- so.markers %>% group_by(cluster) %>% filter(cluster == 4) %>% .$gene
gostres.cluster4 <- gost(query = cluster4.DEG, organism = "mmusculus", ordered_query = TRUE)
gostplot(gostres.cluster4, capped = FALSE, interactive = TRUE)

write_lines(cluster4.DEG, "cluster4-DEG.txt")


cluster5.DEG <- so.markers %>% group_by(cluster) %>% filter(cluster == 5) %>% .$gene
gostres.cluster5 <- gost(query = cluster5.DEG, organism = "mmusculus", ordered_query = TRUE)
gostplot(gostres.cluster5, capped = FALSE, interactive = TRUE)

write_lines(cluster5.DEG, "cluster5-DEG.txt")

cluster6.DEG <- so.markers %>% group_by(cluster) %>% filter(cluster == 6) %>% .$gene
gostres.cluster6 <- gost(query = cluster6.DEG, organism = "mmusculus", ordered_query = TRUE)
gostplot(gostres.cluster6, capped = FALSE, interactive = TRUE)

write_lines(cluster6.DEG, "cluster6-DEG.txt")


cluster7.DEG <- so.markers %>% group_by(cluster) %>% filter(cluster == 7) %>% .$gene
gostres.cluster7 <- gost(query = cluster7.DEG, organism = "mmusculus", ordered_query = TRUE)
gostplot(gostres.cluster7, capped = FALSE, interactive = TRUE)

write_lines(cluster7.DEG, "cluster7-DEG.txt")

cluster8.DEG <- so.markers %>% group_by(cluster) %>% filter(cluster == 8) %>% .$gene
gostres.cluster8 <- gost(query = cluster8.DEG, organism = "mmusculus", ordered_query = TRUE)
gostplot(gostres.cluster8, capped = FALSE, interactive = TRUE)
write_lines(cluster8.DEG, "cluster8-DEG.txt")

# Renamed UMAP

so <- RenameIdents(so, "0" = "Mesenchymal Cells", "1" = "Type I Hair Cells", "2" = "Transitional Epithelial Cells", "3" = "Supporting Cells", "4" = "Type II Hair Cells", "5" = "Glia", "6" = "Roof Cells (CNMD+)", "7" = "Pericytes", "8" = "Endothelial Cells")
table(Idents(so))
plot <- DimPlot(so, reduction = "umap")
LabelClusters(plot = plot, id = "ident")

macrophages.cells <- c("A5_Taha_p4_dt_wt_mes_plt1_180604_29",
                       "B3_Taha_p4_dt_wt_mes_plt1_180604_30",
                       "D10_Taha_p6_dt_wt_mes_plt3_180604_24R",
                       "F6_Taha_p6_dt_wt_mes_plt3_180604_24R",
                       "G11_Taha_p6_dt_wt_mes_plt3_180604_24R",
                       "H1_Taha_p4_dt_wt_mes_plt1_180604_30",
                       "A6_Taha_p6_sl_dtr_mes_plt1",
                       "A9_Taha_p6_sl_dtr_mes_plt2",
                       "B1_Taha_p6_sl_dtr_mes_plt2",
                       "E5_Taha_p6_sl_dtr_mes_plt1")

melanocytes.cells <- c("D12_Taha_p4_dt_wt_mes_plt1_180604_30",
                       "D8_Taha_p6_dt_wt_mes_plt1_180604_22",
                       "G9_Taha_p4_dt_wt_se_plt2_180604_28")

schwann.cells <- c("B3_Taha_p6_dt_wt_mes_plt3_180604_24R",
                   "B9_Taha_p6_dt_wt_mes_plt1_180604_22",
                   "C1_Taha_p6_dt_wt_mes_plt3_180604_24R",
                   "C11_Taha_p6_dt_wt_mes_plt1_180604_22",
                   "D4_Taha_p6_dt_wt_mes_plt3_180604_24R",
                   "E1_Taha_p6_dt_wt_mes_plt1_180604_22",
                   "E12_Taha_p6_dt_wt_mes_plt1_180604_22",
                   "E6_Taha_p6_dt_wt_mes_plt1_180604_22",
                   "F4_Taha_p4_dt_wt_mes_plt1_180604_29",
                   "F6_Taha_p4_dt_wt_mes_plt1_180604_29",
                   "G8_Taha_p4_dt_wt_mes_plt1_180604_29",
                   "H5_Taha_p6_dt_wt_mes_plt1_180604_22",
                   "A10_Taha_p6_sl_dtr_mes_plt2",
                   "C5_Taha_p6_sl_dtr_mes_plt2",
                   "D4_Taha_p6_sl_dtr_mes_plt1",
                   "E10_Taha_p6_sl_dtr_mes_plt1",
                   "F1_Taha_p6_sl_dtr_mes_plt2",
                   "F10_Taha_p6_sl_dtr_mes_plt1",
                   "F5_Taha_p6_sl_dtr_mes_plt2",
                   "G4_Taha_p6_sl_dtr_mes_plt1",
                   "G9_Taha_p6_sl_dtr_mes_plt1")
so2 <- so
Idents(so2, cells = macrophages.cells) <- "Macrophages"
Idents(so2, cells = melanocytes.cells) <- "Melanocytes"
Idents(so2, cells = schwann.cells) <- "Schwann Cells"

plot <- DimPlot(so2, reduction = "umap")
LabelClusters(plot = plot, id = "ident")
LabelClusters(plot = plot, id = "ident") + NoLegend()

# Find DEG with new Clusters

so2.markers <- FindAllMarkers(so2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Heatmap (adj p-value < 0.01)
sig.markers <- so2.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01)
DoHeatmap(so2, features = sig.markers$gene) + NoLegend() + NoAxes()

# Top 10
sig.markers <- so2.markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.01) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(so2, features = sig.markers$gene) + NoLegend()

# SingleCell Haystack

so.hs <- haystack(so, assay = "originalexp", slot = "data", coord = "umap", cutoff = 1, method = "highD") # Cutoff FDR 1%?
Idents(so) %>% table()

str(so.hs) #results: 25695 all genes
gene <- "Mia"
plot_gene_haystack(so, dim1 = 1, dim2 = 2, assay = "originalexp", slot = "data", coord = "umap", gene = gene)

# res.top <- show_result_haystack(res.haystack = so.hs, n = 1000)
res.top <- show_result_haystack(res.haystack = so.hs)
dim(res.top) #25695     4
res.top$bonf <- p.adjust(10^res.top$log.p.vals, method = "bonferroni")
res.top <- res.top %>% filter(bonf < 0.01)
dim(res.top) #2221    5
# res.top
genes.top <- row.names(res.top)

FeaturePlot(so, rownames(res.top)[1:4], order = TRUE) & NoAxes()
FeaturePlot(so, rownames(res.top)[1:12], order = TRUE) & NoAxes()

plot_gene_set_haystack(so, dim1 = 1, dim2 = 2, assay = "originalexp", slot = "data", coord = "umap", gene = genes.top)

features <- rownames(res.top)[1:12]

RidgePlot(so, features = features, ncol = 3)
VlnPlot(so, features = features)
DotPlot(so, features = features) + RotatedAxis()

# Get UMAP Coord
# Detection
dat.expression <- GetAssayData(so, slot = "counts")
dat.detection <- dat.expression > 1
dim(dat.expression) #25695   955
umap.coord <- Embeddings(so, reduction = "umap")
dim(umap.coord) # 955   2
res.hc <- hclust_haystack(umap.coord, detection = dat.detection, genes = genes.top)
res.hc.clusters <- cutree(res.hc, k = 6)
table(res.hc.clusters)
length(res.hc.clusters) #5867

pl <- lapply(1:6, function(cluster) {
  gene.set <- names(res.hc.clusters)[res.hc.clusters == cluster]
  plot.title <- paste0("Cluster ", cluster)
  p <- plot_gene_set_haystack(so, dim1 = 1, dim2 = 2, assay = "originalexp", slot = "data", coord = "umap", gene = gene.set)
  p + ggtitle(plot.title) + theme(legend.title = element_text(size = 8))
})
plot_grid(plotlist = pl, ncol = 2)

Idents(so) <- 'Cell_type'
levels(so)

# Find differentially expressed features between SE and MES Cell_type
so.celltype.de.markers <- FindMarkers(so, ident.1 = "mes", ident.2 = "se")

head(so.celltype.de.markers)

# Subclustering - Mesenchymal Cells

so.nonep <- subset(so2, idents = c("Mesenchymal Cells", "Endothelial Cells", "Pericytes", "Glia", "Schwann Cells"))
so.nonep <- FindNeighbors(nonep, dims = 1:20)
so.nonep <- FindClusters(nonep, resolution = 1)
so.nonep <- RunUMAP(nonep, dims = 1:20)
DimPlot(nonep, reduction = "umap")
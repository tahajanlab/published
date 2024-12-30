
## ----------------------------------------------------------------------------------------------------------------
library(Seurat)
library(DoubletFinder)

## ----------------------------------------------------------------------------------------------------------------
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


# P7cochlear
# read in the matrix files that were gzipped using Read10X from Rose et al 2023
# data link: https://umgear.org/expression.html?gene_symbol=Pou3f4&gene_symbol_exact_match=1&is_multigene=0&layout_id=8c69d179
# data link (alternate): UMgEAR.org/mesenchyme
# data can also be accessed through NCBI GEO (GSE217727, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE217727)
# P7, Mouse, scRNA-seq, whole cochlea (Hertzano/Kelley)
## ----------------------------------------------------------------------------------------------------------------
raw_data <- Read10X(data.dir = '/Users/Desktop/wt-mesenchyme-project/P7cochlear_rdata/')

# read in  metadata
## ----------------------------------------------------------------------------------------------------------------
metadata <- read.csv('/Users/Desktop/wt-mesenchyme-project/metadata.csv')

# create your seurat object
## ----------------------------------------------------------------------------------------------------------------
Cochlear <- CreateSeuratObject(counts = raw_data, meta.data = metadata)
Cochlear
Cochlear[[]]

# store mitochondrial percentage in object meta data and subset mes from P7 cochlea
## ----------------------------------------------------------------------------------------------------------------
Cochlear <- PercentageFeatureSet(Cochlear, pattern = "^mt-", col.name = "percent.mt")
Cochlear <- PercentageFeatureSet(Cochlear, pattern = "^Rps-", col.name = "percent.rp")
Cochlear <- PercentageFeatureSet(Cochlear, pattern = "^Hb[^(p)]", col.name = "percent.hb")

# filter out cells that have more than the specified thresholds following their filtered Mes data 
## ----------------------------------------------------------------------------------------------------------------
Cochlear <- subset(Cochlear, subset = percent.mt < 7 & percent.rp < 30 & percent.hb < 2.5)
Cochlear

# percent.mt check
## ----------------------------------------------------------------------------------------------------------------
FeaturePlot(
  object = Cochlear,
  features = "percent.mt",
  cols = c("blue", "red"),  
  reduction = "umap"        
) +
  ggtitle("UMAP: Mitochondrial Percentage") +
  theme_minimal()

# Normalize data
Cochlear <- NormalizeData(Cochlear, normalization.method = "LogNormalize", scale.factor = 10000)

# Find variable genes
Cochlear <- FindVariableFeatures(Cochlear, selection.method = "vst", nfeatures = 2000)

# scale data
all.genes <- rownames(Cochlear)
Cochlear <- ScaleData(Cochlear, features = all.genes)

# run PCA
Cochlear <- RunPCA(Cochlear, features = VariableFeatures(object = Cochlear))

# Find nearest neighbors and clusters
Cochlear <- FindNeighbors(Cochlear, dims = 1:10)
Cochlear <- FindClusters(Cochlear, resolution = 0.5)

# run UMAP and dimension reduction
Cochlear <- RunUMAP(Cochlear, dims = 1:10)
DimPlot(Cochlear, reduction = "umap")

# check mesenchymal marker genes from Rose et al 2023 (from their figures/plots)
FeaturePlot(Cochlear, features = c('Pou3f4', 'Scara3', 'Car3', 'Tgfbi', 'Kcnk2', 'Dlx5', 'Col15a1', 'Dcn', 'Otos', 'Emilin2', 'Ibsp'))

# reduce clustering resolution to match Rose et al resolution/number of clusters
Cochlear <- FindClusters(Cochlear, resolution = 0.3)

# Show UMAP & violin plot to identify mesenchymal cells 
DimPlot(Cochlear, reduction = "umap")
VlnPlot(Cochlear, features = c('Pou3f4', 'Scara3', 'Car3', 'Tgfbi', 'Kcnk2', 'Dlx5', 'Col15a1', 'Dcn', 'Otos', 'Emilin2', 'Ibsp'))

#  We selected the following clusters as mesenchyme:0,1,2,3,4 
#  subset mes
CochlearMes <- subset(Cochlear, subset = RNA_snn_res.0.3 %in% c("0", "1", "2", "3", "4"))
CochlearMes

#  Show mesenchyme UMAP
DimPlot(CochlearMes, reduction = "umap")

# redo analysis with subseted data: normalization, feature selection, PCA, UMAP 
CochlearMes <- NormalizeData(CochlearMes, normalization.method = "LogNormalize", scale.factor = 10000)
CochlearMes <- FindVariableFeatures(CochlearMes, selection.method = "vst", nfeatures = 2000)
CochlearMesall.genes <- rownames(CochlearMes)
CochlearMes <- ScaleData(CochlearMes, features = CochlearMesall.genes)
CochlearMes <- RunPCA(CochlearMes, features = VariableFeatures(object = CochlearMes))
CochlearMes <- FindNeighbors(CochlearMes, dims = 1:10)
CochlearMes <- FindClusters(CochlearMes, resolution = 0.3)
CochlearMes <- RunUMAP(CochlearMes, dims = 1:10)

# assess UMAP
DimPlot(CochlearMes, reduction = "umap")

# decrease resolution to match Rose et al 2023 manuscript
CochlearMes <- FindClusters(CochlearMes, resolution = 0.05)
DimPlot(CochlearMes, reduction = "umap")
CochlearMes <- FindClusters(CochlearMes, resolution = 0.1)
DimPlot(CochlearMes, reduction = "umap")
CochlearMes[[]]

# label the different sub-clusters of mesenchyme 
CochlearMes$res0.1_named <- as.character(CochlearMes$RNA_snn_res.0.1)
CochlearMes$res0.1_named[CochlearMes$res0.1_named %in% c("0", "4")] <- "type1"
CochlearMes$res0.1_named[CochlearMes$res0.1_named == "1"] <- "type3"
CochlearMes$res0.1_named[CochlearMes$res0.1_named == "2"] <- "type4"
CochlearMes$res0.1_named[CochlearMes$res0.1_named %in% c("5", "3")] <- "type2"

# check percent.mt
FeaturePlot(
  object = CochlearMes,
  features = "percent.mt",
  cols = c("blue", "red"),  
  reduction = "umap"        
) +
  ggtitle("UMAP: Mitochondrial Percentage") +
  theme_minimal()

# doublet detection 
sweep.res.list <- paramSweep(CochlearMes, PCs = 1:10, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
nExp <- round(0.075 * nrow(CochlearMes@meta.data))  # set 7.5% are doublet
CochlearMes <- doubletFinder(CochlearMes, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp)
table(CochlearMes@meta.data$DF.classifications_0.25_0.09_75)  

CochlearMes[[]]

# UMAP of singlets/doublets
DimPlot(
  CochlearMes, 
  reduction = "umap", 
  group.by = "DF.classifications_0.25_0.09_345",  
  cols = c("Singlet" = "blue", "Doublet" = "red")
)

# filter out doublet
CochlearMes <- subset(CochlearMes, subset = DF.classifications_0.25_0.09_345 == "Singlet")

# check the number of cell after filtering
dim(CochlearMes)

# gene expression of markers
FeaturePlot(CochlearMes, features = c('Pou3f4', 'Scara3', 'Car3', 'Tgfbi', 'Kcnk2', 'Dlx5', 'Col15a1', 'Dcn', 'Otos', 'Emilin2', 'Ibsp'))

# check umap with new labels
## ----------------------------------------------------------------------------------------------------------------
DimPlot(
  CochlearMes,
  reduction = "umap",
  group.by = c("res0.1_named"),
  combine = FALSE, label.size = 2
)

# we now have the cleaned up data as type 1 - 4 mesenchyme from p7 (Rose et al. 2023)
# we will now combine this p7 data (using Seurat's CCA integration) with only mesenchymal cells from our paper (David et al.)
# so = seurat object of only mesenchymal cells (smart-Seq2)

merged_obj <- merge(x = CochlearMes, y = so)
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)

# merge  & integrate 
merged_obj <- IntegrateLayers(object = merged_obj, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE)

# check percent.mt
FeaturePlot(
  object = merged_obj,
  features = "percent.mt",
  cols = c("blue", "red"), 
  reduction = "umap.cca"        
) +
  ggtitle("UMAP: Mitochondrial Percentage") +
  theme_minimal()

# find nearest neighbor and cluster
merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.cca", dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 0.3, cluster.name = "cca_clusters")
merged_obj[[]]

## ----------------------------------------------------------------------------------------------------------------
unique(merged_obj@meta.data$Age)

# label NA for any CochlearP7 category that doesn't match
merged_obj@meta.data$Age[is.na(merged_obj@meta.data$Age)] <- "CochlearP7"
unique(merged_obj@meta.data$Age)

# run UMAP and plot
merged_obj <- RunUMAP(merged_obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 <- DimPlot(
  merged_obj,
  reduction = "umap.cca",
  group.by = c("cca_clusters"),
  combine = FALSE, label.size = 2
)

DimPlot(
  merged_obj,
  reduction = "umap.cca",
  group.by = c("cca_clusters"),
  combine = FALSE, label.size = 2
)

merged_obj

# look at age distribution
DimPlot(
  merged_obj,
  reduction = "umap.cca",
  group.by = c("Age"),
  combine = FALSE, label.size = 2
)

# create a new column 'TissueType' based on conditions in the 'Age' column
## ----------------------------------------------------------------------------------------------------------------
meta_data <- merged_obj@meta.data

meta_data <- meta_data %>%
  mutate(
    TissueType = case_when(
      Age == "CochlearP7" ~ "Cochlear",
      Age %in% c("p4", "p6") ~ "Utricle",
      TRUE ~ "Other" # For any other values in the 'Age' column
    )
  )

# Assign the updated metadata back to the Seurat object
merged_obj@meta.data <- meta_data

# keep Utricle's color the same
## ----------------------------------------------------------------------------------------------------------------
custom_colors_TissueType<- c(
  "Cochlear" = "#F8766D", 
  "Utricle" = "#3366CC" 
)

# plot UMAP
DimPlot(
  merged_obj,
  reduction = "umap.cca",
  group.by = c("TissueType"),
  cols = custom_colors_TissueType,
  combine = FALSE, label.size = 2
)


# save pdf
plots <- DimPlot(
  merged_obj,
  reduction = "umap.cca",
  group.by = "res0.1_named",
  cols = custom_colors,
  label = TRUE
)

output_file <- "Type1234_Utricle.pdf"

ggsave(
  filename = output_file,
  plot = cowplot::plot_grid(plotlist = plots), 
  width = 8, 
  height = 6, 
  device = "pdf" 
)

# UMAP split by TissueType
## ----------------------------------------------------------------------------------------------------------------
DimPlot(
  merged_obj,
  reduction = "umap.cca",
  group.by = c("TissueType"),
  cols = custom_colors_TissueType,
  combine = FALSE, label.size = 2,
  split.by = "TissueType"
)

unique(merged_obj@meta.data$res0.1_named)

merged_obj@meta.data$res0.1_named[is.na(merged_obj@meta.data$res0.1_named)] <- "Utricle"

DimPlot(
  merged_obj,
  reduction = "umap.cca",
  group.by = c("res0.1_named"),
  combine = FALSE, label.size = 2
)

## ----------------------------------------------------------------------------------------------------------------
# save plot data
plot_obj <- DimPlot(
  merged_obj,
  reduction = "umap.cca",
  group.by = c("res0.1_named"),
  combine = FALSE, label.size = 2
)

## ----------------------------------------------------------------------------------------------------------------

group_colors <- ggplot2::ggplot_build(plot_obj[[1]])$data[[1]][, c("colour", "group")]

# check color
unique_colors <- unique(group_colors)
print(unique_colors)


# define colors
custom_colors <- c(
  "type1" = "#F8766D",  
  "type2" = "#7CAE00", 
  "type3" = "#00BA38",  
  "type4" = "#C77CFF",  
  "Utricle" = "#3366CC" 
)

plot <- DimPlot(
  merged_obj,
  reduction = "umap.cca",
  group.by = "res0.1_named",
  cols = custom_colors,
  split.by = "res0.1_named"
)

output_file <- "Type1234_Utricle_separate2.pdf"

ggsave(
  filename = output_file,
  plot = plot,#cowplot::plot_grid(plotlist = plots), 
  width = 40, 
  height = 6, 
  device = "pdf" 
)


## ----------------------------------------------------------------------------------------------------------------
DimPlot(
  merged_obj,
  reduction = "umap.cca",
  group.by = "res0.1_named",
  cols = custom_colors,
  split.by = "res0.1_named"
)

# Increasing the resolution, will not sub-divide the utricle cells in the integrated space
merged_obj <- FindClusters(merged_obj, resolution = 1.0, cluster.name = "1.0") 
DimPlot(merged_obj, reduction = "umap.cca", label = TRUE)
Idents(merged_obj) <- merged_obj$res0.1_named
merged_obj[[]]
DimPlot(merged_obj, reduction = "umap.cca", group.by = "new_label", label = TRUE) 


# iterate through all Utricle cells and assign them to the closest main cochlear cell type
# Extract UMAP coordinates and meta.data
umap_data <- Embeddings(merged_obj, reduction = "umap.cca")
meta_data <- merged_obj@meta.data

# Add UMAP coordinates to meta.data
meta_data$umapcca_1 <- umap_data[, 1]
meta_data$umapcca_2 <- umap_data[, 2]

# Initialize a new label column with the current cluster identities
meta_data$new_label <- as.character(meta_data$res0.1_named)

# Calculate the center coordinates of each main type
type_centers <- lapply(c("type1", "type2", "type3", "type4"), function(type) {
  # Get cells belonging to the current type
  type_cells <- rownames(meta_data[meta_data$res0.1_named == type, ])
  # Calculate the mean coordinates (center) for the type
  type_center <- colMeans(meta_data[type_cells, c("umapcca_1", "umapcca_2")])
  return(type_center)
})
names(type_centers) <- c("type1", "type2", "type3", "type4")

# Iterate through all Utricle cells and assign them to the closest main cochlear cell type
utricle_cells <- rownames(meta_data[meta_data$res0.1_named == "Utricle", ])
for (utricle_cell in utricle_cells) {
  # Get the UMAP coordinates of the current Utricle cell
  utricle_coord <- meta_data[utricle_cell, c("umapcca_1", "umapcca_2")]
  
  # Compute Euclidean distances to each main type center
  distances <- sapply(type_centers, function(center) {
    sqrt(sum((utricle_coord - center)^2))
  })
  
  # Find the closest main type
  closest_type <- names(which.min(distances))
  
  # Update the label for this Utricle cell to include the closest main type
  #meta_data[utricle_cell, "new_label"] <- paste0(closest_type, "_with_utricle")
  meta_data[utricle_cell, "new_label"] <- paste0("utricle_like_", closest_type)}

# Update the Seurat object with the modified meta.data
merged_obj@meta.data <- meta_data

# Set the new labels as the active identities for the Seurat object
Idents(merged_obj) <- meta_data$new_label

# Check the updated label distribution
table(Idents(merged_obj))

DimPlot(merged_obj, reduction = "umap.cca", group.by = "new_label", label = TRUE) 
DimPlot(merged_obj, reduction = "umap.cca", group.by = "new_label", split.by = "new_label") 

# calculate percentage: Utricle / (Utricle + Type)
# Calculate the total number of cells in each cluster
label_counts <- table(Idents(merged_obj))

# Extract the counts for utricle-related labels
utricle_counts <- label_counts[grep("_with_utricle", names(label_counts))]

# Extract the counts for the main types
type_counts <- label_counts[!grepl("_with_utricle", names(label_counts))]

# Calculate the percentage of Utricle cells in each type
utricle_percentages <- sapply(names(utricle_counts), function(label) {
  # Extract the corresponding type name (e.g., type1 from type1_with_utricle)
  type_name <- sub("_with_utricle", "", label)
  # Calculate percentage: Utricle / (Utricle + Type)
  percentage <- utricle_counts[label] / (utricle_counts[label] + type_counts[type_name]) * 100
  return(percentage)
})

# Create a summary table with the results
utricle_summary <- data.frame(
  Type = names(utricle_percentages),
  Utricle_Count = utricle_counts,
  Type_Count = type_counts[names(utricle_counts)],
  Percentage = utricle_percentages
)

# Display the summary table
print(utricle_summary)


## ----------------------------------------------------------------------------------------------------------------
merged_obj[[]]

# Check for NA values in cluster_labels
sum(is.na(meta_data$cluster_labels))

# Replace NA values in cluster_labels with corresponding values from res0.1_named
meta_data$cluster_labels[is.na(meta_data$cluster_labels)] <- meta_data$res0.1_named[is.na(meta_data$cluster_labels)]

# Verify if NA values are replaced
sum(is.na(meta_data$cluster_labels))

# Update the Seurat object if necessary
merged_obj@meta.data <- meta_data

## ----------------------------------------------------------------------------------------------------------------
DimPlot(merged_obj, reduction = "umap.cca", group.by = "cluster_labels", label = TRUE) 

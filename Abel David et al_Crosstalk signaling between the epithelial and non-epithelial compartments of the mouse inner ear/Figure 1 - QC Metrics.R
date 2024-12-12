library(Seurat)
library(Matrix)
#library(DoubletFinder)
library(magrittr)
library(dplyr)
library(tidyverse)
library(scater)
library(SingleCellExperiment)

library(ggpubr)

# Load Seurat Object
so <- readRDS("undmg.p4p6.RDS")
so.raw <- readRDS("undmg.p4p6.RDS")
dim(so) #52628  1155

# Undamaged cells, treated only with saline (sl)
so <- subset(so, subset = state == "undamaged")
# Idents(so) <- 'state'

# Make SCE Object from Seurat Object
sce <- as.SingleCellExperiment(so)

# Rename Plate
table(sce$plate_name) #p4-dt-wt-mes-plt1 173; p4-dt-wt-mes-plt2 0
sce[,sce$plate_code == "180604-30"]$plate_name<- "p4-dt-wt-mes-plt2"
table(sce$plate_name) #p4-dt-wt-mes-plt1 87; p4-dt-wt-mes-plt2 86

dim(sce)

# Calculate QC Metrics
sce <- addPerCellQC(sce, 
                    subsets = list(
                      Mito = grep("^mt-", rownames(sce)),
                      ERCC = grep("^ERCC-", rownames(sce))
                    )
)

# Hard Code
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

# Adatpive threshold cutoffs
attr(qc.lib2, "thresholds") # lower 106440.9, higher 2671637.6 
attr(qc.nexprs2, "thresholds") # lower 907.1321, higher 8679.7008 

# Aggregate
discard2 <- qc.lib2 | qc.nexprs2 | qc.spike | qc.mito
summary(discard2) # Outliers (aggregate) = 148
sce$outliers <- discard2

pal <- c("black", "#ff2a00")
ggplot2::update_geom_defaults("point",list(size = 3, stroke = 0.5, alpha = 0.7))

# Checking for plate quality 
df <- as.data.frame(cbind(plate_name=as.character(sce$plate_name), Outlier= as.character(sce$outliers)))
df[,1] <- as.factor(df[,1])
df[,2] <- as.factor(df[,2])
pct <- df %>% group_by(plate_name) %>% dplyr::count(Outlier, sort = TRUE) %>% mutate(Percent = n / sum(n)*100)
pdf("plate_percentage.pdf", width = 16, height = 2)
ggplot(pct, aes(x = plate_name, y = Percent, fill = Outlier)) + geom_bar(stat = "identity") + labs(x = "Plt_name", y = "Percentage", fill = "Outlier") + scale_fill_manual(values = pal) + theme_pubr() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "right")
dev.off()
pdf("plate_percentage-noaxis.pdf", width = 16, height = 2)
ggplot(pct, aes(x = plate_name, y = Percent, fill = Outlier)) + geom_bar(stat = "identity") + labs(x = "Plt_name", y = "Percentage", fill = "Outlier") + scale_fill_manual(values = pal) + theme_pubr() + theme(axis.text.x = element_blank(), legend.position = "right")
dev.off()

# Removing plate "p6-dt-wt-mes-plt2" had greater than 25% (nearly 50%) outliers
sce <- sce[,sce$plate_name != "p6-dt-wt-mes-plt2"]

# QC Plots

# Stats
sce[, sce$outliers == TRUE] %>% dim()
sce[, sce$outliers == FALSE] %>% dim()
###############################

# QC violin plots highlighting outliers
pdf("readcounts-celltype.pdf", width = 6, height = 6)
plot.x <- "cell_type"
plot.y <- "nCount_RNA"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + scale_color_manual(values = pal) + theme_pubr() + theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust=0.5))
dev.off()

pdf(".detected-celltype.pdf", width = 6, height = 6)
plot.x <- "cell_type"
plot.y <- "detected"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + scale_color_manual(values = pal) + theme_pubr() + theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust=0.5))
dev.off()

pdf("readcount-celltype.pdf", width = 6, height = 6)
plot.x <- "cell_type"
plot.y <- "detected" # readcunt
color.plot <- "outliers"
title <- "Read Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + scale_color_manual(values = pal) + theme_pubr() + theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust=0.5))
dev.off()

pdf("mitopct-celltype.pdf", width = 6, height = 6)
plot.x <- "cell_type"
plot.y <- "subsets_Mito_percent"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + scale_color_manual(values = pal) + theme_pubr() +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust=0.5))
dev.off()

pdf("sum-vs-mitopct-celltype.pdf", width = 10, height = 6)
plotColData(sce, x = "sum", y = "subsets_Mito_percent", colour_by = "outliers", other_fields="cell_type") + scale_color_manual(values= pal)  + facet_grid(~ cell_type) + theme_pubr()
dev.off()

pdf("sum-genotype-supp.pdf", width = 6, height = 6)
plot.x <- "Genotype"
plot.y <- "sum"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + scale_color_manual(values = pal) + theme_pubr() +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust=0.5))
dev.off()

pdf("sum-age-supp.pdf", width = 6, height = 6)
plot.x <- "Age"
plot.y <- "sum"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + scale_color_manual(values = pal) + theme_pubr() +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust=0.5))
dev.off()

pdf("sum-seqpool-supp.pdf", width = 6, height = 6)
plot.x <- "seq_pool"
plot.y <- "sum"
color.plot <- "outliers"
title <- "Total Count"
plotColData(sce,x= plot.x, y=plot.y, colour_by="outliers") +  scale_y_log10() + scale_color_manual(values = pal) + theme_pubr() +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust=0.5))
dev.off()

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
pdf("plate_position.pdf", width = 10, height = 6)
sce$plate_position <- sce$Cell_position
plotPlatePosition(sce, size_by = "sum",  colour_by = NULL, point_alpha = 1.0) + scale_colour_manual(values= c("black")) + scale_fill_manual(values= c("black"))
dev.off()

# Outliers Sum versus nExpres
plotColData(sce, x = "sum", y="detected", colour_by="outliers") + scale_color_manual(values= pal)

# [X] highlight outliers with color, add pal and the colour_by outliers

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

summary(sce[, sce$cell_type == "mes"]$nCount_RNA)

summary(sce[, sce$cell_type == "mes"]$nFeature_RNA)

summary(sce[, sce$cell_type == "se"]$nCount_RNA)

summary(sce[, sce$cell_type == "se"]$nFeature_RNA)

# Convert SCE Object to Seurat Object
so <- as.Seurat(sce)
dim(so) #25815  1007
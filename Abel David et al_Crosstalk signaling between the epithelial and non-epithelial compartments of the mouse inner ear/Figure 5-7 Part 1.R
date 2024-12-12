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

# Load Seurat Object, Visualize Cluster Groups
so <- readRDS("analyzed_utricle_mes_se_p4p6.rds")


scPalette <- function(n) {
  # colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
  #                 '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  
  colorSpace <- c("#3366CC","#DC3912","#FF9900","#109618","#990099","#0099C6","#DD4477","#64C3B7","#B82E2E","#66AA00","#F46FC7","#94594D")
  
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}

rlang::env_unlock(env = asNamespace('CellChat'))
rlang::env_binding_unlock(env = asNamespace('CellChat'))
assign('scPalette', scPalette, envir = asNamespace('CellChat'))
rlang::env_binding_lock(env = asNamespace('CellChat'))
rlang::env_lock(asNamespace('CellChat'))

# Relevel in descending order
clust.desc.order <- Idents(so) %>% table() %>% as.data.frame() %>% arrange(desc(Freq)) %>% select(-Freq)
colnames(clust.desc.order) <- "clust.desc.order"
clust.desc.order
levels(so) <- clust.desc.order %>% as.matrix()
levels(so) # Confirm reorder

# Extract CellChat input files from Seurat object
assay <- "originalexp"
data.input <- GetAssayData(so, assay = assay, slot = "data")
labels <- Idents(so)
meta <- data.frame(group = labels, row.names = names(labels))

# Create CellChat Object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")

# Check Cell Information
levels(cellchat@idents)
cellchat@idents %>% table()

# Set the ligand-receptor interaction database, mouse
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

## Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

## Use a subset of CellChatDB for cell-cell communication analysub
db.subset <- "Secreted Signaling" # Secreted Signaling, Cell-Cell Contact, NULL
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
## To use all CellChatDB for cell-cell communication anlysis
# CellChatDB.use <- CellChatDB

## Set the use database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
# The number of highly variable ligand-receptor pairs used for signaling inference is 661 

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups

cellchat <- filterCommunication(cellchat, min.cells = 10) # This is is the filter step.

# Extracting the inferred cellular communication network as a dataframe

## All inferred cell-cell comms, slot.name = "netP" to access teh inferred comms at the level of signaling pathways
df.net <- subsetCommunication(cellchat)

# Infer the cell cell signaling at a pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network

cellchat <- aggregateNet(cellchat)

# Visualize the aggregated cell-cell communication network (shows number of interactions and total interaction strength) with a circle plot

groupSize <- as.numeric(table(cellchat@idents))

pdf("./figure-5/aggregate_circle_count_labels.pdf", width = 6, height = 6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()
pdf("./figure-5/aggregate_circle_strength_labels.pdf", width = 6, height = 6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

# No Labels

pdf("./figure-5/aggregate_circle_count_nolabels.pdf", width = 6, height = 6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 0.0000000001)
dev.off()
pdf("./figure-5/aggregate_circle_strength_nolabels.pdf", width = 6, height = 6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", vertex.label.cex = 0.0000000001)
dev.off()

# examine the signaling sent from each cell group, controlling the parameter edge.weight.max so that edge wiehgts can be compared betweeen different networks

# calculate number of significant interactions
extractEnrichedLR(cellchat, signaling = cellchat@netP$pathways, geneLR.return = FALSE) # %>% dim() # 94  1

cellchat@netP$pathways

write_csv(as.data.frame(cellchat@netP$pathways), "./figure-5/cellchat-pathways-enriched.csv")
pdf("./figure-5/individual_circle_labels.pdf", width = 18, height = 18)
mat <- cellchat@net$count
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

# Group Numbers
cellchat@meta$group %>% table()
# 1 - Mesenchymal Cells
# 2 - Type I Hair Cells
# 3 - Transitional Epithelial Cells
# 4 - Supporting Cells
# 5 - Type II Hair Cells
# 6 - Glia
# 7 - Roof Cells (CNMD+)
# 8 - Pericytes
# 9 - Schwann Cells
# 10 - Endothelial Cells
# 11 - Macrophages
# 12 - Melanocytes


# Hiearchy plot
vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
pdf("./figure-5/hierarchy_PTN_labels.tiff", width = 8, height = 6)
pathways.show <- c("PTN")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

pdf("./figure-5/hierarchy_PTN_nolabels.tiff", width = 8, height = 6)
pathways.show <- c("PTN")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy", vertex.label.cex = 0.0000000001)
dev.off()

pdf("./figure-5/hierarchy_TGFb_labels.pdf", width = 8, height = 6)
pathways.show <- c("TGFb")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

pdf("./figure-5/hierarchy_TGFb_nolabels.pdf", width = 8, height = 6)
pathways.show <- c("TGFb")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy", vertex.label.cex = 0.0000000001)
dev.off()

pdf("./figure-5/hierarchy_MK_labels.pdf", width = 8, height = 6)
pathways.show <- c("MK")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

pdf("./figure-5/hierarchy_MK_nolabels.pdf", width = 8, height = 6)
pathways.show <- c("MK")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy", vertex.label.cex = 0.0000000001)
dev.off()

pdf("./figure-5/hierarchy_WNT_labels.pdf", width = 8, height = 6)
pathways.show <- c("WNT")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

pdf("./figure-5/hierarchy_WNT_nolabels.pdf", width = 8, height = 6)
pathways.show <- c("WNT")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy", vertex.label.cex = 0.0000000001)
dev.off()

pdf("./figure-5/contribution_PTN.pdf", width = 4, height = 4)
pathways.show <- c("PTN") 
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

pdf("./figure-5/contribution_TGFb.pdf", width = 4, height = 4)
pathways.show <- c("TGFb") 
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

pdf("./figure-5/contribution_MK.pdf", width = 4, height = 4)
pathways.show <- c("MK")
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
pdf("./figure-5/contribution_WNT.pdf", width = 4, height = 4)
pathways.show <- c("WNT")
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

pdf("./figure-5/contribution_PTN_nolabels.pdf", width = 4, height = 4)
pathways.show <- c("PTN") 
netAnalysis_contribution(cellchat, signaling = pathways.show,  font.size = 0)
dev.off()

pdf("./figure-5/contribution_TGFb_nolabels.pdf", width = 4, height = 4)
pathways.show <- c("TGFb") 
netAnalysis_contribution(cellchat, signaling = pathways.show,  font.size = 0)
dev.off()

pdf("./figure-5/contribution_MK_nolabels.pdf", width = 4, height = 4)
pathways.show <- c("MK")
netAnalysis_contribution(cellchat, signaling = pathways.show,  font.size = 0)
dev.off()
pdf("./figure-5/contribution_WNT_nolabels.pdf", width = 4, height = 4)
pathways.show <- c("WNT")
netAnalysis_contribution(cellchat, signaling = pathways.show, font.size = 0)
dev.off()

pdf("./figure-5/mesenchymal_sender_bubble.pdf", width = 4, height = 8)
netVisual_bubble(cellchat, sources.use = "Mesenchymal Cells", targets.use = c("Type I Hair Cells", "Type II Hair Cells", "Supporting Cells", "Transitional Epithelial Cells", "Roof Cells (CNMD+)"), remove.isolate = FALSE) + ggtitle("Mesenchymal Cell Sender")
dev.off()

pdf("./figure-5/mesenchymal_receiver_bubble.pdf", width = 4, height = 8)
netVisual_bubble(cellchat, targets.use = "Mesenchymal Cells", sources.use = c("Type I Hair Cells", "Type II Hair Cells", "Supporting Cells", "Transitional Epithelial Cells", "Roof Cells (CNMD+)"), remove.isolate = FALSE) + ggtitle("Mesenchymal Cell Receiver")
dev.off()

netVisual_bubble(cellchat, sources.use = "Supporting Cells", targets.use = c("Type I Hair Cells", "Type II Hair Cells", "Transitional Epithelial Cells", "Roof Cells (CNMD+)", "Mesenchymal Cells"), remove.isolate = FALSE) + ggtitle("Supporting Cells Sender")

netVisual_bubble(cellchat, sources.use = "Pericytes", targets.use = c("Type I Hair Cells", "Type II Hair Cells", "Supporting Cells", "Transitional Epithelial Cells", "Roof Cells (CNMD+)", "Mesenchymal Cells"), remove.isolate = FALSE) + ggtitle("Pericytes Sender")

netVisual_bubble(cellchat, sources.use = c( "Mesenchymal Cells", "Supporting Cells", "Transitional Epithelial Cells", "Roof Cells (CNMD+)"), targets.use = c("Type I Hair Cells", "Type II Hair Cells"), remove.isolate = FALSE) + ggtitle("Hair Cells Receiver")

netVisual_bubble(cellchat, sources.use = c("Type I Hair Cells", "Type II Hair Cells"), targets.use = c("Pericytes", "Mesenchymal Cells", "Supporting Cells"), remove.isolate = FALSE) + ggtitle("Type I Hair Cells Targetting Pericytes")

par(mfrow=c(1,1))
netVisual_chord_gene(cellchat, sources.use = "Mesenchymal Cells", targets.use = c("Type I Hair Cells", "Type II Hair Cells", "Supporting Cells", "Transitional Epithelial Cells", "Roof Cells (CNMD+)"), lab.cex = 0.5, legend.pos.y = 30)

par(mfrow=c(1,1))
netVisual_chord_gene(cellchat, sources.use = "Mesenchymal Cells", targets.use = c("Type I Hair Cells", "Type II Hair Cells", "Supporting Cells", "Transitional Epithelial Cells", "Roof Cells (CNMD+)"), slot.name = "netP", legend.pos.y = 10)

plotGeneExpression(cellchat, signaling = "PTN", enriched.only = FALSE) # enriched.only = FALSE shows all signaling genes related to one signaling pathway
plotGeneExpression(cellchat, signaling = "MK", enriched.only = FALSE) # enriched.only = FALSE shows all signaling genes related to one signaling pathway
plotGeneExpression(cellchat, signaling = "WNT", enriched.only = FALSE) # enriched.only = FALSE shows all signaling genes related to one signaling pathway

pdf("./figure-5/PTN - Violin Plots.pdf", width = 11, height = 8.5)
plotGeneExpression(cellchat, signaling = "PTN", enriched.only = FALSE) # enriched.only = FALSE shows all signaling genes related to one signaling pathway
dev.off()

pdf("./figure-5/TGFb - Violin Plots.pdf", width = 11, height = 8.5)
plotGeneExpression(cellchat, signaling = "TGFb", enriched.only = FALSE) # enriched.only = FALSE shows all signaling genes related to one signaling pathway
dev.off()

pdf("./figure-5/MK - Violin Plots.pdf", width = 11, height = 8.5)
plotGeneExpression(cellchat, signaling = "MK", enriched.only = FALSE) # enriched.only = FALSE shows all signaling genes related to one signaling pathway
dev.off()

pdf("./figure-5/WNT - Violin Plots.pdf", width = 11, height = 8.5)
plotGeneExpression(cellchat, signaling = "WNT", enriched.only = FALSE) # enriched.only = FALSE shows all signaling genes related to one signaling pathway
dev.off()

pathways.interest <- cellchat@netP$pathways
pathway <- "PTN"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "MK"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "PSAP"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "ncWNT"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "SPP1"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "SEMA3"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "NT"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "KIT"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "EGF"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "FGF"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "PDGF"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "PROS"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "VEGF"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "IGF"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "MIF"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "GAS"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "GALECTIN"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "CXCL"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "GRN"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "ANGPT"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "TGFb"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "BMP"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "ENHO"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "WNT"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off()

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "VISFATIN"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

pathway <- "PERIOSTIN"
dir.create(pathway)
tiff(filename = paste0("./", pathway, "/", pathway, "-violin.tiff"), width = 11, height = 8.5, units = "in", res = 300)
plotGeneExpression(cellchat, signaling = pathway, enriched.only = FALSE) + ggtitle(pathway) 
dev.off

vertex.receiver = c(2,3,4,5,7) # a numeric vector. 
tiff(filename = paste0("./", pathway, "/", pathway, "-hierarchy.tiff"), res = 300, units = "in", width = 8, height = 6)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathway,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

tiff(filename = paste0("./", pathway, "/", pathway, "-contribution.tiff"), res = 300, units = "in", width = 4, height = 4)
netAnalysis_contribution(cellchat, signaling = pathway)
dev.off()

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups

pdf("./figure-5/network_PTN.pdf", width = 8, height = 4)
pathways.show <- "PTN"
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 3, font.size = 10)
dev.off()

pdf("./figure-5/network_TGFb.pdf", width = 8, height = 4)
pathways.show <- "TGFb"
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 3, font.size = 10)
dev.off()


pdf("./figure-5/network_MK.pdf", width = 8, height = 4)
pathways.show <- "MK"
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 3, font.size = 10)
dev.off()

pdf("./figure-5/network_WNT.pdf", width = 8, height = 4)
pathways.show <- "WNT"
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 3, font.size = 10)
dev.off()

# Signaling  role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("PTN"))
#> Signaling role analysis on the cell-cell communication network from user's input

gg1 + gg2
pdf("./figure-5/signaling role scatter_PTN.pdf", width = 4, height = 4)
netAnalysis_signalingRole_scatter(cellchat, signaling = c("PTN"))
dev.off()

pdf("./figure-5/signaling role scatter_TGFb.pdf", width = 4, height = 4)
netAnalysis_signalingRole_scatter(cellchat, signaling = c("TGFb"))
dev.off()


pdf("./figure-5/signaling role scatter_MK.pdf", width = 4, height = 4)
netAnalysis_signalingRole_scatter(cellchat, signaling = c("MK"))
dev.off()

pdf("./figure-5/signaling role scatter_WNT.pdf", width = 4, height = 4)
netAnalysis_signalingRole_scatter(cellchat, signaling = c("WNT"))
dev.off()

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
par(mfrow=c(1,1))

pdf("./figure-5/outgoing heatmap.pdf", width = 8, height = 8)
print(ht1)
dev.off()

pdf("./figure-5/incoming heatmap.pdf", width = 8, height = 8)
print(ht2)
dev.off()

# Signaling role analysis on the cell-cell communication networks of interest
# ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("PTN"))
# ht

library(NMF)
library(ggalluvial)

# selectK(cellchat, pattern = "outgoing")

nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
pdf("./figure-5/River_outgoing_labels.pdf", width = 12, height = 8)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

pdf("./figure-5/River_outgoing_nolabels.pdf", width = 12, height = 8)
netAnalysis_river(cellchat, pattern = "outgoing", font.size = 0)
dev.off()


# dot plot
pdf("./figure-5/Dot_outgoing.pdf", width = 10, height = 3)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

# selectK(cellchat, pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot

pdf("./figure-5/River_incoming_labels.pdf", width = 12, height = 8)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
pdf("./figure-5/River_incoming_nolabels.pdf", width = 12, height = 8)
netAnalysis_river(cellchat, pattern = "incoming", font.size = 0)
dev.off()
# dot plot
pdf("./figure-5/Dot_incoming.pdf", width = 10, height = 3)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()

set.seed(32)
future_options(seed = TRUE)
options(future.rng.onMisuse="ignore")
source("CellChat_issue167_netClusteringFix.R")

# Functional Similiarity
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
tiff("./figure-5/functional_sim_labels.tiff", res = 300, units = "in", width = 8, height = 6)
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
dev.off()
tiff("./figure-5/functional_sim_nolabels.tiff", res = 300, units = "in", width = 8, height = 6)
netVisual_embedding(cellchat, type = "functional", label.size = 0)
dev.off()

netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2)

# Structural Similarity
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering_GitHub(cellchat, type = "structural")
tiff("./figure-5/structural_sim_labels.tiff", res = 300, units = "in", width = 8, height = 6)
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
dev.off()
tiff("./figure-5/structural_sim_nolabels.tiff", res = 300, units = "in", width = 8, height = 6)
netVisual_embedding(cellchat, type = "structural", label.size = 0)
dev.off()

netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)











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

scPalette <- function(n) {
  # colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
  #                 '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  
  colorSpace <- c("#3366CC","#DC3912","#FF9900","#109618","#990099","#DD4477","#64C3B7","#B82E2E","#66AA00","#F46FC7","#94594D")
  
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
freq <- Idents(so) %>% table() %>% as.data.frame() %>% arrange(desc(Freq))
colnames(freq) <- c("Cell Type", "Freq")
str(freq)
freq %>% 
  mutate(labels = Freq) %>%
  ggplot(aes(x="", y = Freq, fill=`Cell Type`, colour = `Cell Type`)) + geom_col() + coord_polar("y", start=0) + blank_theme + ggtitle("Cell Type Frequency") + labs_pubr() + theme(plot.title = element_text(hjust = 0.5)) + geom_text(aes(label = labels), colour = "black", position = position_stack(vjust = 0.5), show.legend = FALSE) + theme(legend.position="bottom", axis.text.x=element_blank())

# UMAP
plot <- DimPlot(so, reduction = "umap")
LabelClusters(plot = plot, id = "ident") + NoLegend()

selected.cells <- c("Mesenchymal Cells", "Transitional Epithelial Cells", "Type I Hair Cells", "Supporting Cells", "Type II Hair Cells", "Roof Cells (CNMD+)")
so <- subset(so, idents = selected.cells)

plot <- DimPlot(so, reduction = "umap")
LabelClusters(plot = plot, id = "ident") + NoLegend()

# Extract CellChat input files from Seurat object
assay <- "originalexp"

# Extract Data
data.input <- GetAssayData(so, assay = assay, slot = "data")
labels <- Idents(so)
meta <- data.frame(group = labels, row.names = names(labels))

# Create CellChat Objects
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
# db.subset <- "Secreted Signaling" # Secreted Signaling, Cell-Cell Contact, NULL
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
## To use all CellChatDB for cell-cell communication anlysis
# CellChatDB.use <- CellChatDB

## Set the use database in the object
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI network (optional)
#cellchat <- projectData(cellchat, PPI.mouse) # mouse PPI network

cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups

cellchat <- filterCommunication(cellchat, min.cells = 10) # This is is the filter step.

# Extracting the inferred cellular communication network as a dataframe

## All inferred cell-cell comms, slot.name = "netP" to access teh inferred comms at the level of signaling pathways
df.net <- subsetCommunication(cellchat)
## from groups 1 and 2 to cell groups 5 and 5
# df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
## mediate by signaling WNT and TGFb
# df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))


# Infer the cell cell signaling at a pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell-cell communication network

cellchat <- aggregateNet(cellchat)

# Visualize the aggregated cell-cell communication network (shows number of interactions and total interaction strength) with a circle plot

groupSize <- as.numeric(table(cellchat@idents))
#par(mfrow = c(1,2), xpd=TRUE)
pdf("./figure-5-2/aggregate_circle_count_labels.pdf", width = 6, height = 6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")
dev.off()

pdf("./figure-5-2/aggregate_circle_count_nolabels.pdf", width = 6, height = 6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 0.0000000001)
dev.off()

pdf("./figure-5-2/aggregate_circle_count_nolabels small.pdf", width = 6, height = 6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, edge.width.max = 4, label.edge= F, title.name = "Number of interactions", vertex.label.cex = 0.0000000001)
dev.off()

pdf("./figure-5-2/individual_circle_labels no labels.pdf", width = 18, height = 18)
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize/4, weight.scale = T, edge.weight.max = 30, label.edge = TRUE, title.name = rownames(mat)[i])
}
dev.off()

pdf("./figure-5-2/individual_circle_labels.pdf", width = 8, height = 8)
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize/4, weight.scale = T, edge.weight.max = 30, label.edge = TRUE, title.name = rownames(mat)[i])
}
dev.off()

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")

# calculate number of significant interactions
extractEnrichedLR(cellchat, signaling = cellchat@netP$pathways, geneLR.return = FALSE) # %>% dim() # 94  1
# extractEnrichedLR(cellchat, signaling = cellchat@netP$pathways, geneLR.return = TRUE) # %>% dim() # 94  1
cellchat@netP$pathways

pathways.show <- "WNT"
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netAnalysis_contribution(cellchat, signaling = pathways.show)
par(mfrow=c(1,1))
plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = FALSE)

pathways.show <- "PTN"
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netAnalysis_contribution(cellchat, signaling = pathways.show)
par(mfrow=c(1,1))
plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = FALSE)

pathways.show <- "MK"
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netAnalysis_contribution(cellchat, signaling = pathways.show)
par(mfrow=c(1,1))
plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = FALSE)

pathways.show <- "TGFb"
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netAnalysis_contribution(cellchat, signaling = pathways.show)
par(mfrow=c(1,1))
plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = FALSE)

pathways.show <- "ENHO"
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netAnalysis_contribution(cellchat, signaling = pathways.show)
par(mfrow=c(1,1))
plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = FALSE)

netVisual_bubble(cellchat, sources.use = "Mesenchymal", targets.use = c("Type I Hair Cells", "Type II Hair Cells", "Supporting Cells", "Transitional Epithelial Cells", "Roof Cells (CNMD+)"), remove.isolate = FALSE) + ggtitle("Supporting Cells Sender")

# # Network Analysis
par(mfrow=c(1,1))
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
# selectK(cellchat, pattern = "outgoing") # drops just after 3
# 
nPatterns = 3
par(mfrow=c(1,1))
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns) # ERROR: Error in U[, i] : subscript out of bounds
par(mfrow=c(1,1))
netAnalysis_river(cellchat, pattern = "outgoing")
par(mfrow=c(1,1))

# river plot
pdf("./figure-5-2/River_outgoing_labels.pdf", width = 12, height = 8)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()

pdf("./figure-5-2/River_outgoing_nolabels.pdf", width = 12, height = 8)
netAnalysis_river(cellchat, pattern = "outgoing", font.size = 0)
dev.off()


# dot plot
pdf("./figure-5-2/Dot_outgoing.pdf", width = 10, height = 3)
netAnalysis_dot(cellchat, pattern = "outgoing")
dev.off()

# 
# selectK(cellchat, pattern = "incoming") # drops just after 3
nPatterns = 3
par(mfrow=c(1,1))
par(mfrow=c(1,1))

cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns) # ERROR: Error in U[, i] : subscript out of bounds
par(mfrow=c(1,1))
netAnalysis_river(cellchat, pattern = "incoming")

pdf("./figure-5-2/River_incoming_labels.pdf", width = 12, height = 8)
netAnalysis_river(cellchat, pattern = "incoming")
dev.off()
pdf("./figure-5-2/River_incoming_nolabels.pdf", width = 12, height = 8)
netAnalysis_river(cellchat, pattern = "incoming", font.size = 0)
dev.off()
# dot plot
pdf("./figure-5-2/Dot_incoming.pdf", width = 10, height = 3)
netAnalysis_dot(cellchat, pattern = "incoming")
dev.off()




# Hierarchy Plots

# Hiearchy plot
vertex.receiver = c(1,3,4,6) # a numeric vector. 
pdf("./figure-5-2/hierarchy_PTN_labels.pdf", width = 8, height = 6)
pathways.show <- c("PTN")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

pdf("./figure-5-2/hierarchy_PTN_nolabels.pdf", width = 8, height = 6)
pathways.show <- c("PTN")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy", vertex.label.cex = 0.0000000001)
dev.off()

pdf("./figure-5-2/hierarchy_TGFb_labels.pdf", width = 8, height = 6)
pathways.show <- c("TGFb")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()

pdf("./figure-5-2/hierarchy_TGFb_nolabels.pdf", width = 8, height = 6)
pathways.show <- c("TGFb")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy", vertex.label.cex = 0.0000000001)
dev.off()

pdf("./figure-5-2/contribution_PTN.pdf", width = 4, height = 4)
pathways.show <- c("PTN") 
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

pdf("./figure-5-2/contribution_TGFb.pdf", width = 4, height = 4)
pathways.show <- c("TGFb") 
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

pdf("./figure-5-2/contribution_PTN_nolabels.pdf", width = 4, height = 4)
pathways.show <- c("PTN") 
netAnalysis_contribution(cellchat, signaling = pathways.show,  font.size = 0)
dev.off()

pdf("./figure-5-2/contribution_TGFb_nolabels.pdf", width = 4, height = 4)
pathways.show <- c("TGFb") 
netAnalysis_contribution(cellchat, signaling = pathways.show,  font.size = 0)
dev.off()

pdf("./figure-5-2/PTN - Violin Plots enriched.pdf", width = 11, height = 8.5)
plotGeneExpression(cellchat, signaling = "PTN", enriched.only = TRUE) # enriched.only = FALSE shows all signaling genes related to one signaling pathway
dev.off()

pdf("./figure-5-2/TGFb - Violin Plots enriched.pdf", width = 11, height = 8.5)
plotGeneExpression(cellchat, signaling = "TGFb", enriched.only = TRUE) # enriched.only = FALSE shows all signaling genes related to one signaling pathway
dev.off()

pdf("./figure-5-2/PTN - Violin Plots.pdf", width = 11, height = 8.5)
plotGeneExpression(cellchat, signaling = "PTN", enriched.only = FALSE) # enriched.only = FALSE shows all signaling genes related to one signaling pathway
dev.off()

pdf("./figure-5-2/TGFb - Violin Plots.pdf", width = 11, height = 8.5)
plotGeneExpression(cellchat, signaling = "TGFb", enriched.only = FALSE) # enriched.only = FALSE shows all signaling genes related to one signaling pathway
dev.off()

pdf("./figure-5-2/network_PTN.pdf", width = 8, height = 4)
pathways.show <- "PTN"
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 3, font.size = 10)
dev.off()

pdf("./figure-5-2/network_TGFb.pdf", width = 8, height = 4)
pathways.show <- "TGFb"
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 3, font.size = 10)
dev.off()

pdf("./figure-5-2/signaling role scatter_PTN.pdf", width = 4, height = 4)
netAnalysis_signalingRole_scatter(cellchat, signaling = c("PTN"))
dev.off()

pdf("./figure-5-2/signaling role scatter_TGFb.pdf", width = 4, height = 4)
netAnalysis_signalingRole_scatter(cellchat, signaling = c("TGFb"))
dev.off()

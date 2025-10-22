######################################################################################
###Code for Cell communication analysis in adipose
######################################################################################
merged_seurat_colon=readRDS("merged_seurat_colon.rds")


library(CellChat)
table(merged_seurat_colon$celllineage)
unique(merged_seurat_colon@meta.data$group)


groups <- c("ND_KO", "HFD_KO", "HFD_WT")


#standard process for cellchat
cellchat_list <- list()
for (group_name in groups) {
  subset_seurat <- subset(merged_seurat_colon, subset = group == group_name)
  cellchat <- createCellChat(
    object = subset_seurat[["RNA"]]$data,
    meta = subset_seurat@meta.data,
    group.by = "celllineage" 
  )
  cellchat <- setIdent(cellchat, ident.use = "celllineage")
  CellChatDB <- CellChatDB.mouse 
  cellchat@DB <- CellChatDB
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "truncatedMean", trim = 0.1)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05)?
    cellchat_list[[group_name]] <- cellchat
}
saveRDS(cellchat_list,file = "cellchat_list.colon.RDS")

#standard process for cellchat
cellchat.colon$HFD_KO = computeCommunProbPathway(cellchat.colon$HFD_KO, thresh = 0.05)
cellchat.colon$HFD_WT = computeCommunProbPathway(cellchat.colon$HFD_WT, thresh = 0.05)
cellchat.colon$ND_KO = computeCommunProbPathway(cellchat.colon$ND_KO, thresh = 0.05)
cellchat.colon$HFD_KO <- aggregateNet(cellchat.colon$HFD_KO)
cellchat.colon$HFD_WT <- aggregateNet(cellchat.colon$HFD_WT)
cellchat.colon$ND_KO <- aggregateNet(cellchat.colon$ND_KO)
cellchat.colon$HFD_KO <- netAnalysis_computeCentrality(cellchat.colon$HFD_KO, slot.name = "netP")
cellchat.colon$HFD_WT <- netAnalysis_computeCentrality(cellchat.colon$HFD_WT, slot.name = "netP")
cellchat.colon$ND_KO <- netAnalysis_computeCentrality(cellchat.colon$ND_KO, slot.name = "netP")

#save csv results
df.net.HFD_KO <- subsetCommunication(cellchat.colon$HFD_KO)
df.netp.HFD_KO <- subsetCommunication(cellchat.colon$HFD_KO, slot.name = "netP")
write.csv(df.net.HFD_KO,file="df.net.HFD_KO.csv")
write.csv(df.netp.HFD_KO,file="df.netp.HFD_KO.csv")

df.net.HFD_WT <- subsetCommunication(cellchat.colon$HFD_WT)
df.netp.HFD_WT <- subsetCommunication(cellchat.colon$HFD_WT, slot.name = "netP")
write.csv(df.net.HFD_WT,file="df.net.HFD_WT.csv")
write.csv(df.netp.HFD_WT,file="df.netp.HFD_WT.csv")

df.net.ND_KO <- subsetCommunication(cellchat.colon$ND_KO)
df.netp.ND_KO <- subsetCommunication(cellchat.colon$ND_KO, slot.name = "netP")
write.csv(df.net.ND_KO,file="df.net.ND_KO.csv")
write.csv(df.netp.ND_KO,file="df.netp.ND_KO.csv")



##############colon cell chat all in one
object.list <- list(
  HFD_KO = cellchat.colon$HFD_KO,
  HFD_WT = cellchat.colon$HFD_WT,
  ND_KO = cellchat.colon$ND_KO
)


cellchat.merged <- mergeCellChat(
  object.list, 
  add.names = names(object.list),
  cell.prefix = TRUE 
)

identical(
  levels(object.list[[1]]@idents), 
  levels(object.list[[2]]@idents)
)

identical(
  levels(object.list[[1]]@idents), 
  levels(object.list[[3]]@idents)
)

#### Compare LR pairs change in Macrophage 

pdf(file="colon.mac.detailLR.HFDKO.NDKO.pdf",width=8,height=6)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merged, idents.use = "Macrophage",comparison = c(3,1))
gg1 + scale_color_manual(values = rep("black", 4)) +
  scale_shape_manual(values = rep(21, 4))
dev.off()

pdf(file="colon.mac.detailLR.HFDKO.HFDWT.pdf",width=8,height=6)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merged, idents.use = "Macrophage",comparison = c(2,1))
gg1 + scale_color_manual(values = rep("black", 4)) +
  scale_shape_manual(values = rep(21, 4))
dev.off()




##########################################################
##Identify different regulated LRs between groups
##########################################################

#HFD_KO vs ND_KO

pos.dataset = "HFD_KO"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")
# perform differential expression analysis 
cellchat.merged.HFD_KO.ND_KO <- identifyOverExpressedGenes(cellchat.merged.HFD_KO.ND_KO, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

net.HFD_KO.ND_KO <- netMappingDEG(cellchat.merged.HFD_KO.ND_KO, features.name = features.name, variable.all = TRUE)
net.up.HFD_KO.ND_KO <- subsetCommunication(cellchat.merged.HFD_KO.ND_KO, net = net.HFD_KO.ND_KO, datasets = "HFD_KO",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down.HFD_KO.ND_KO <- subsetCommunication(cellchat.merged.HFD_KO.ND_KO, net = net.HFD_KO.ND_KO, datasets = "ND_KO",ligand.logFC = -0.05, receptor.logFC = NULL)


net.HFD_KO.ND_KO$index=paste(net.HFD_KO.ND_KO$source, net.HFD_KO.ND_KO$target, net.HFD_KO.ND_KO$interaction_name,sep = ".")
net.up.HFD_KO.ND_KO$index=paste(net.up.HFD_KO.ND_KO$source, net.up.HFD_KO.ND_KO$target, net.up.HFD_KO.ND_KO$interaction_name,sep = ".")
net.down.HFD_KO.ND_KO$index=paste(net.down.HFD_KO.ND_KO$source, net.down.HFD_KO.ND_KO$target, net.down.HFD_KO.ND_KO$interaction_name,sep = ".")
net.HFD_KO.ND_KO$group="HFD_KO.ND_KO"
net.up.HFD_KO.ND_KO$group="HFD_KO.ND_KO"
net.down.HFD_KO.ND_KO$group="HFD_KO.ND_KO"

write.csv(net.HFD_KO.ND_KO,file="net.HFD_KO.ND_KO.csv")
write.csv(net.up.HFD_KO.ND_KO,file="net.up.HFD_KO.ND_KO.csv")
write.csv(net.down.HFD_KO.ND_KO,file="net.down.HFD_KO.ND_KO.csv")


#HFD_KO vs HFD_WT
pos.dataset = "HFD_KO"
features.name = paste0(pos.dataset, ".merged")
cellchat.merged.HFD_KO.HFD_WT <- identifyOverExpressedGenes(cellchat.merged.HFD_KO.HFD_WT, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

net.HFD_KO.HFD_WT <- netMappingDEG(cellchat.merged.HFD_KO.HFD_WT, features.name = features.name, variable.all = TRUE)
net.up.HFD_KO.HFD_WT <- subsetCommunication(cellchat.merged.HFD_KO.HFD_WT, net = net.HFD_KO.HFD_WT, datasets = "HFD_KO",ligand.logFC = 0.05, receptor.logFC = NULL)
net.down.HFD_KO.HFD_WT <- subsetCommunication(cellchat.merged.HFD_KO.HFD_WT, net = net.HFD_KO.HFD_WT, datasets = "HFD_WT",ligand.logFC = -0.05, receptor.logFC = NULL)


net.HFD_KO.HFD_WT$index=paste(net.HFD_KO.HFD_WT$source, net.HFD_KO.HFD_WT$target, net.HFD_KO.HFD_WT$interaction_name,sep = ".")
net.up.HFD_KO.HFD_WT$index=paste(net.up.HFD_KO.HFD_WT$source, net.up.HFD_KO.HFD_WT$target, net.up.HFD_KO.HFD_WT$interaction_name,sep = ".")
net.down.HFD_KO.HFD_WT$index=paste(net.down.HFD_KO.HFD_WT$source, net.down.HFD_KO.HFD_WT$target, net.down.HFD_KO.HFD_WT$interaction_name,sep = ".")
net.HFD_KO.HFD_WT$group="HFD_KO.HFD_WT"
net.up.HFD_KO.HFD_WT$group="HFD_KO.HFD_WT"
net.down.HFD_KO.HFD_WT$group="HFD_KO.HFD_WT"

write.csv(net.HFD_KO.HFD_WT,file="net.HFD_KO.HFD_WT.csv")
write.csv(net.up.HFD_KO.HFD_WT,file="net.up.HFD_KO.HFD_WT.csv")
write.csv(net.down.HFD_KO.HFD_WT,file="net.down.HFD_KO.HFD_WT.csv")


net.up.HFD_KO.HFD_WT.and.HFD_KO.ND_KO=intersect(net.up.HFD_KO.HFD_WT$index, net.up.HFD_KO.ND_KO$index)
tmp1=net.up.HFD_KO.ND_KO[match(net.up.HFD_KO.HFD_WT.and.HFD_KO.ND_KO, net.up.HFD_KO.ND_KO$index),]
tmp2=net.up.HFD_KO.HFD_WT[match(net.up.HFD_KO.HFD_WT.and.HFD_KO.ND_KO, net.up.HFD_KO.HFD_WT$index),]
net.up.HFD_KO.HFD_WT.and.HFD_KO.ND_KO=cbind(tmp1,tmp2)
write.csv(net.up.HFD_KO.HFD_WT.and.HFD_KO.ND_KO,file="net.up.HFD_KO.HFD_WT.and.HFD_KO.ND_KO.csv")

net.down.HFD_KO.HFD_WT.and.HFD_KO.ND_KO=intersect(net.down.HFD_KO.HFD_WT$index, net.down.HFD_KO.ND_KO$index)
tmp1=net.down.HFD_KO.ND_KO[match(net.down.HFD_KO.HFD_WT.and.HFD_KO.ND_KO, net.down.HFD_KO.ND_KO$index),]
tmp2=net.down.HFD_KO.HFD_WT[match(net.down.HFD_KO.HFD_WT.and.HFD_KO.ND_KO, net.down.HFD_KO.HFD_WT$index),]
net.down.HFD_KO.HFD_WT.and.HFD_KO.ND_KO=cbind(tmp1,tmp2)
write.csv(net.down.HFD_KO.HFD_WT.and.HFD_KO.ND_KO,file="net.down.HFD_KO.HFD_WT.and.HFD_KO.ND_KO.csv")



####################################################################################################################
###########################CXCL network plot
pdf(file="colon.netVisual_aggregate.CXCL.pdf",width=12,height=4)
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

####################################################################################################################
# Cxcl L-R ppair visualize
pairLR.use.up=net.up.HFD_KO.ND_KO[net.up.HFD_KO.ND_KO$pathway_name=="CXCL",]
pairLR.use.up = pairLR.use.up[, "interaction_name", drop = F]

pdf(file="colon.netVisual_bubble.hfdko.ndko.cxcl.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()

pdf(file="colon.netVisual_bubble.hfdko.hfdwt.cxcl.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(2, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()


#return results for visualize the comparison of hfdko.ndko and hfdko.hfdwt in the same figure
tmp <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(2, 3, 1),angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]),return.data = T)

write.csv(tmp$communication,file="cxcl_allgroups.csv")
tmp=tmp$communication
tmp$ygroup=paste(tmp$interaction_name,tmp$dataset,sep="(")
tmp$ygroup=paste(tmp$ygroup,")",sep="")
tmp=tmp[,c("prob","group.names","ygroup","pval")]

library(ggplot2)

tmp$group.names <- factor(tmp$group.names, levels = unique(tmp$group.names))
tmp$ygroup <- factor(tmp$ygroup, levels = unique(tmp$ygroup))

############plot regroup, macrophage first
macrophage_cells <- grep("^Macrophage", unique(tmp$group.names), value = TRUE)

#Construct a new order: "Macrophage" at the beginning + remaining categories
new_order <- c(
  macrophage_cells,  
  setdiff(unique(tmp$group.names), macrophage_cells)  
)

# group into factor
tmp$group.names <- factor(tmp$group.names, levels = new_order)


#####
pdf(file="colon.netVisual_bubble.allgroup.cxcl.all.sigdif.pdf",width=13,height=6)

ggplot(tmp, aes(
  x = group.names, 
  y = ygroup, 
  color = prob, 
  size = pval
)) +
  geom_point(alpha = 0.8, stroke = 1.5) +
  scale_color_gradientn(
    colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
    values = scales::rescale(c(0, 0.2, 0.5, 0.8, 1))
  )+
  scale_size_continuous(
    name = "P-value", 
    range = c(4, 6),
    breaks = c(3)
  ) +
  labs(
    x = "Cell-Cell Interaction", 
    y = "Signaling Pathway Group",
    title = "Cell Interaction Probability and Significance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    legend.position = "right"
  )
dev.off()








library(CellChat)
table(merged_seurat_int$celllineage)
unique(merged_seurat_int@meta.data$group)

# ??ȡ?????б?
groups <- c("ND_KO", "HFD_KO", "HFD_WT")  # ????????????????

# ????CellChat?????б?
cellchat_list <- list()
for (group_name in groups) {
  # ??ȡ?Ӽ?
  subset_seurat <- subset(merged_seurat_int, subset = group == group_name)
  
  # ????CellChat????
  cellchat <- createCellChat(
    object = subset_seurat[["RNA"]]$data,  # ??׼??????????
    meta = subset_seurat@meta.data,
    group.by = "celllineage" 
  )
  
  # ????ϸ?????ͱ?ǩ
  cellchat <- setIdent(cellchat, ident.use = "celllineage")
  
  # ???????ݿ⣨??????С????
  CellChatDB <- CellChatDB.mouse  # ?? CellChatDB.mouse
  cellchat@DB <- CellChatDB
  
  # Ԥ????
  cellchat <- subsetData(cellchat)  # ???????źŻ???
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # ?ƶ?ͨѶ???磨?ؼ????裬??ʱ??
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, type = "truncatedMean", trim = 0.1)
  cellchat <- filterCommunication(cellchat, min.cells = 10)  # ???˵?????ͨѶ
  cellchat <- computeCommunProbPathway(cellchat, thresh = 0.05)
  # ???浽?б?
  cellchat_list[[group_name]] <- cellchat
}
saveRDS(cellchat_list,file = "cellchat_list.int.RDS")


###############################
cellchat.int=readRDS("cellchat_list.int.RDS")
cellchat.int$HFD_KO = computeCommunProbPathway(cellchat.int$HFD_KO, thresh = 0.05)
cellchat.int$HFD_WT = computeCommunProbPathway(cellchat.int$HFD_WT, thresh = 0.05)
cellchat.int$ND_KO = computeCommunProbPathway(cellchat.int$ND_KO, thresh = 0.05)
cellchat.int$HFD_KO <- aggregateNet(cellchat.int$HFD_KO)
cellchat.int$HFD_WT <- aggregateNet(cellchat.int$HFD_WT)
cellchat.int$ND_KO <- aggregateNet(cellchat.int$ND_KO)
cellchat.int$HFD_KO <- netAnalysis_computeCentrality(cellchat.int$HFD_KO, slot.name = "netP")
cellchat.int$HFD_WT <- netAnalysis_computeCentrality(cellchat.int$HFD_WT, slot.name = "netP")
cellchat.int$ND_KO <- netAnalysis_computeCentrality(cellchat.int$ND_KO, slot.name = "netP")


df.net.HFD_KO <- subsetCommunication(cellchat.int$HFD_KO)
df.netp.HFD_KO <- subsetCommunication(cellchat.int$HFD_KO, slot.name = "netP")
write.csv(df.net.HFD_KO,file="df.net.HFD_KO.csv")
write.csv(df.netp.HFD_KO,file="df.netp.HFD_KO.csv")

df.net.HFD_WT <- subsetCommunication(cellchat.int$HFD_WT)
df.netp.HFD_WT <- subsetCommunication(cellchat.int$HFD_WT, slot.name = "netP")
write.csv(df.net.HFD_WT,file="df.net.HFD_WT.csv")
write.csv(df.netp.HFD_WT,file="df.netp.HFD_WT.csv")

df.net.ND_KO <- subsetCommunication(cellchat.int$ND_KO)
df.netp.ND_KO <- subsetCommunication(cellchat.int$ND_KO, slot.name = "netP")
write.csv(df.net.ND_KO,file="df.net.ND_KO.csv")
write.csv(df.netp.ND_KO,file="df.netp.ND_KO.csv")


###
object.list.HFD_KO.ND_KO <- list(
  HFD_KO = cellchat.int$HFD_KO,
  ND_KO = cellchat.int$ND_KO
)


cellchat.merged.HFD_KO.ND_KO <- mergeCellChat(
  object.list.HFD_KO.ND_KO, 
  add.names = names(object.list.HFD_KO.ND_KO),
  cell.prefix = TRUE  # ????ϸ?????ظ?
)

pdf(file="int.pathway.information.flow.HFD_KO.ND_KO.pdf",width=12,height=12)
gg1 <- rankNet(cellchat.merged.HFD_KO.ND_KO, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat.merged.HFD_KO.ND_KO, mode = "comparison",measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()


####cell type information flow
celltype_idents=levels(cellchat.int$HFD_KO@idents)
for (i in celltype_idents) {
  pdf(file=paste0("int.pathway.information.flow.HFD_KO.ND_KO.",i,".source.pdf"),width=12,height=8)
  #  pdf(file="int.pathway.information.flow.HFD_KO.ND_KO.pdf",width=12,height=8)
  gg1 <- rankNet(cellchat.merged.HFD_KO.ND_KO, mode = "comparison", measure = "weight", sources.use = i, targets.use = NULL, stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat.merged.HFD_KO.ND_KO, mode = "comparison", measure = "weight", sources.use = i, targets.use = NULL, stacked = F, do.stat = TRUE)
  print(gg1 + gg2)
  dev.off()
  
}
for (i in celltype_idents) {
  pdf(file=paste0("int.pathway.information.flow.HFD_KO.ND_KO.",i,".target.pdf"),width=12,height=8)
  #  pdf(file="int.pathway.information.flow.HFD_KO.ND_KO.pdf",width=12,height=8)
  gg1 <- rankNet(cellchat.merged.HFD_KO.ND_KO, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = i, stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat.merged.HFD_KO.ND_KO, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = i, stacked = F, do.stat = TRUE)
  print(gg1 + gg2)
  dev.off()
  
}

#####################
########HFD_KO HFD_WT
object.list.HFD_KO.HFD_WT <- list(
  HFD_KO = cellchat.int$HFD_KO,
  HFD_WT = cellchat.int$HFD_WT
)

cellchat.merged.HFD_KO.HFD_WT <- mergeCellChat(
  object.list.HFD_KO.HFD_WT, 
  add.names = names(object.list.HFD_KO.HFD_WT),
  cell.prefix = TRUE
)

#### each LR pairs in each celltye compare
library(CellChat)
pdf(file="int.mac.detailLR.HFDKO.HFDWT.pdf",width=16,height=6)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merged.HFD_KO.HFD_WT, idents.use = "Macrophage")
gg2 <- netAnalysis_signalingChanges_scatter(cellchat.merged.HFD_KO.HFD_WT, idents.use = "T_cell")
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()


pdf(file="int.pathway.information.flow.HFD_KO.HFD_WT.pdf",width=12,height=12)
gg1 <- rankNet(cellchat.merged.HFD_KO.HFD_WT, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat.merged.HFD_KO.HFD_WT, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE)
gg1 + gg2
dev.off()

####each cell type information flow
celltype_idents=levels(cellchat.int$HFD_KO@idents)
for (i in celltype_idents) {
  pdf(file=paste0("int.pathway.information.flow.HFD_KO.HFD_WT.",i,".source.pdf"),width=12,height=8)
  #  pdf(file="int.pathway.information.flow.HFD_KO.HFD_WT.pdf",width=12,height=8)
  gg1 <- rankNet(cellchat.merged.HFD_KO.HFD_WT, mode = "comparison", measure = "weight", sources.use = i, targets.use = NULL, stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat.merged.HFD_KO.HFD_WT, mode = "comparison", measure = "weight", sources.use = i, targets.use = NULL, stacked = F, do.stat = TRUE)
  print(gg1 + gg2)
  dev.off()
  
}
for (i in celltype_idents) {
  pdf(file=paste0("int.pathway.information.flow.HFD_KO.HFD_WT.",i,".target.pdf"),width=12,height=8)
  #  pdf(file="int.pathway.information.flow.HFD_KO.HFD_WT.pdf",width=12,height=8)
  gg1 <- rankNet(cellchat.merged.HFD_KO.HFD_WT, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = i, stacked = T, do.stat = TRUE)
  gg2 <- rankNet(cellchat.merged.HFD_KO.HFD_WT, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = i, stacked = F, do.stat = TRUE)
  print(gg1 + gg2)
  dev.off()
  
}

###############Compare outgoing (or incoming) signaling patterns associated with each cell population
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list.HFD_KO.ND_KO[[i]]@netP$pathways, object.list.HFD_KO.ND_KO[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list.HFD_KO.ND_KO[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list.HFD_KO.ND_KO)[i], width = 5, height = 24)
ht2 = netAnalysis_signalingRole_heatmap(object.list.HFD_KO.ND_KO[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list.HFD_KO.ND_KO)[i+1], width = 5, height = 24)

pdf(file="int.pathway.signaling patterns.HFD_KO.ND.KO.outgoing.pdf",width=12,height=12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list.HFD_KO.ND_KO[[i]]@netP$pathways, object.list.HFD_KO.ND_KO[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list.HFD_KO.ND_KO[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list.HFD_KO.ND_KO)[i], width = 5, height = 24)
ht2 = netAnalysis_signalingRole_heatmap(object.list.HFD_KO.ND_KO[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list.HFD_KO.ND_KO)[i+1], width = 5, height = 24)

pdf(file="int.pathway.signaling patterns.HFD_KO.ND.KO.incoming.pdf",width=12,height=12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


####
i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list.HFD_KO.HFD_WT[[i]]@netP$pathways, object.list.HFD_KO.HFD_WT[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list.HFD_KO.HFD_WT[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list.HFD_KO.HFD_WT)[i], width = 5, height = 24)
ht2 = netAnalysis_signalingRole_heatmap(object.list.HFD_KO.HFD_WT[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list.HFD_KO.HFD_WT)[i+1], width = 5, height = 24)

pdf(file="int.pathway.signaling patterns.HFD_KO.HFD_WT.outgoing.pdf",width=12,height=12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

i = 1
# combining all the identified signaling pathways from different datasets 
pathway.union <- union(object.list.HFD_KO.HFD_WT[[i]]@netP$pathways, object.list.HFD_KO.HFD_WT[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list.HFD_KO.HFD_WT[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list.HFD_KO.HFD_WT)[i], width = 5, height = 24)
ht2 = netAnalysis_signalingRole_heatmap(object.list.HFD_KO.HFD_WT[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list.HFD_KO.HFD_WT)[i+1], width = 5, height = 24)

pdf(file="int.pathway.signaling patterns.HFD_KO.HFD_WT.incoming.pdf",width=12,height=12)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()



####
df.net.HFD_KO <- subsetCommunication(cellchat.int$HFD_KO,thresh = 1)
df.netp.HFD_KO <- subsetCommunication(cellchat.int$HFD_KO, slot.name = "netP",thresh = 1)
write.csv(df.net.HFD_KO,file="df.net.HFD_KO.csv")
write.csv(df.netp.HFD_KO,file="df.netp.HFD_KO.csv")

df.net.HFD_WT <- subsetCommunication(cellchat.int$HFD_WT,thresh = 1)
df.netp.HFD_WT <- subsetCommunication(cellchat.int$HFD_WT, slot.name = "netP",thresh = 1)
write.csv(df.net.HFD_WT,file="df.net.HFD_WT.csv")
write.csv(df.netp.HFD_WT,file="df.netp.HFD_WT.csv")

df.net.ND_KO <- subsetCommunication(cellchat.int$ND_KO,thresh = 1)
df.netp.ND_KO <- subsetCommunication(cellchat.int$ND_KO, slot.name = "netP",thresh = 1)
write.csv(df.net.ND_KO,file="df.net.ND_KO.csv")
write.csv(df.netp.ND_KO,file="df.netp.ND_KO.csv")



##############int cell chat all in one
object.list <- list(
  HFD_KO = cellchat.int$HFD_KO,
  HFD_WT = cellchat.int$HFD_WT,
  ND_KO = cellchat.int$ND_KO
)


cellchat.merged <- mergeCellChat(
  object.list, 
  add.names = names(object.list),
  cell.prefix = TRUE  # ????ϸ?????ظ?
)

identical(
  levels(object.list[[1]]@idents), 
  levels(object.list[[2]]@idents)
)

identical(
  levels(object.list[[1]]@idents), 
  levels(object.list[[3]]@idents)
)


#### each LR pairs in each celltye compare
library(CellChat)
pdf(file="int.mac.detailLR.HFDKO.NDKO.pdf",width=16,height=6)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merged, idents.use = "Macrophage",comparison = c(3,1))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat.merged, idents.use = "T_cell",comparison = c(3,1))
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()

pdf(file="int.mac.detailLR.HFDKO.HFDWT.pdf",width=16,height=6)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merged, idents.use = "Macrophage",comparison = c(2,1))
gg2 <- netAnalysis_signalingChanges_scatter(cellchat.merged, idents.use = "T_cell",comparison = c(2,1))
patchwork::wrap_plots(plots = list(gg1,gg2))
dev.off()

##################one type one color
library(CellChat)
pdf(file="int.mac.detailLR.HFDKO.NDKO1.pdf",width=8,height=6)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merged, idents.use = "Macrophage",comparison = c(3,1))
gg1 + scale_color_manual(values = rep("black", 4)) +
  scale_shape_manual(values = rep(21, 4))
dev.off()

pdf(file="int.mac.detailLR.HFDKO.HFDWT1.pdf",width=8,height=6)
gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merged, idents.use = "Macrophage",comparison = c(2,1))
gg1 + scale_color_manual(values = rep("black", 4)) +
  scale_shape_manual(values = rep(21, 4))
dev.off()


#############number
pdf(file="int.compareInteractions.pdf",width=6,height=4)
gg1=compareInteractions(cellchat.merged,show.legend = F,group=c("HFD_KO","ND_KO","HFD_WT"),group.levels=as.factor(c("HFD_KO","HFD_WT","ND_KO")),color.use=c("#66B99F", "#EC8254","#8C9EC6"))
gg2=compareInteractions(cellchat.merged,show.legend = F,group=c("HFD_KO","ND_KO","HFD_WT"),group.levels=as.factor(c("HFD_KO","HFD_WT","ND_KO")),color.use=c("#66B99F", "#EC8254","#8C9EC6"), measure = "weight")
gg1 + gg2
dev.off()


pdf(file="int.netVisual_diffInteraction.hfdko.ndko.pdf",width=12,height=12)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat.merged, weight.scale = T,comparison = c("ND_KO","HFD_KO"),color.use=mypal_3)
netVisual_diffInteraction(cellchat.merged, weight.scale = T,comparison = c("ND_KO","HFD_KO"),color.use=mypal_3, measure = "weight")
dev.off()

pdf(file="int.netVisual_diffInteraction.hfdko.hfdwt.pdf",width=12,height=12)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat.merged, weight.scale = T,comparison = c("HFD_WT","HFD_KO"),color.use=mypal_3)
netVisual_diffInteraction(cellchat.merged, weight.scale = T,comparison = c("HFD_WT","HFD_KO"),color.use=mypal_3, measure = "weight")
dev.off()


####
# ʾ??ϸ??????
int.groups <- levels(object.list[[1]]@idents)
# ????????3????ɫ
mypal_cellchat.int <- mypal_3[1:length(int.groups)]
# ?????????ֵ???ɫ??��
cols_named <- setNames(mypal_cellchat.int,int.groups)

# ?????б??ｻ??????ʹ??
col_list <- list(int.groups = cols_named)


pdf(file="int.netVisual_diffInteraction.heatmap.hfdko.ndko.pdf",width=12,height=6)
gg1 <- netVisual_heatmap(cellchat.merged,comparison = c("ND_KO","HFD_KO"),color.use = cols_named)
gg2 <- netVisual_heatmap(cellchat.merged,comparison = c("ND_KO","HFD_KO"), measure = "weight",color.use = cols_named)
gg1 + gg2
dev.off()

pdf(file="int.netVisual_diffInteraction.heatmap.hfdko.hfdwt.pdf",width=12,height=6)
gg1 <- netVisual_heatmap(cellchat.merged,comparison = c("HFD_WT","HFD_KO"),color.use = cols_named)
gg2 <- netVisual_heatmap(cellchat.merged,comparison = c("HFD_WT","HFD_KO"), measure = "weight",color.use = cols_named)
gg1 + gg2
dev.off()


###only strength
pdf(file="int.netVisual_diffInteraction.heatmap.hfdko.ndko.strength.pdf",width=7,height=6)
gg2 <- netVisual_heatmap(cellchat.merged,comparison = c("ND_KO","HFD_KO"), measure = "weight",color.use = cols_named)
gg2
dev.off()

pdf(file="int.netVisual_diffInteraction.heatmap.hfdko.hfdwt.strength.pdf",width=7,height=6)
gg2 <- netVisual_heatmap(cellchat.merged,comparison = c("HFD_WT","HFD_KO"), measure = "weight",color.use = cols_named)
gg2
dev.off()


#interactions 
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
pdf(file="int.netVisual_diffInteraction.numbers.circleplot.pdf",width=18,height=6)
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]),color.use = cols_named)
}
dev.off()


weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
pdf(file="int.netVisual_diffInteraction.weights.circleplot.pdf",width=18,height=6)
par(mfrow = c(1,3), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Strength of interactions - ", names(object.list)[i]),color.use = cols_named)
}
dev.off()



tmp=c("BMP","Adenosine","FLRT","PDGF","CD200","LIFR","SELL","NCAM","IFN???II","NECTIN","SIRP","PECAM2","GAS","GALECTIN","CD23","CLEC","TGFb","SELPLG","ncWNT")
tmp1=c("CXCL","CD52","MHC???I","12oxoLTB4","CD39","LAMININ","FN1","CNTN","LCK","CypA","IL6","PD???L1","WNT","Histamine","COLLAGEN")

#############################
#Identify cell populations with significant changes in sending or receiving signals

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax,color.use = cols_named)
  gg[[i]] = gg[[i]] + coord_cartesian(ylim = c(0, 7),xlim = c(0, 7))
}
pdf(file="int.netAnalysis_signalingRole_scatter.pdf",width=12,height=4)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
dev.off()

##################
###Identi ###Identify the signaling changes of specific cell populations
# pdf(file="int.netAnalysis_signalingChanges_scatter.macrophage.pdf",width=8,height=4)
# gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merged, idents.use = "Macrophage", comparison = c("ND_KO","HFD_KO"))
# gg2 <- netAnalysis_signalingChanges_scatter(cellchat.merged, idents.use = "Macrophage",comparison = c("HFD_WT","HFD_KO"))
# patchwork::wrap_plots(plots = list(gg1,gg2))
# dev.off()
# 
# pdf(file="int.netAnalysis_signalingChanges_scatter.macrophage.pdf",width=12,height=6)
# gg1 <- netAnalysis_signalingChanges_scatter(cellchat.merged.HFD_KO.ND_KO, idents.use = "Macrophage")
# gg2 <- netAnalysis_signalingChanges_scatter(cellchat.merged.HFD_KO.HFD_WT, idents.use = "Macrophage")
# patchwork::wrap_plots(plots = list(gg1,gg2))
# )

###########################
###L-R pairs
###########################
###L-R pairs
levels(object.list[[1]]@idents)
#[1] "Ma# > levels(object.list[[1]]@idents)
# [1] "Macrophage" "B_cell"     "T_cell"     "DC"         "NK"        
# [6] "Mast"       "ILC2"       "Neutrophil" "Plasma" ###
levels(cellchat.merged.HFD_KO.ND_KO@meta$datasets)
#[1] "HFD_KO" "ND_KO" 

levels(cellchat.merged@meta$datasets)
#[1] "HFD_KO" "HFD_WT" "ND_KO"

pdf(file="int.netVisual_bubble.hfdko.hfdwt.Macrophage.out.pdf",width=12,height=12)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = 1, targets.use = c(1:9),  comparison = c(2,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = 1, targets.use = c(1:9),  comparison = c(2,1), max.dataset = 2, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


pdf(file="int.netVisual_bubble.hfdko.hfdwt.Macrophage.in.pdf",width=12,height=12)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 1,  comparison = c(2,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 1,  comparison = c(2,1), max.dataset = 2, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

####DC
pdf(file="int.netVisual_bubble.hfdko.hfdwt.DC.out.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = 4, targets.use = c(1:9),  comparison = c(2,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = 4, targets.use = c(1:9),  comparison = c(2,1), max.dataset = 2, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


pdf(file="int.netVisual_bubble.hfdko.hfdwt.DC.in.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 4,  comparison = c(2,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 4,  comparison = c(2,1), max.dataset = 2, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


####B
pdf(file="int.netVisual_bubble.hfdko.hfdwt.B.out.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = 2, targets.use = c(1:9),  comparison = c(2,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = 2, targets.use = c(1:9),  comparison = c(2,1), max.dataset = 2, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


pdf(file="int.netVisual_bubble.hfdko.hfdwt.B.in.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 2,  comparison = c(2,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 2,  comparison = c(2,1), max.dataset = 2, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

####T
pdf(file="int.netVisual_bubble.hfdko.hfdwt.T.out.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = 3, targets.use = c(1:9),  comparison = c(2,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = 3, targets.use = c(1:9),  comparison = c(2,1), max.dataset = 2, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


pdf(file="int.netVisual_bubble.hfdko.hfdwt.T.in.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 3,  comparison = c(2,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 3,  comparison = c(2,1), max.dataset = 2, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

#############################
######hfd ko hfd wt
pdf(file="int.netVisual_bubble.hfdko.ndko.Macrophage.out.pdf",width=12,height=12)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = 1, targets.use = c(1:9),  comparison = c(3,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = 1, targets.use = c(1:9),  comparison = c(3,1), max.dataset = 3, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


pdf(file="int.netVisual_bubble.hfdko.ndko.Macrophage.in.pdf",width=12,height=12)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 1,  comparison = c(3,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 1,  comparison = c(3,1), max.dataset = 3, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

####DC
pdf(file="int.netVisual_bubble.hfdko.ndko.DC.out.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = 4, targets.use = c(1:9),  comparison = c(3,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = 4, targets.use = c(1:9),  comparison = c(3,1), max.dataset = 3, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


pdf(file="int.netVisual_bubble.hfdko.ndko.DC.in.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 4,  comparison = c(3,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 4,  comparison = c(3,1), max.dataset = 3, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


####B
pdf(file="int.netVisual_bubble.hfdko.ndko.B.out.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = 2, targets.use = c(1:9),  comparison = c(3,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = 2, targets.use = c(1:9),  comparison = c(3,1), max.dataset = 3, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


pdf(file="int.netVisual_bubble.hfdko.ndko.B.in.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 2,  comparison = c(3,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 2,  comparison = c(3,1), max.dataset = 3, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

####T
pdf(file="int.netVisual_bubble.hfdko.ndko.T.out.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = 3, targets.use = c(1:9),  comparison = c(3,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = 3, targets.use = c(1:9),  comparison = c(3,1), max.dataset = 3, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


pdf(file="int.netVisual_bubble.hfdko.ndko.T.in.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 3,  comparison = c(3,1),max.dataset = 1, title.name = "Increased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat.merged, sources.use = c(1:9), targets.use = 3,  comparison = c(3,1), max.dataset = 3, title.name = "Decreased signaling in HFD_KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


##########################################################
############
##########################################################
pos.dataset = "HFD_KO"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")
# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 
cellchat.merged.HFD_KO.ND_KO <- identifyOverExpressedGenes(cellchat.merged.HFD_KO.ND_KO, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net.HFD_KO.ND_KO <- netMappingDEG(cellchat.merged.HFD_KO.ND_KO, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up.HFD_KO.ND_KO <- subsetCommunication(cellchat.merged.HFD_KO.ND_KO, net = net.HFD_KO.ND_KO, datasets = "HFD_KO",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
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



pos.dataset = "HFD_KO"
# define a char name used for storing the results of differential expression analysis
features.name = paste0(pos.dataset, ".merged")
# perform differential expression analysis 
# Of note, compared to CellChat version < v2, CellChat v2 now performs an ultra-fast Wilcoxon test using the presto package, which gives smaller values of logFC. Thus we here set a smaller value of thresh.fc compared to the original one (thresh.fc = 0.1). Users can also provide a vector and dataframe of customized DEGs by modifying the cellchat@var.features$LS.merged and cellchat@var.features$LS.merged.info. 
cellchat.merged.HFD_KO.HFD_WT <- identifyOverExpressedGenes(cellchat.merged.HFD_KO.HFD_WT, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 
#> Use the joint cell labels from the merged CellChat object

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net.HFD_KO.HFD_WT <- netMappingDEG(cellchat.merged.HFD_KO.HFD_WT, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up.HFD_KO.HFD_WT <- subsetCommunication(cellchat.merged.HFD_KO.HFD_WT, net = net.HFD_KO.HFD_WT, datasets = "HFD_KO",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in NL, i.e.,downregulated in LS
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


#############################
levels(object.list[[1]]@idents)
#[1] "Macrophage" "B_cell"     "T_cell"     "DC"         "NK"        
#[6] "Mast"       "ILC2"       "Neutrophil" "Plasma" 

pairLR.use.up = net.up.HFD_KO.ND_KO[, "interaction_name", drop = F]
pdf(file="int.netVisual_bubble.hfdko.ndko.Macrophage.out.sigdif.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = 1, targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.down, sources.use = 1, targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

pairLR.use.up = net.up.HFD_KO.ND_KO[, "interaction_name", drop = F]
pdf(file="int.netVisual_bubble.hfdko.ndko.Macrophage.in.sigdif.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = 1, comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.down, sources.use = c(1:9), targets.use = 1 , comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


pairLR.use.up = net.up.HFD_KO.ND_KO[, "interaction_name", drop = F]
pdf(file="int.netVisual_bubble.hfdko.ndko.B.in.sigdif.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = 2, comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.down, sources.use = c(1:9), targets.use = 2 , comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

pairLR.use.up = net.up.HFD_KO.ND_KO[, "interaction_name", drop = F]
pdf(file="int.netVisual_bubble.hfdko.ndko.B.out.sigdif.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = 2, targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.down, sources.use = 2, targets.use = c(1:9) , comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2
dev.off()




pairLR.use.up = net.up.HFD_KO.ND_KO[, "interaction_name", drop = F]
pdf(file="int.netVisual_bubble.hfdko.ndko.T.in.sigdif.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = 3, comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.down, sources.use = c(1:9), targets.use = 3 , comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

pairLR.use.up = net.up.HFD_KO.ND_KO[, "interaction_name", drop = F]
pdf(file="int.netVisual_bubble.hfdko.ndko.T.out.sigdif.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = 3, targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.down, sources.use = 3, targets.use = c(1:9) , comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


###DC
pairLR.use.up = net.up.HFD_KO.ND_KO[, "interaction_name", drop = F]
pdf(file="int.netVisual_bubble.hfdko.ndko.DC.in.sigdif.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = 4, comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.down, sources.use = c(1:9), targets.use = 4 , comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2
dev.off()

pairLR.use.up = net.up.HFD_KO.ND_KO[, "interaction_name", drop = F]
pdf(file="int.netVisual_bubble.hfdko.ndko.DC.out.sigdif.pdf",width=12,height=18)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = 4, targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.down, sources.use = 4, targets.use = c(1:9) , comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
gg1 + gg2
dev.off()


##########
pairLR.use.up=net.up.HFD_KO.ND_KO[net.up.HFD_KO.ND_KO$pathway_name=="CXCL",]
pairLR.use.up = pairLR.use.up[, "interaction_name", drop = F]

pdf(file="int.netVisual_bubble.hfdko.ndko.cxcl.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()

pdf(file="int.netVisual_bubble.hfdko.hfdwt.cxcl.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(2, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()

###
pairLR.use.up=net.up.HFD_KO.ND_KO[net.up.HFD_KO.ND_KO$pathway_name=="TNF",]
pairLR.use.up = pairLR.use.up[, "interaction_name", drop = F]

pdf(file="int.netVisual_bubble.hfdko.ndko.tnf.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()

pdf(file="int.netVisual_bubble.hfdko.hfdwt.tnf.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(2, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()

###########################IL6 dotplot
pairLR.use.up=net.up.HFD_KO.ND_KO[net.up.HFD_KO.ND_KO$pathway_name=="IL6",]
pairLR.use.up = pairLR.use.up[, "interaction_name", drop = F]

pdf(file="int.netVisual_bubble.hfdko.ndko.Il6.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()

pdf(file="int.netVisual_bubble.hfdko.hfdwt.Il6.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(2, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()


###########################IL6 dotplot
pairLR.use.up=net.up.HFD_KO.ND_KO[net.up.HFD_KO.ND_KO$pathway_name=="IL6" | net.up.HFD_KO.ND_KO$pathway_name=="CXCL",]
pairLR.use.up = pairLR.use.up[, "interaction_name", drop = F]

pdf(file="int.netVisual_bubble.hfdko.ndko.cxcl.Il6.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()

pdf(file="int.netVisual_bubble.hfdko.hfdwt.cxcl.Il6.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(2, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()


# gg2 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.down, sources.use = c(1:9), targets.use = c(1:9) , comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
# #> Comparing communications on a merged object
# gg1 + gg2
# )


pairLR.use.up = pairLR.use.up[pairLR.use.up[,1] == "CXCL2_CXCR2", "interaction_name",drop=F]
pdf(file="int.netVisual_bubble.hfdko.ndko.cxcl.all.sigdif.pdf",width=8,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
gg1
dev.off()

pdf(file="int.netVisual_bubble.hfdko.hfdwt.cxcl.all.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(2, 3, 1),angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]),return.data = T)
gg1
dev.off()

pairLR.use.up=net.up.HFD_KO.ND_KO[net.up.HFD_KO.ND_KO$pathway_name=="CXCL",]
pairLR.use.up = pairLR.use.up[, "interaction_name", drop = F]
tmp <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(2, 3, 1),angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]),return.data = T)

write.csv(tmp$communication,file="cxcl_allgroups.csv")
tmp=tmp$communication
tmp$ygroup=paste(tmp$interaction_name,tmp$dataset,sep="(")
tmp$ygroup=paste(tmp$ygroup,")",sep="")
tmp=tmp[,c("prob","group.names","ygroup","pval")]

library(ggplot2)

# 确保 group.names 和 ygroup 是因子，能控制排序
tmp$group.names <- factor(tmp$group.names, levels = unique(tmp$group.names))
tmp$ygroup <- factor(tmp$ygroup, levels = unique(tmp$ygroup))

############plot regroup

macrophage_cells <- grep("^Macrophage", unique(tmp$group.names), value = TRUE)
# 2. 构建新顺序：Macrophage开头 + 剩余类别
new_order <- c(
  macrophage_cells,  
  setdiff(unique(tmp$group.names), macrophage_cells)  # 剩余类别按原始顺序
)

# 3. 转换为因子
tmp$group.names <- factor(tmp$group.names, levels = new_order)
#####
pdf(file="int.netVisual_bubble.allgroup.cxcl.all.sigdif.pdf",width=12,height=6)

ggplot(tmp, aes(
  x = group.names, 
  y = ygroup, 
  color = prob, 
  size = pval
)) +
  geom_point(alpha = 0.8, stroke = 1.5) +  # 气泡：透明度+描边
  scale_color_gradientn(
    colours = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
    values = scales::rescale(c(0, 0.2, 0.5, 0.8, 1))  # 非线性映射
  )+
  scale_size_continuous(
    name = "P-value", 
    range = c(4, 6),             # 大小范围：pval=3对应直径4-10
    breaks = c(3)                 # 图例仅显示pval=3
  ) +
  labs(
    x = "Cell-Cell Interaction", 
    y = "Signaling Pathway Group",
    title = "Cell Interaction Probability and Significance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),  # X轴标签倾斜30°
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    legend.position = "right"
  )
dev.off()








pdf(file="int.netVisual_aggregate.IFN2.pdf",width=12,height=4)
pathways.show <- c("IFN-II") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


pdf(file="int.netVisual_aggregate.TGFb.pdf",width=12,height=4)
pathways.show <- c("TGFb") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pdf(file="int.netVisual_aggregate.CXCL.pdf",width=12,height=4)
pathways.show <- c("CXCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


pdf(file="int.netVisual_aggregate.CCL.pdf",width=12,height=4)
pathways.show <- c("CCL") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


pdf(file="int.netVisual_aggregate.TNF.pdf",width=12,height=4)
pathways.show <- c("TNF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pdf(file="int.netVisual_aggregate.TNFlr.pdf",width=12,height=4)
pathways.show <- c("TNF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show,pairLR.use = "Tnf  - Tnfrsf1a", layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

pdf(file="int.netVisual_aggregate.IL6.pdf",width=12,height=4)
pathways.show <- c("IL6") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()

# pdf(file="int.netVisual_aggregate.IFN2.heatmap.pdf",width=12,height=4)
# pathways.show <- c("IFN-II"), 
# par(mfrow = c(1,2), xpd=TRUE)
# ht <- list()
# for (i in 1:length(object.list)) {
#   ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
# }
# #> Do heatmap based on a single object 
# #> 
# #> Do heatmap based on a single object
# ComplexHeatmap::draw(ht[[1]] + ht[[2]] + ht[[3]], ht_gap = unit(0.5, "cm"))
# dev.off()



# # Chord diagram each LR pairs
# pathways.show <- c("CXCL")
# pdf(file="int.Chord.diagram.cxcl.sigdif.pdf",width=12,height=4)
# par(mfrow = c(1,2), xpd=TRUE)
# for (i in 1:length(object.list)) {
#   netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
# }
# dev.off()


###
pairLR.use.up=net.up.HFD_KO.ND_KO[net.up.HFD_KO.ND_KO$pathway_name=="TNF",]
pairLR.use.up = pairLR.use.up[, "interaction_name", drop = F]

pdf(file="int.netVisual_bubble.hfdko.ndko.tnf.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()

pdf(file="int.netVisual_bubble.hfdko.hfdwt.tnf.sigdif.pdf",width=12,height=4)
gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(2, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
gg1
dev.off()

# pairLR.use.up=net.up.HFD_KO.ND_KO[net.up.HFD_KO.ND_KO$pathway_name=="IL1",]
# pairLR.use.up = pairLR.use.up[, "interaction_name", drop = F]
# 
# pdf(file="int.netVisual_bubble.hfdko.ndko.il1.sigdif.pdf",width=12,height=4)
# gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(3, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
# #> Comparing communications on a merged object
# pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
# gg1
# dev.off()
# 
# pdf(file="int.netVisual_bubble.hfdko.hfdwt.il1.sigdif.pdf",width=12,height=4)
# gg1 <- netVisual_bubble(cellchat.merged, pairLR.use = pairLR.use.up, sources.use = c(1:9), targets.use = c(1:9), comparison = c(2, 1),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
# #> Comparing communications on a merged object
# pairLR.use.down = net.down.HFD_KO.ND_KO[, "interaction_name", drop = F]
# gg1
# dev.off()
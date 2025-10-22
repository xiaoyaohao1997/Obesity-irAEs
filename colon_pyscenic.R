###################################################################################
###################################################################################
######code for pyscenic analysis and visualization
###################################################################################
###################################################################################
loom <- open_loom("./colon.aucell.loom")
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")

regulons_incidMat[1:4,1:4]

regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)

######################## visualization
sub_regulonAUC <- regulonAUC[,match(colnames(merged_seurat_int),colnames(regulonAUC))]
dim(sub_regulonAUC)
merged_seurat_int

identical(colnames(sub_regulonAUC), colnames(merged_seurat_int))

cellClusters <- data.frame(row.names = colnames(merged_seurat_int), 
                           seurat_clusters = as.character(merged_seurat_int$celllineage))
cellTypes <- data.frame(row.names = colnames(merged_seurat_int), 
                        celltype = merged_seurat_int$celllineage)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

auc_mat <- assay(sub_regulonAUC)
regulon_assay <- CreateAssayObject(counts = auc_mat)
merged_seurat_int[["RegulonAUC"]] <- regulon_assay
DefaultAssay(merged_seurat_int) <- "RegulonAUC"


#save the regulon and target list
regulons_targetgene.list=matrix(,ncol = 2,nrow = nrow(regulons_incidMat))
regulons_targetgene.list[,1]=row.names(regulons_incidMat)
for (i in 1:nrow(regulons_incidMat)) {
  regulons_targetgene.list[i,2]=paste(colnames(regulons_incidMat)[regulons_incidMat[i,]!=0],collapse = ",")
}
write.csv(regulons_targetgene.list,file="int.regulons_targetgene.list.csv")



##TFs significant in  hfdko vs ndko
results_list <- list()
meta_data <- merged_seurat_int@meta.data
all_celltypes <- unique(meta_data$celllineage)
regulons.to.plot=rownames(auc_mat)


for (ct in all_celltypes) {
  cells_ct <- rownames(meta_data[meta_data$celllineage == ct, ])
  
  df <- FetchData(merged_seurat_int, vars = c(regulons.to.plot, "group")) %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::filter(cell %in% cells_ct, group %in% c("ND_KO", "HFD_KO"))
  
  df$group <- factor(df$group, levels = c("ND_KO", "HFD_KO"))
  
  for (r in regulons.to.plot) {
    if (length(unique(df$group)) == 2) {
      avg_by_group <- df %>%
        group_by(group) %>%
        summarise(mean_auc = mean(.data[[r]], na.rm = TRUE)) %>%
        pivot_wider(names_from = group, values_from = mean_auc, names_prefix = "mean_")
      
      delta <- avg_by_group$mean_HFD_KO - avg_by_group$mean_ND_KO
      direction <- ifelse(delta > 0, "up", "down")
      
      test <- wilcox.test(df[[r]] ~ df$group)
      
      results_list[[paste(ct, r, sep = "_")]] <- data.frame(
        celltype = ct,
        regulon = r,
        p.value = test$p.value,
        mean_ND_KO = avg_by_group$mean_ND_KO,
        mean_HFD_KO = avg_by_group$mean_HFD_KO,
        delta = delta,
        direction = direction
      )
    }
  }
}

#adjust p
results_ct.hfdko.ndko <- bind_rows(results_list) %>%
  group_by(celltype) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

print(head(results_ct.hfdko.ndko, 10))

#save
write.csv(results_ct.hfdko.ndko,file="int.results_tf.hfdko.ndko.csv")
results_ct.hfdko.ndko=read.csv("int.results_tf.hfdko.ndko.csv")

#add sig label
results_ct.hfdko.ndko$sig="not_sig"
results_ct.hfdko.ndko$sig[results_ct.hfdko.ndko$p.adj<0.05]="sig"
head(results_ct.hfdko.ndko)
results_ct.hfdko.ndko=results_ct.hfdko.ndko[results_ct.hfdko.ndko$p.adj<0.05,]
results_ct.hfdko.ndko$index=paste(results_ct.hfdko.ndko$celltype,results_ct.hfdko.ndko$regulon,results_ct.hfdko.ndko$direction,sep="_")



##TFs significant in  HFD_KO vs HFD_WT
results_list <- list()
meta_data <- merged_seurat_int@meta.data
all_celltypes <- unique(meta_data$celllineage)
regulons.to.plot=rownames(auc_mat)


for (ct in all_celltypes) {
  cells_ct <- rownames(meta_data[meta_data$celllineage == ct, ])
  
  df <- FetchData(merged_seurat_int, vars = c(regulons.to.plot, "group")) %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::filter(cell %in% cells_ct, group %in% c("HFD_WT", "HFD_KO"))
  
  df$group <- factor(df$group, levels = c("HFD_WT", "HFD_KO"))
  
  for (r in regulons.to.plot) {
    if (length(unique(df$group)) == 2) {
      avg_by_group <- df %>%
        group_by(group) %>%
        summarise(mean_auc = mean(.data[[r]], na.rm = TRUE)) %>%
        pivot_wider(names_from = group, values_from = mean_auc, names_prefix = "mean_")
      
      delta <- avg_by_group$mean_HFD_KO - avg_by_group$mean_HFD_WT
      direction <- ifelse(delta > 0, "up", "down")
      
      test <- wilcox.test(df[[r]] ~ df$group)
      
      results_list[[paste(ct, r, sep = "_")]] <- data.frame(
        celltype = ct,
        regulon = r,
        p.value = test$p.value,
        mean_HFD_WT = avg_by_group$mean_HFD_WT,
        mean_HFD_KO = avg_by_group$mean_HFD_KO,
        delta = delta,
        direction = direction
      )
    }
  }
}

#adjust p
results_ct.hfdko.hfdwt <- bind_rows(results_list) %>%
  group_by(celltype) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()


print(head(results_ct.hfdko.hfdwt, 10))

#save
write.csv(results_ct.hfdko.hfdwt,file="int.results_tf.hfdko.hfdwt.csv")
results_ct.hfdko.hfdwt=read.csv("int.results_tf.hfdko.hfdwt.csv")

#add sig label
results_ct.hfdko.hfdwt$sig="not_sig"
results_ct.hfdko.hfdwt$sig[results_ct.hfdko.hfdwt$p.adj<0.05]="sig"
results_ct.hfdko.hfdwt=results_ct.hfdko.hfdwt[results_ct.hfdko.hfdwt$p.adj<0.05,]
head(results_ct.hfdko.hfdwt)
results_ct.hfdko.hfdwt$index=paste(results_ct.hfdko.hfdwt$celltype,results_ct.hfdko.hfdwt$regulon,results_ct.hfdko.hfdwt$direction,sep="_")



#########common in both comparison
commonids=intersect(results_ct.hfdko.ndko$index,results_ct.hfdko.hfdwt$index)
results_ct.hfdko.ndko.common=results_ct.hfdko.ndko[match(commonids, results_ct.hfdko.ndko$index),]
results_ct.hfdko.hfdwt.common=results_ct.hfdko.hfdwt[match(commonids, results_ct.hfdko.hfdwt$index),]
results_ct.hfdko.ndko.and.hfdko.hfdwt.common=cbind(results_ct.hfdko.ndko.common,results_ct.hfdko.hfdwt.common)

#save common results
write.csv(results_ct.hfdko.ndko.and.hfdko.hfdwt.common,file="int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.csv")



#####################################################################################
## visualize results  (Prepare input for Gephi v0.10)
#####################################################################################
int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common=read.csv("../int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.csv")
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

##read the regulon and target list
regulons_targetgene.list=read.csv("./int.regulons_targetgene.list.csv")
regulons_targetgene.list=regulons_targetgene.list[,-1]

int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$target=regulons_targetgene.list[match(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$regulon, regulons_targetgene.list[,1]),2]


int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common=na.omit(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common)
int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up=int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common[int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="Macrophage" & int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up",]

int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2=gsub("\\(\\+\\)","",int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$regulon)

##TF significant 
tfs=c("Ets2(+)","Nfe2l2(+)","Spi1(+)","Pparg(+)","Irf5(+)")
int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up.tfs=int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up[int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up$regulon %in% tfs,]

###dif TFs Ets2(+), Nfe2l2(+)
target.list=list()
for (i in 1:nrow(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up.tfs)) {
  target.list[[i]]=strsplit(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up.tfs$target[i],split = ",")
}


###tf all, each cell each tf and it's target in a data.frame, for visualize
target.list.all=list()
for (i in 1:nrow(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common)) {
  target.list.all[[i]]=strsplit(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$target[i],split = ",")
}
tmp=length(unlist(target.list.all))

target.df.all=int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common[rep(1, tmp),]
target.df.all$target.gene=""
target.df.all$target.gene=unlist(target.list.all)

####replace with true data
index=1
for (i in 1:nrow(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common)) {
  target.df.all[index:(index+length(target.list.all[[i]][[1]])-1),1:24]=int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common[i,1:24]
  index=index+length(target.list.all[[i]][[1]])
}

####Data process for regulon network visualize
head(target.df.all)
target.df.all$source.gene=target.df.all$regulon
target.df.all$source.gene=gsub("\\(\\+\\)","",target.df.all$source.gene)

target.df.all$relation="regulate"

###add regulon and mac (cell types)
tmp=int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common
tmp$target.gene=int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$regulon
tmp$target.gene=gsub("\\(\\+\\)","",tmp$target.gene)
tmp$source.gene=tmp$celltype
tmp$relation="Direct"
target.df.all=rbind(tmp,target.df.all)


length(unique(unlist(target.list.all)))

#########to node and edges
#all edge data with both source and targt and edge type
edges <- target.df.all %>% 
  rename(from = source.gene, to = target.gene, relation = relation)

#all node data with target genes
nodes <- data.frame(
  id = unique(c(target.df.all$source.gene, target.df.all$target.gene)),
  label = unique(c(target.df.all$source.gene, target.df.all$target.gene))
)
edges=edges[,c(1:24,26,25,27)]

edges.mac.up=edges[edges$celltype=="Macrophage" & edges$direction=="up",]
write.table(edges.mac.up[,],file="edges.mac.up.edge.txt",sep="\t",row.names = F,col.names = T)

edges.mac=edges[edges$celltype=="Macrophage",]
which(colnames(edges.mac)=="direction")
tmp=edges.mac[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
tmp[tmp[,3]=="regulate",4]=""

##save edge results
write.csv(tmp,file="int.mac.edge.csv",fileEncoding = "UTF-8",row.names = F)

#node data process
nodes.mac <- data.frame(
  id = unique(c(edges.mac$from, edges.mac$to)),
  label = unique(c(edges.mac$from, edges.mac$to))
)
nodes.mac$type="genes"
nodes.mac$type[1]="Celltype"

###match up tfs in nodes.mac, then update this data
nodes.mac$type[match(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="Macrophage" & int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.mac[,1])]="up"
nodes.mac$type[match(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="Macrophage" & int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.mac[,1])]="down"

#save node results
write.csv(nodes.mac,file="int.mac.node.csv",fileEncoding = "UTF-8",row.names = F)

#####################################################################################
# replicate for bcell
#####################################################################################
edges.bcell=edges[edges$celltype=="B_cell",]
which(colnames(edges.bcell)=="direction")
tmp=edges.bcell[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
tmp[tmp[,3]=="regulate",4]=""
write.csv(tmp,file="int.bcell.edge.csv",fileEncoding = "UTF-8",row.names = F)
nodes.bcell <- data.frame(
  id = unique(c(edges.bcell$from, edges.bcell$to)),
  label = unique(c(edges.bcell$from, edges.bcell$to))
)
nodes.bcell$type="genes"
nodes.bcell$type[1]="Celltype"

nodes.bcell$type[match(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="B_cell" & int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.bcell[,1])]="up"
nodes.bcell$type[match(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="B_cell" & int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.bcell[,1])]="down"
write.csv(nodes.bcell,file="int.bcell.node.csv",fileEncoding = "UTF-8",row.names = F)

#####################################################################################
# replicate for Tcell
#####################################################################################
edges.tcell=edges[edges$celltype=="T_cell",]
which(colnames(edges.tcell)=="direction")
tmp=edges.tcell[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
tmp[tmp[,3]=="regulate",4]=""
write.csv(tmp,file="int.tcell.edge.csv",fileEncoding = "UTF-8",row.names = F)
nodes.tcell <- data.frame(
  id = unique(c(edges.tcell$from, edges.tcell$to)),
  label = unique(c(edges.tcell$from, edges.tcell$to))
)
nodes.tcell$type="genes"
nodes.tcell$type[1]="Celltype"

nodes.tcell$type[match(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="T_cell" & int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.tcell[,1])]="up"
nodes.tcell$type[match(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="T_cell" & int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.tcell[,1])]="down"
write.csv(nodes.tcell,file="int.tcell.node.csv",fileEncoding = "UTF-8",row.names = F)

#####################################################################################
# replicate for DC
#####################################################################################
edges.dc=edges[edges$celltype=="DC",]
which(colnames(edges.dc)=="direction")
tmp=edges.dc[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
tmp[tmp[,3]=="regulate",4]=""
write.csv(tmp,file="int.dc.edge.csv",fileEncoding = "UTF-8",row.names = F)
nodes.dc <- data.frame(
  id = unique(c(edges.dc$from, edges.dc$to)),
  label = unique(c(edges.dc$from, edges.dc$to))
)
nodes.dc$type="genes"
nodes.dc$type[1]="Celltype"

nodes.dc$type[match(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="DC" & int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.dc[,1])]="up"
nodes.dc$type[match(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="DC" & int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.dc[,1])]="down"
write.csv(nodes.dc,file="int.dc.node.csv",fileEncoding = "UTF-8",row.names = F)





#####################################################################################
# vlnplot of TFs and corresponding genes
#####################################################################################
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )+
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2)))+
    stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6,step.increase = 0.1)
  
  return(p)
}
## main function
StackedVlnPlot <- function(obj, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

target.group=c("HFD_KO","HFD_WT","ND_KO")
my_comp=list(c("HFD_KO","ND_KO"),c("HFD_KO","HFD_WT"))
pdf(file="int.regulon.mac.vlnplot.pdf",width=3,height=6)
StackedVlnPlot(subset(merged_seurat_int,subset = group %in% target.group),ident="Macrophage", group.by="group",c("Ets2","Nfe2l2","Irf5"), pt.size=0, cols=c("#FF8080", "#FDCABF","#9B8EB6"))
dev.off()

merged_seurat_int@active.assay="RegulonAUC"
pdf(file="int.regulon.mac.vlnplot.activaty.pdf",width=4,height=6)
p=StackedVlnPlot(subset(merged_seurat_int,subset = group %in% target.group),ident="Macrophage",features=c("Ets2(+)","Nfe2l2(+)","Irf5(+)"),pt.size=0, cols=c("#FF8080", "#FDCABF","#9B8EB6"),group.by="group")+
  ylab("Activaty Level")
print(p)
dev.off()

merged_seurat_int@active.assay="RegulonAUC"
target.group=c("HFD_KO","HFD_WT","ND_KO")
my_comp=list(c("HFD_KO","ND_KO"),c("HFD_KO","HFD_WT"))


pdf(file="int.regulon.mac.vlnplot.activity.ets2.pdf",width=4,height=6)
VlnPlot(subset(merged_seurat_int,subset = group %in% target.group),ident="Macrophage", group.by="group",features=c("Ets2(+)"), pt.size=0, cols=c("#FF8080", "#FDCABF","#9B8EB6"))+
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()


#Visualize Ets2 expression
ETS2_expression <- GetAssayData(merged_seurat_int, assay = "RNA",slot="data")[c("Ets2"), ,drop=F]
merged_seurat_int=AddMetaData(merged_seurat_int,t(ETS2_expression))

Ets2_HFD_KO_mac <- merged_seurat_int@meta.data[which(merged_seurat_int@meta.data$group == "HFD_KO"&merged_seurat_int@meta.data$celllineage =="Macrophage"),]$Ets2
Ets2_HFD_WT_mac <- merged_seurat_int@meta.data[which(merged_seurat_int@meta.data$group == "HFD_WT"&merged_seurat_int@meta.data$celllineage =="Macrophage"),]$Ets2
Ets2_ND_KO_mac <- merged_seurat_int@meta.data[which(merged_seurat_int@meta.data$group == "ND_KO"&merged_seurat_int@meta.data$celllineage =="Macrophage"),]$Ets2
wilcox.test(Ets2_HFD_KO_mac, Ets2_HFD_WT_mac)
wilcox.test(Ets2_HFD_KO_mac, Ets2_ND_KO_mac) 
mean(Ets2_HFD_KO_mac)
mean(Ets2_HFD_WT_mac)
mean(Ets2_ND_KO_mac)
wilcox.test(Ets2_HFD_KO_mac, Ets2_HFD_WT_mac)$p.value
wilcox.test(Ets2_HFD_KO_mac, Ets2_ND_KO_mac)$p.value

#Visualize Ets2 activity in mac
ETS2_activaty <- GetAssayData(merged_seurat_int, assay = "RegulonAUC")[c("Ets2(+)"), ,drop=F]
merged_seurat_int=AddMetaData(merged_seurat_int,t(ETS2_activaty))

Ets2_activaty_HFD_KO_mac <- merged_seurat_int@meta.data[which(merged_seurat_int@meta.data$group == "HFD_KO"&merged_seurat_int@meta.data$celllineage =="Macrophage"),]$`Ets2(+)`
Ets2_activaty_HFD_WT_mac <- merged_seurat_int@meta.data[which(merged_seurat_int@meta.data$group == "HFD_WT"&merged_seurat_int@meta.data$celllineage =="Macrophage"),]$`Ets2(+)`
Ets2_activaty_ND_KO_mac <- merged_seurat_int@meta.data[which(merged_seurat_int@meta.data$group == "ND_KO"&merged_seurat_int@meta.data$celllineage =="Macrophage"),]$`Ets2(+)`
wilcox.test(Ets2_activaty_HFD_KO_mac, Ets2_activaty_HFD_WT_mac)
wilcox.test(Ets2_activaty_HFD_KO_mac, Ets2_activaty_ND_KO_mac) 
mean(Ets2_activaty_HFD_KO_mac)
mean(Ets2_activaty_HFD_WT_mac)
mean(Ets2_activaty_ND_KO_mac)
wilcox.test(Ets2_activaty_HFD_KO_mac, Ets2_activaty_HFD_WT_mac)$p.value
wilcox.test(Ets2_activaty_HFD_KO_mac, Ets2_activaty_ND_KO_mac)$p.value



#####################################################################################
#Dotplot of TFs changes
#####################################################################################

#combine two comparison into a data.frame
#using the delta (difference between comparison)
for (i in unique(int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype)) {
  eval(parse(text = paste0("int.tf.common.",i,"=int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common[int.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype==\"",i,"\",]")))
  eval(parse(text = paste0("int.tf.common.",i,"$delta.mean=(int.tf.common.",i,"$delta + int.tf.common.",i,"$delta.1) / 2")))
  eval(parse(text = paste0("int.tf.common.",i,"$logfoldchange1=log2(int.tf.common.",i,"$mean_HFD_KO /int.tf.common.",i,"$mean_ND_KO)")))
  eval(parse(text = paste0("int.tf.common.",i,"$logfoldchange2=log2(int.tf.common.",i,"$mean_HFD_KO / int.tf.common.",i,"$mean_HFD_WT)")))
  eval(parse(text = paste0("int.tf.common.",i,"$logfc.mean=(int.tf.common.",i,"$logfoldchange1 + int.tf.common.",i,"$logfoldchange2) / 2")))
  
}

library(ggplot2)
library(ggrepel)  
library(dplyr)    

#data for plot
df <- data.frame(
  Regulon = int.tf.common.T_cell$regulon,
  Delta_KO_ND = int.tf.common.T_cell$delta,       
  Delta_KO_WT = int.tf.common.T_cell$delta.1,     
  Direction = factor(int.tf.common.T_cell$direction, levels = c("up", "down")),
  Delta_Mean = int.tf.common.T_cell$delta.mean,     
  logfc.mean = int.tf.common.T_cell$logfc.mean
)

data_range <- range(c(df$Delta_KO_ND, df$Delta_KO_WT))
max_abs <- max(abs(data_range)) * 1.05  
limit <- c(-max_abs, max_abs)

#color mapping
color_mapping <- c("up" = "#AA0000", "down" = "#0066FF")


# top 5 highest and lowest delta.mean 
top_high <- df %>%
  arrange(desc(Delta_Mean)) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "High")

top_low <- df %>%
  arrange(Delta_Mean) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "Low")

extreme_points <- bind_rows(top_high, top_low)

# scatter plot
p <- ggplot(df, aes(x = Delta_KO_ND, y = Delta_KO_WT)) +
  geom_point(aes(color = Direction), size = 3, alpha = 0.7) +
  
  scale_color_manual(values = color_mapping) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # add Axis
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  coord_equal(ratio = 1, xlim = limit, ylim = limit) +
  
  # Axis labels
  labs(
    x = "Delta in HFD.KO vs ND.KO",
    y = "Delta in HFD.KO vs HFD.WT",
    color = "Regulation Direction"
  ) +
  
  # theme optimazation
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   
    panel.grid.minor = element_blank(),   
    panel.background = element_blank(),  
    axis.line = element_line(color = "black")
  )

# Add extreme point labels
p_labeled <- p +
  geom_point(data = extreme_points, aes(color = Direction), size = 4, shape = 1, stroke = 1.5) +
  ggrepel::geom_text_repel(
    data = extreme_points,
    aes(label = paste0(Regulon, "\n(", round(logfc.mean, 2), ")")),
    size = 3.5,
    color = "black",
    box.padding = 0.6,
    max.overlaps = 20,
    segment.color = "gray50"
  )

# save plot
pdf(file="T_cell_tf.Delta.dotplot.pdf", width = 8, height = 8)
print(p_labeled)
dev.off()


#####################################################################################
#Macrophage
df <- data.frame(
  Regulon = int.tf.common.Macrophage$regulon,
  Delta_KO_ND = int.tf.common.Macrophage$delta, 
  Delta_KO_WT = int.tf.common.Macrophage$delta.1,
  Direction = factor(int.tf.common.Macrophage$direction, levels = c("up", "down")),
  Delta_Mean = int.tf.common.Macrophage$delta.mean,
  logfc.mean = int.tf.common.Macrophage$logfc.mean
)
data_range <- range(c(df$Delta_KO_ND, df$Delta_KO_WT))
max_abs <- max(abs(data_range)) * 1.05
limit <- c(-max_abs, max_abs)
color_mapping <- c("up" = "#AA0000", "down" = "#0066FF")
top_high <- df %>%
  arrange(desc(Delta_Mean)) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "High")

top_low <- df %>%
  arrange(Delta_Mean) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "Low")

extreme_points <- bind_rows(top_high, top_low)

p <- ggplot(df, aes(x = Delta_KO_ND, y = Delta_KO_WT)) +
  geom_point(aes(color = Direction), size = 3, alpha = 0.7) +
  
  scale_color_manual(values = color_mapping) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  coord_equal(ratio = 1, xlim = limit, ylim = limit) +
  labs(
    x = "Delta in HFD.KO vs ND.KO",
    y = "Delta in HFD.KO vs HFD.WT",
    color = "Regulation Direction"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )
p_labeled <- p +
  geom_point(data = extreme_points, aes(color = Direction), size = 4, shape = 1, stroke = 1.5) +
  ggrepel::geom_text_repel(
    data = extreme_points,
    aes(label = paste0(Regulon, "\n(", round(logfc.mean, 2), ")")),
    size = 3.5,
    color = "black",
    box.padding = 0.6,
    max.overlaps = 20,
    segment.color = "gray50"
  )


pdf(file="Macrophage_tf.Delta.dotplot.pdf", width = 8, height = 8)
print(p_labeled)
dev.off()


#####################################################################################
#DC

df <- data.frame(
  Regulon = int.tf.common.DC$regulon,
  Delta_KO_ND = int.tf.common.DC$delta,
  Delta_KO_WT = int.tf.common.DC$delta.1,
  Direction = factor(int.tf.common.DC$direction, levels = c("up", "down")),
  Delta_Mean = int.tf.common.DC$delta.mean,
  logfc.mean = int.tf.common.DC$logfc.mean
)

data_range <- range(c(df$Delta_KO_ND, df$Delta_KO_WT))
max_abs <- max(abs(data_range)) * 1.05
limit <- c(-max_abs, max_abs)

color_mapping <- c("up" = "#AA0000", "down" = "#0066FF")

top_high <- df %>%
  arrange(desc(Delta_Mean)) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "High")

top_low <- df %>%
  arrange(Delta_Mean) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "Low")

extreme_points <- bind_rows(top_high, top_low)

p <- ggplot(df, aes(x = Delta_KO_ND, y = Delta_KO_WT)) +
  geom_point(aes(color = Direction), size = 3, alpha = 0.7) +
  
  scale_color_manual(values = color_mapping) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  coord_equal(ratio = 1, xlim = limit, ylim = limit) +
  
  labs(
    x = "Delta in HFD.KO vs ND.KO",
    y = "Delta in HFD.KO vs HFD.WT",
    color = "Regulation Direction"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )

p_labeled <- p +
  geom_point(data = extreme_points, aes(color = Direction), size = 4, shape = 1, stroke = 1.5) +
  ggrepel::geom_text_repel(
    data = extreme_points,
    aes(label = paste0(Regulon, "\n(", round(logfc.mean, 2), ")")),
    size = 3.5,
    color = "black",
    box.padding = 0.6,
    max.overlaps = 20,
    segment.color = "gray50"
  )
pdf(file="DC_tf.Delta.dotplot.pdf", width = 8, height = 8)
print(p_labeled)
dev.off()


#####################################################################################
#B cell

df <- data.frame(
  Regulon = int.tf.common.B_cell$regulon,
  Delta_KO_ND = int.tf.common.B_cell$delta,
  Delta_KO_WT = int.tf.common.B_cell$delta.1,
  Direction = factor(int.tf.common.B_cell$direction, levels = c("up", "down")),
  Delta_Mean = int.tf.common.B_cell$delta.mean,
  logfc.mean = int.tf.common.B_cell$logfc.mean
)

data_range <- range(c(df$Delta_KO_ND, df$Delta_KO_WT))
max_abs <- max(abs(data_range)) * 1.05
limit <- c(-max_abs, max_abs)


color_mapping <- c("up" = "#AA0000", "down" = "#0066FF")

top_high <- df %>%
  arrange(desc(Delta_Mean)) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "High")

top_low <- df %>%
  arrange(Delta_Mean) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "Low")

extreme_points <- bind_rows(top_high, top_low)

p <- ggplot(df, aes(x = Delta_KO_ND, y = Delta_KO_WT)) +
  geom_point(aes(color = Direction), size = 3, alpha = 0.7) +
  scale_color_manual(values = color_mapping) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  coord_equal(ratio = 1, xlim = limit, ylim = limit) +
  
  labs(
    x = "Delta in HFD.KO vs ND.KO",
    y = "Delta in HFD.KO vs HFD.WT",
    color = "Regulation Direction"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )

p_labeled <- p +
  geom_point(data = extreme_points, aes(color = Direction), size = 4, shape = 1, stroke = 1.5) +
  ggrepel::geom_text_repel(
    data = extreme_points,
    aes(label = paste0(Regulon, "\n(", round(logfc.mean, 2), ")")),
    size = 3.5,
    color = "black",
    box.padding = 0.6,
    max.overlaps = 20,
    segment.color = "gray50"
  )


pdf(file="B_cell_tf.Delta.dotplot.pdf", width = 8, height = 8)
print(p_labeled)
dev.off()


##############################################################################################################
#TF numbers across cell lineages
##############################################################################################################

library(dplyr)
library(tidyr)
#extract non-duplicated information
tmp=results_ct.hfdko.ndko.and.hfdko.hfdwt.common[,1:10]
tmp=tmp[!is.na(tmp[,1]),]

#summary the numbers
summary_df.tf.common <- tmp %>%
  group_by(celltype, direction) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    celltype = reorder(celltype, count, FUN = sum)  # Rank by total genes (descending)
  )
summary_df.tf.common

#bar plot
pdf(file="int.pancelltype.overlap.detfs.number.counts.pdf",width=6,height=6)
ggplot(summary_df.tf.common, aes(y = celltype, x = count, fill = direction)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(
    aes(label = count), 
    position = position_stack(vjust = 0.5),
    color = "black",  # Changed from white to black
    size = 4,        # Increased from 3.5
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c("down" = "#99CCFF", "up" = "#FFCCCC"),
    labels = c("Downregulated", "Upregulated")
  ) +
  labs(
    y = "Cell Lineage", 
    x = "Number of TFs", 
    title = "Differentially Regulated TFs by Cell Lineage",
    fill = "Regulation"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 11, color = "black"),  # Larger text for both axes
    axis.title = element_text(size = 12, color = "black"), # Larger axis titles
    axis.text.y = element_text(size = 11, color = "black"), # Specific y-axis text
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13, color = "black"),
    legend.position = "top",
    legend.text = element_text(size = 11, color = "black"),  # Larger legend text
    legend.title = element_text(size = 12, color = "black"), # Larger legend title
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.5)
  )
dev.off()
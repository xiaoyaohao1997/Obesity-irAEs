library(Seurat)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(ggpubr)
library(BiocParallel)
library(harmony) 
library(DoubletFinder)
library(SCP)
library(ggpubr)
library(ggrepel)
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
mypal_1 <- union(pal_npg("nrc", alpha = 0.7)(10),my36colors)

mypal_2<-c('#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#B53E2B','#FDBF6F','#FF7F00','#CAB2D6','#6A3D9A','#FFFF99','#B15928'
           ,"#BBBDDC","#C1E6F3","#FEE1D2","#7DADC6","#F7A073","#E6754B","#AD7F7B")
mypal_3<-union(mypal_2,mypal_1)
mypal_4<-unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))[-c(6,17,18,19)][-c(7)]
mypal_5<-c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
mypal_6<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')



#############################################################################################
#Code for Human PBMC scRNA-seq data
###Human PBMC scRNA-seq data was obtained from GSE285888
#load data
seurat_matrix <- fread("GSM8712033_matrix.txt.gz", sep = "\t", header = TRUE,fill = T,data.table = F) 
seurat.meta=read.csv("GSM8712033_metadata.csv")

#extract severe_irAE and non irAE samples
seurat.meta=seurat.meta[seurat.meta$irAE %in% c("severe_irAE","non"),]

colnames(seurat_matrix)=c("gene",colnames(seurat_matrix)[1:222144])
row.names(seurat_matrix)=seurat_matrix[,1]
seurat_matrix=seurat_matrix[,-1]
seurat_matrix=seurat_matrix[,match(seurat.meta[,1],colnames(seurat_matrix))]

#convert to seurat object
pbmc <- CreateSeuratObject(
  counts = seurat_matrix,     
  project = "pbmc"
)

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#standard process,following the origin article
pbmc <- NormalizeData(pbmc,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)
library(future)
plan("multisession", workers = 4)
options(future.globals.maxSize = 200 * 1024^3)
pbmc <- ScaleData(pbmc,
                  vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                  block.size = 1000)
plan(sequential)
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",
                             nfeatures = 2000)
pbmc <- RunPCA(pbmc, verbose = FALSE)
set.seed(100)
pbmc <- RunHarmony(pbmc,orig.reduction = "pca", new.reduction = "harmony",
                         group.by.vars = "orig.ident",
                         plot_convergence = F,
                         max.iter.harmony = 50)


pbmc <- FindNeighbors(pbmc,
                      reduction = "harmony",
                      dims = 1:30)
pbmc <- RunUMAP(pbmc, dims = 1:30,reduction = "harmony")
pbmc <- FindClusters(pbmc, resolution = 0.2,verbose = TRUE)#16 clusters

#Dimplot
pdf(file="pbmc.dimplot.pdf",width=12,height=12)
DimPlot(pbmc,reduction="umap",raster=F,label=T,cols = mypal_3)
dev.off()

png(file="int.feature.plot.png",width=1600,height=1200)
FeaturePlot(pbmc, features = c("CD3E","CD4","CD8A","MKI67",
                               "NKG7","CD79A","C1QA","JCHAIN",
                               "CD14","CLEC10A","HBA1"),raster=FALSE,min.cutoff =0,ncol=4,reduction = "umap",cols= c("gray90", "red"))
dev.off()

###remove cluster 4 low quality, and then re run the standard process
pbmc=subset(pbmc,idents="4",invert=T)
pbmc <- NormalizeData(pbmc,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)
library(future)
plan("multisession", workers = 4)
options(future.globals.maxSize = 200 * 1024^3)
pbmc <- ScaleData(pbmc,
                  vars.to.regress = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                  block.size = 1000)
plan(sequential)
pbmc <- FindVariableFeatures(pbmc,
                             selection.method = "vst",
                             nfeatures = 2000)
pbmc <- RunPCA(pbmc, npcs = 50, verbose = FALSE)

set.seed(100)
pbmc <- RunHarmony(pbmc,orig.reduction = "pca", new.reduction = "harmony",
                   group.by.vars = "orig.ident",
                   plot_convergence = F,
                   max.iter.harmony = 50)


pbmc <- FindNeighbors(pbmc,
                      reduction = "harmony",
                      dims = 1:30)
pbmc <- RunUMAP(pbmc, dims = 1:30,reduction = "harmony")

pbmc <- FindClusters(pbmc, resolution = 0.25,verbose = TRUE)

pdf(file="pbmc.dimplot.afterfilter.rs025.pdf",width=12,height=12)
DimPlot(pbmc,reduction="umap",raster=F,label=T,cols = mypal_3)
dev.off()

library(presto)
library(future)

#Markers 
plan(multisession, workers=4)
markers_pbmc <- FindAllMarkers(pbmc,only.pos = TRUE)
write.csv(markers_pbmc,"./markers_pbmc.afterfilter.res025.csv")
plan(sequential)

#Markers dotplot
png(file="PBMC.dotplot.pdf",width=800,height=800)
DotPlot(pbmc, features = c("CD3E","CD4","CD8A","MKI67",
                           "NKG7","CD79A","C1QA","JCHAIN",
                           "CD14","CLEC10A","HBA1"))
dev.off()


#cell annotation
pbmc <- RenameIdents(pbmc, 
                     `0` = "Monocyte",
                     `1` = "NK",
                     `2` = "CD8T",
                     `3` = "CD4T_INPP4B",
                     `4` = "CD4T_IL7R",
                     `5` = "B",
                     `6` = "Macrophage",
                     `7` = "RBC",
                     `8` = "CD4T_RTKN2",
                     `9` = "DC",
                     `10` = "CD8_pro",
                     `11` = "Platelet",
                     `12` = "pDC" )
pbmc$celltype=pbmc@active.ident
pbmc$celllineage=pbmc@active.ident
pbmc$celllineage<-gsub("CD4T_RTKN2", "CD4T", pbmc$celllineage)
pbmc$celllineage<-gsub("CD4T_INPP4B", "CD4T", pbmc$celllineage)
pbmc$celllineage<-gsub("CD4T_IL7R", "CD4T", pbmc$celllineage)
pbmc$celllineage<-gsub("CD8_pro", "CD8T", pbmc$celllineage)
table(pbmc$celllineage)


#Dimplot
library(SCP)
pdf(file="pbmc.CellDimPlot.celllineage.pdf",width=8,height=7)
CellDimPlot(pbmc, group.by = "celllineage",palcolor = mypal_4, reduction = "UMAP",raster=F)
dev.off()


#Cell marker Featureplot
pdf(file="pbmc.feature.plot1.pdf",width=17.5,height=12)
FeaturePlot(pbmc, features = c("CD3E","CD4","CD8A","MKI67",
                               "NKG7","CD79A","C1QA","JCHAIN",
                               "CD14","CLEC10A","HBA1","PPBP"),raster=F,min.cutoff =0,ncol=4,reduction = "umap",cols= c("gray90", "red"))
dev.off()

#Svae loom data for aucell
library(SeuratDisk)
SaveLoom(
  object = pbmc,
  filename = "pbmc.loom",
  overwrite = F
)

##########################################################################################################
#AUCell for AMPK
library(AUCell)
library(GSEABase)
library(msigdbr)
library(KEGGREST)
expr_matrix <- as.matrix(GetAssayData(pbmc, assay = "RNA", slot = "data"))
pathway_id <- "hsa04152"
pathway_info <- keggGet(pathway_id)
genes <- pathway_info[[1]]$GENE
gene_symbols <- genes[seq(2, length(genes), by = 2)]  
gene_symbols <- gsub(";.*", "", gene_symbols)  
gene_symbols <- trimws(gene_symbols)           
print(gene_symbols)
ampk_pathway <- list("hsa04152_AMPK_SIGNALING_PATHWAY" = gene_symbols)
ampk_geneset <- GeneSet(
  geneIds = ampk_pathway[[1]], 
  geneIdType = SymbolIdentifier(), 
  setName = names(ampk_pathway)[1]
)
geneSets <- GeneSetCollection(ampk_geneset)  

set.seed(100)
pdf(file="AUCell_buildRankings.pdf",width=9,height = 6)
cells_rankings <- AUCell_buildRankings(
  expr_matrix, 
  plotStats = T, 
  nCores = 16)
dev.off()

set.seed(100)
cells_AUC <- AUCell_calcAUC(
  geneSets, 
  cells_rankings,
  aucMaxRank = 1227 #Using AUCell_buildRankings
)

ampk_scores <- getAUC(cells_AUC)["hsa04152_AMPK_SIGNALING_PATHWAY", ]
pbmc$AMPK_AUC <- ampk_scores


library(ggpubr)
my_comp = list(c("non","severe_irAE"))

#AMPK Vlnplot Macrophage
pdf(file="pbmc.aucell.boxplot.irAEgroup.mac.pdf",width=3,height=6)
VlnPlot(pbmc,idents="Macrophage",group.by="irAE",c("AMPK_AUC"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

#AMPK Vlnplot All cells
pdf(file="pbmc.aucell.boxplot.irAEgroup.all.pdf",width=3,height=6)
VlnPlot(pbmc,group.by="irAE",c("AMPK_AUC"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

#AMPK Vlnplot CD4T
pdf(file="pbmc.aucell.boxplot.irAEgroup.CD4T.sig.1227.pdf",width=3,height=6)
VlnPlot(pbmc,idents="CD4T",group.by="irAE",c("AMPK_AUC"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

#AMPK Vlnplot CD8T
pdf(file="pbmc.aucell.boxplot.irAEgroup.CD8T.sig.1227.pdf",width=3,height=6)
VlnPlot(pbmc,idents="CD8T",group.by="irAE",c("AMPK_AUC"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

#############################################################################################
#ETS2 expression comparison
pdf(file="pbmc.ETS2.boxplot.irAEgroup.allcells.pdf",width=3,height=6)
VlnPlot(pbmc,group.by="irAE",c("ETS2"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

pdf(file="pbmc.ETS2.boxplot.irAEgroup.mac.sig.pdf",width=3,height=6)
VlnPlot(pbmc,idents="Macrophage",group.by="irAE",c("ETS2"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

pdf(file="pbmc.ETS2.boxplot.irAEgroup.CD8T.sig.pdf",width=3,height=6)
VlnPlot(pbmc,idents="CD8T",group.by="irAE",c("ETS2"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

pdf(file="pbmc.ETS2.boxplot.irAEgroup.CD4T.sig.pdf",width=3,height=6)
VlnPlot(pbmc,idents="CD4T",group.by="irAE",c("ETS2"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()


#####wilcox test
###################################################
AMPK_irae <- pbmc@meta.data[which(pbmc@meta.data$irAE == "severe_irAE"),]$AMPK_AUC
AMPK_non <- pbmc@meta.data[which(pbmc@meta.data$irAE == "non"),]$AMPK_AUC
wilcox.test(AMPK_irae, AMPK_non)
mean(AMPK_irae)
mean(AMPK_non)
wilcox.test(AMPK_irae, AMPK_non)$p.value

AMPK_irae_mac <- pbmc.mac@meta.data[which(pbmc.mac@meta.data$irAE == "severe_irAE"),]$AMPK_AUC
AMPK_non_mac <- pbmc.mac@meta.data[which(pbmc.mac@meta.data$irAE == "non"),]$AMPK_AUC
wilcox.test(AMPK_irae_mac, AMPK_non_mac)
mean(AMPK_irae_mac)
mean(AMPK_non_mac)


AMPK_irae_cd4 <- pbmc@meta.data[which(pbmc@meta.data$irAE == "severe_irAE" & pbmc@meta.data$celllineage == "CD4T"),]$AMPK_AUC
AMPK_non_cd4 <- pbmc@meta.data[which(pbmc@meta.data$irAE == "non"& pbmc@meta.data$celllineage == "CD4T"),]$AMPK_AUC
wilcox.test(AMPK_irae_cd4, AMPK_non_cd4)
mean(AMPK_irae_cd4)
mean(AMPK_non_cd4)
wilcox.test(AMPK_irae_cd4, AMPK_non_cd4)$p.value


AMPK_irae_cd8 <- pbmc@meta.data[which(pbmc@meta.data$irAE == "severe_irAE" & pbmc@meta.data$celllineage == "CD8T"),]$AMPK_AUC
AMPK_non_cd8 <- pbmc@meta.data[which(pbmc@meta.data$irAE == "non"& pbmc@meta.data$celllineage == "CD8T"),]$AMPK_AUC
wilcox.test(AMPK_irae_cd8, AMPK_non_cd8)
mean(AMPK_irae_cd8)
mean(AMPK_non_cd8)
wilcox.test(AMPK_irae_cd8, AMPK_non_cd8)$p.value


ETS2_expression <- GetAssayData(pbmc, assay = "RNA",slot="data")[c("ETS2"), ,drop=F]
pbmc=AddMetaData(pbmc,t(ETS2_expression))

Ets2_irae_mac <- pbmc@meta.data[which(pbmc@meta.data$irAE == "severe_irAE"&pbmc@meta.data$celllineage =="Macrophage"),]$ETS2
Ets2_non_ac <- pbmc@meta.data[which(pbmc@meta.data$irAE == "non"&pbmc@meta.data$celllineage =="Macrophage"),]$ETS2
wilcox.test(Ets2_irae_mac, Ets2_non_mac)
mean(Ets2_irae_mac)
mean(Ets2_non_mac)
wilcox.test(Ets2_irae_mac, Ets2_non_mac)$p.value



#########################################################################################
################# Violin plot for cytokines
modify_vlnplot <- function(obj, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
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


pdf(file="pbmc.pancelltype.group.stack.vlnplot.ccls.macrophage.pdf",width=4,height=16)
StackedVlnPlot(pbmc,ident="Macrophage", group.by="irAE",c("IL1B","CXCL2","CCL3"), pt.size=0, cols=c("#66B99F", "#EC8254","#8C9EC6", "#ADCE6D"))
dev.off()

pdf(file="pbmc.dotplot.ccls.macrophage.pdf",width=12,height=4)
DotPlot(pbmc,idents="Macrophage", group.by="irAE",c("IL1B","CXCL2","CCL3"))
dev.off()



####################################################################################################################################################
#######DIF GENES
#################Identify DEGs 
# Generate a table of cell counts per group and lineage
cell_counts <- pbmc@meta.data %>%
  group_by(celllineage, irAE) %>%
  summarise(n_cells = n()) %>%
  tidyr::pivot_wider(names_from = irAE, values_from = n_cells, values_fill = 0)
print(cell_counts)
# Identify valid lineages (where both groups have at least 3 cells)
valid_lineages <- cell_counts %>%
  filter(severe_irAE >= 3 & non >= 3) %>%
  pull(celllineage)
print(valid_lineages)


plan(multisession, workers=6)
deg_results_iraevsnon <- list()
for (lineage in valid_lineages) {
  subset_cells <- subset(pbmc, celllineage == lineage)
  degs <- FindMarkers(
    subset_cells,
    ident.1 = "severe_irAE",
    ident.2 = "non",
    group.by = "irAE",
    test.use = "wilcox",  
    min.pct = 0.1,
    logfc.threshold = 0.25)
  degs$genesymbol<-row.names(degs)
  degs$celllineage <- lineage
  deg_results_iraevsnon[[lineage]] <- degs
}
plan(sequential)

degs_iraevsnon <- bind_rows(deg_results_iraevsnon, .id = "celllineage")
degs_iraevsnon$comparison<-rep("iraevsnon",nrow(degs_iraevsnon))

#extract degs
degssig_iraevsnon<-degs_iraevsnon[which(degs_iraevsnon$p_val_adj<0.05&abs(degs_iraevsnon$avg_log2FC)>0.25),]
degssig_iraevsnon <- degssig_iraevsnon %>%
  mutate(direction = case_when(
    avg_log2FC > 0~ "up", 
    avg_log2FC < 0~ "down"))
degssig_iraevsnon <- degssig_iraevsnon %>%
  mutate(direction = factor(direction, levels = c("up", "down"))) %>%
  arrange(celllineage, direction)
head(degssig_iraevsnon)
table(degssig_iraevsnon$celllineage)

#add cell lineage information
degssig_iraevsnon$celllineage<-factor(degssig_iraevsnon$celllineage,levels=degssig_iraevsnon$celllineage)
degssig_iraevsnon <- degssig_iraevsnon %>%
  mutate(log10p_val_adj = -log10(replace(p_val_adj, p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm = TRUE)))
  )
head(degssig_iraevsnon)
write.csv(degssig_iraevsnon,file="degssig_iraevsnon.csv")

top10 <- degssig_iraevsnon %>% 
  group_by(celllineage) %>%
  top_n(10, abs(avg_log2FC)) %>%
  ungroup()

######################################################################################################
#overall degs using all cells
######################################################################################################
degs.overall <- FindMarkers(
  pbmc,
  ident.1 = "severe_irAE",
  ident.2 = "non",
  group.by = "irAE",
  test.use = "wilcox",  
  min.pct = 0.1,
  logfc.threshold = 0.25)
degs.overall$genesymbol<-row.names(degs.overall)

degs.overall$comparison<-rep("iraevsnon",nrow(degs.overall))

degs.overall.sig<-degs.overall[which(degs.overall$p_val_adj<0.05&abs(degs.overall$avg_log2FC)>0.25),]
degs.overall.sig <- degs.overall.sig %>%
  mutate(direction = case_when(
    avg_log2FC > 0~ "up", 
    avg_log2FC < 0~ "down"))
degs.overall.sig <- degs.overall.sig %>%
  mutate(direction = factor(direction, levels = c("up", "down"))) %>%
  arrange(direction)
head(degs.overall.sig)

degs.overall.sig <- degs.overall.sig %>%
  mutate(log10p_val_adj = -log10(replace(p_val_adj, p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm = TRUE)))
  )


####cytokines visualize
pdf(file="pbmc.overall.degs.iraenon.volcano.pdf",width=5,height=5)
ggplot(degs.overall.sig, 
       aes(avg_log2FC, log10p_val_adj, color = direction)) +
  geom_point(size = 1) +
  theme_bw() +
  ggtitle("Irae vs non") +
  xlim(c(-4, 4)) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dotted") +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dotted") +
  scale_color_manual(values = c("#FF3333", "#0066CC")) +

  geom_point(data = degs.overall.sig %>% 
               filter(genesymbol %in% c("IL1B","CXCL2","CCL3","CCL4","NFKB2")),
             size = 3, color = "black") + 

  geom_point(data = degs.overall.sig %>% 
               filter(genesymbol %in% c("IL1B","CXCL2","CCL3","CCL4","NFKB2")),
             size = 2) + 

  geom_text_repel(data = degs.overall.sig %>% 
                    filter(genesymbol %in% c("IL1B","CXCL2","CCL3","CCL4","NFKB2")),
                  aes(label = genesymbol),
                  size = 4,
                  col = 'black',
                  box.padding = 0.5,  
                  max.overlaps = Inf) 
dev.off()


#####################################################################################
# degs in Macrophage
#####################################################################################
library(ggrepel)
library(ggpubr)
degssig_iraevsnon_mac=degssig_iraevsnon[degssig_iraevsnon$celllineage=="Macrophage",]
pdf(file="pbmc.mac.degs.iraenon.volcano.pdf",width=5,height=5)
ggplot(degssig_iraevsnon_mac, 
       aes(avg_log2FC, log10p_val_adj, color = direction)) +
  geom_point(size = 1) +
  theme_bw() +
  ggtitle("Irae vs non") +
  xlim(c(-4, 4)) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dotted") +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dotted") +
  scale_color_manual(values = c("#FF3333", "#0066CC")) +
 
  geom_point(data = degssig_iraevsnon_mac %>% 
               filter(genesymbol %in% c("IL1B","CXCL2","CCL3")),
             size = 3, color = "black") + 

  geom_point(data = degssig_iraevsnon_mac %>% 
               filter(genesymbol %in% c("IL1B","CXCL2","CCL3")),
             size = 2) + 

  geom_text_repel(data = degssig_iraevsnon_mac %>% 
                    filter(genesymbol %in% c("IL1B","CXCL2","CCL3")),
                  aes(label = genesymbol),
                  size = 4,
                  col = 'black',
                  box.padding = 0.5,  
                  max.overlaps = Inf)  
dev.off()


############################################################################################################
##metascape results of mac up-regulated genes 
############################################################################################################
tmp=readxl::read_xlsx("./Mac_up/metascape_result.xlsx",sheet = 2)
tmp=as.data.frame(tmp)
tmp=tmp[grep("Summary",tmp[,1]),]
tmp$Pathways=paste(tmp$Term,tmp$Description,sep=": ")
tmp$neglogP=-tmp$LogP
tmp$neglogQ=-tmp$`Log(q-value)`
tmp=head(tmp,20)
head(tmp)

ggplot(data=tmp, aes(x=reorder(Pathways,neglogQ), y=neglogQ)) + #reorder 
  geom_bar(stat="identity", width=0.8,fill="#B53E2B") + coord_flip() + 
  #  scale_x_discrete(labels=labels) +
  xlab("") + ylab("")+
  labs(title = "")+
  theme_bw()+
  theme(axis.text=element_text(face = "bold", color="black",size = 12),
        axis.text.x = element_text(size = 12),  # X???̶?
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.copy2pdf(file="pbmc.mac.metascape.up.kegg.go.top20.pdf",width=9,height=0.3*(nrow(tmp)+1))





###############################################################################################
#pyscenic visualize
library(SCopeLoomR)  
library(Seurat)       
library(AUCell)       
library(dplyr)        
library(SCENIC)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2) 
library(stringr)
library(circlize)
library(tidyr)
loom <- open_loom("./aucell.loom")
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

table(pbmc$celllineage)

#regulons target genes
regulons_targetgene.list=matrix(,ncol = 2,nrow = nrow(regulons_incidMat))
regulons_targetgene.list[,1]=row.names(regulons_incidMat)
for (i in 1:nrow(regulons_incidMat)) {
  regulons_targetgene.list[i,2]=paste(colnames(regulons_incidMat)[regulons_incidMat[i,]!=0],collapse = ",")
}
write.csv(regulons_targetgene.list,file="pbmc.regulons_targetgene.list.csv")



######################## Visualize
sub_regulonAUC <- regulonAUC[,match(colnames(pbmc),colnames(regulonAUC))]
dim(sub_regulonAUC)
pbmc
identical(colnames(sub_regulonAUC), colnames(pbmc))
cellClusters <- data.frame(row.names = colnames(pbmc), 
                           seurat_clusters = as.character(pbmc$celllineage))
cellTypes <- data.frame(row.names = colnames(pbmc), 
                        celltype = pbmc$celllineage)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 


auc_mat <- assay(sub_regulonAUC)
regulon_assay <- CreateAssayObject(counts = auc_mat)
pbmc[["RegulonAUC"]] <- regulon_assay


DefaultAssay(pbmc) <- "RegulonAUC"


#differentially regulated TFs
results_list <- list()
meta_data <- pbmc@meta.data
all_celltypes <- unique(meta_data$celllineage)
regulons.to.plot=rownames(auc_mat)


for (ct in all_celltypes) {
  cells_ct <- rownames(meta_data[meta_data$celllineage == ct, ])
  
  df <- FetchData(pbmc, vars = c(regulons.to.plot, "irAE")) %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::filter(cell %in% cells_ct, irAE %in% c("non", "severe_irAE"))
  
  df$group <- factor(df$irAE, levels = c("non", "severe_irAE"))
  
  for (r in regulons.to.plot) {
    if (length(unique(df$group)) == 2) {

      avg_by_group <- df %>%
        group_by(group) %>%
        summarise(mean_auc = mean(.data[[r]], na.rm = TRUE)) %>%
        pivot_wider(names_from = group, values_from = mean_auc, names_prefix = "mean_")
      
      delta <- avg_by_group$mean_severe_irAE - avg_by_group$mean_non
      direction <- ifelse(delta > 0, "up", "down")
      
      # Wilcoxon
      test <- wilcox.test(df[[r]] ~ df$group)
      
      results_list[[paste(ct, r, sep = "_")]] <- data.frame(
        celltype = ct,
        regulon = r,
        p.value = test$p.value,
        mean_non = avg_by_group$mean_non,
        mean_severe_irAE = avg_by_group$mean_severe_irAE,
        delta = delta,
        direction = direction
      )
    }
  }
}

# p-adjust
results_ct<- bind_rows(results_list) %>%
  group_by(celltype) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()


#Vlnplot for ETS2 activity in allcells
pdf(file="pbmc.ETS2.boxplot.irAEgroup.allcells.regulon.pdf",width=3,height=6)
VlnPlot(pbmc,group.by="irAE",c("ETS2(+)"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  ylab("Regulon Activity")+
  NoLegend()
dev.off()

#Vlnplot for ETS2 activity in Macrophage
pdf(file="pbmc.ETS2.boxplot.irAEgroup.mac.regulon.pdf",width=3,height=6)
VlnPlot(pbmc,idents="Macrophage",group.by="irAE",c("ETS2(+)"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  ylab("Regulon Activity")+
  NoLegend()
dev.off()


#wilcox.test for ETS2 expression and activity
ETS2_expression <- GetAssayData(pbmc, assay = "RNA",slot="data")[c("ETS2"), ,drop=F]
pbmc=AddMetaData(pbmc,t(ETS2_expression))

#Macrophage
Ets2_non_mac <- pbmc@meta.data[which(pbmc@meta.data$irAE == "non"&pbmc@meta.data$celllineage =="Macrophage"),]$ETS2
Ets2_severe_irAE_mac <- pbmc@meta.data[which(pbmc@meta.data$irAE == "severe_irAE"&pbmc@meta.data$celllineage =="Macrophage"),]$ETS2
wilcox.test(Ets2_severe_irAE_mac, Ets2_non_mac)
mean(Ets2_non_mac)
mean(Ets2_severe_irAE_mac)
wilcox.test(Ets2_severe_irAE_mac, Ets2_non_mac)$p.value #0.02997921


#####ETS2 activaty
ETS2_activaty <- GetAssayData(pbmc, assay = "RegulonAUC")[c("ETS2(+)"), ,drop=F]
pbmc=AddMetaData(pbmc,t(ETS2_activaty))

Ets2_activaty_non_mac <- pbmc@meta.data[which(pbmc@meta.data$irAE == "non"&pbmc@meta.data$celllineage =="Macrophage"),]$`ETS2(+)`
Ets2_activaty_severe_irAE_mac <- pbmc@meta.data[which(pbmc@meta.data$irAE == "severe_irAE"&pbmc@meta.data$celllineage =="Macrophage"),]$`ETS2(+)`
wilcox.test(Ets2_activaty_severe_irAE_mac, Ets2_activaty_non_mac)
mean(Ets2_activaty_non_mac)
mean(Ets2_activaty_severe_irAE_mac)
wilcox.test(Ets2_activaty_severe_irAE_mac, Ets2_activaty_non_mac)$p.value



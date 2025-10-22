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
library(Cairo)
library(SCP)
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




###########################################################################Code for Human BAF data analysis (EGAS0000100676)

#load meta data
pneumonitis_metadata=read.table("../3086-ICI-pneumonitis_metadata.csv",sep = ",",header = T,quote = "")
head(pneumonitis_metadata)
pneumonitis_metadata[,1]=gsub("\"","",pneumonitis_metadata[,1])
pneumonitis_metadata[,2]=gsub("\"","",pneumonitis_metadata[,2])
table(pneumonitis_metadata[,1])
length(table(pneumonitis_metadata[,1]))
row.names(pneumonitis_metadata)=gsub("\"","",row.names(pneumonitis_metadata))
colnames(pneumonitis_metadata)=c("sample","group")

#load counts
BALF=Read10X("../CIP_counts")
BALF = CreateSeuratObject(counts = BALF)

BALF <- NormalizeData(BALF,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)

BALF[["percent.mt"]] <- PercentageFeatureSet(BALF, pattern = "^MT-")
BALF <- ScaleData(BALF)
BALF <- FindVariableFeatures(BALF,
                             selection.method = "vst",
                             nfeatures = 2000)
BALF <- RunPCA(BALF, verbose = FALSE)

BALF@meta.data=cbind(BALF@meta.data, pneumonitis_metadata)


#quality plot
pdf(file="./lung.qc.vlnplot.pdf",width=16,height=6)
VlnPlot(BALF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()


#harmony integration
BALF <- JoinLayers(BALF)
BALF[["RNA"]] <- split(BALF[["RNA"]], f = BALF$sample)
set.seed(100)
BALF <- IntegrateLayers(
  object = BALF, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

BALF <- FindNeighbors(BALF,  reduction = "harmony",dims = 1:30)
BALF <- RunUMAP(BALF, dims = 1:30, reduction = "harmony",return.model = T, verbose = F)
BALF <- FindClusters(BALF, resolution = 0.25)


#Feature plots
png(file="lung.feature.plot.harmony.png",width=1500,height=900)
FeaturePlot(BALF, features = c("PTPRC","CD3E","CD4","CD8A","TRDC","NKG7","CD79A","JCHAIN","C1QA","ADGRE1",
                               "CD86","CD163","S100A8","FLT3","KIT","GATA3"),raster=FALSE,min.cutoff =0,ncol=5,reduction = "umap",cols= c("gray90", "red"))
dev.off()

#Dimplot and Dimplot split by sample or group
png(file="lung.dimplot.res025.harmony.png",width=900,height=900)
DimPlot(BALF,raster=FALSE,label = T,reduction = "umap",cols= mypal_4)
dev.off()
png(file="lung.dimplot.res025.samplesplit.harmony.png",width=2100,height=1600)
DimPlot(BALF,raster=FALSE,label = T,reduction = "umap",cols= mypal_4,ncol=5,split.by = "sample")
dev.off()
png(file="lung.dimplot.res025.groupsplit.harmony.png",width=600,height=1200)
DimPlot(BALF,raster=FALSE,label = T,reduction = "umap",cols= mypal_4,ncol=1,split.by = "group")
dev.off()

#Marker dotplot
png(file="lung.dotplot.dim30res025.harmony.png",width=1200,height=1500)
DotPlot(BALF, features = c("CD3D","CD3E","CD8A","CD79A","CD79B","MS4A1","LYZ","CD68","FCGR3A",
                           "CLEC4C","IL3RA","LILRA4","FCGR3B","PI3","G0S2","MS4A2","TPSAB1",
                           "TPSB2","CD24","KRT13","KRT18"))
dev.off()





BALF <- JoinLayers(BALF)
markers_BALF <- FindAllMarkers(BALF, only.pos = TRUE)
write.csv(markers_BALF,file="markers_BALF.dim30.res0.25.csv")

#remove low quality clusters,then re-run the standard process
BALF = subset(BALF,idents=c(2,4,8),invert=T)
BALF <- JoinLayers(BALF)
BALF <- NormalizeData(BALF,
                      normalization.method = "LogNormalize",
                      scale.factor = 10000)
BALF <- FindVariableFeatures(BALF, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BALF)
BALF <- ScaleData(BALF, features = all.genes)
BALF <- RunPCA(BALF, features = VariableFeatures(object = BALF))


BALF[["RNA"]] <- split(BALF[["RNA"]], f = BALF$sample)
set.seed(100)
BALF <- IntegrateLayers(
  object = BALF, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

BALF <- FindNeighbors(BALF,  reduction = "harmony",dims = 1:30)
BALF <- RunUMAP(BALF, dims = 1:30, reduction = "harmony",return.model = T, verbose = F)
BALF <- FindClusters(BALF, resolution = 0.25)

png(file="lung.dimplot.res01.harmony.afterfilter.png",width=900,height=900)
DimPlot(BALF,raster=FALSE,label = T,reduction = "umap",cols= mypal_4)
dev.off()
png(file="lung.dimplot.res01.samplesplit.harmony.afterfilter.png",width=2100,height=1600)
DimPlot(BALF,raster=FALSE,label = T,reduction = "umap",cols= mypal_4,ncol=5,split.by = "sample")
dev.off()
png(file="lung.dimplot.res01.groupsplit.harmony.afterfilter.png",width=600,height=1200)
DimPlot(BALF,raster=FALSE,label = T,reduction = "umap",cols= mypal_4,ncol=1,split.by = "group")
dev.off()


#Feature and Dotplot and Dimplot

png(file="lung.feature.plot.harmony.afterfilter.png",width=1500,height=900)
FeaturePlot(BALF, features = c("PTPRC","CD3E","CD4","CD8A","TRDC","NKG7","CD79A","JCHAIN","C1QA","ADGRE1",
                               "CD86","CD163","S100A8","FLT3","KIT","GATA3"),raster=FALSE,min.cutoff =0,ncol=5,reduction = "umap",cols= c("gray90", "red"))
dev.off()
png(file="lung.dimplot.res025.harmony.afterfilter.png",width=900,height=900)
DimPlot(BALF,raster=FALSE,label = T,reduction = "umap",cols= mypal_4)
dev.off()
png(file="lung.dimplot.res025.samplesplit.harmony.afterfilter.png",width=2100,height=1600)
DimPlot(BALF,raster=FALSE,label = T,reduction = "umap",cols= mypal_4,ncol=5,split.by = "sample")
dev.off()
png(file="lung.dimplot.res025.groupsplit.harmony.afterfilter.png",width=600,height=1200)
DimPlot(BALF,raster=FALSE,label = T,reduction = "umap",cols= mypal_4,ncol=1,split.by = "group")
dev.off()

png(file="lung.dotplot.dim30res025.harmony.afterfilter.png",width=1200,height=1500)
DotPlot(BALF, features = c("CD3D","CD3E","CD8A","CD79A","CD79B","MS4A1","LYZ","CD68","FCGR3A",
                           "CLEC4C","IL3RA","LILRA4","CLEC10A","CD1C","CLEC9A","FCGR3B","PI3","G0S2","MS4A2","TPSAB1",
                           "TPSB2","CD24","KRT13","KRT18"))
dev.off()

#re-run markers for BALF
BALF <- JoinLayers(BALF)
BALF.markers <- FindAllMarkers(BALF, only.pos = TRUE)
write.csv(BALF.markers,file = './BALF.markers.dim30.res0.25.harmony.afterfilter.csv')


pdf(file="./lung.qc.vlnplot.afterfilter.pdf",width=16,height=6)
VlnPlot(BALF, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()

library(SeuratDisk)
SaveLoom(
  object = BALF,
  filename = "BALF.harmony.dim30.res0.25.loom",
  overwrite = TRUE
)




###########################################################################
# extract Myeloid for down-stream analysis
###########################################################################
#extract and re-run the standard process
BALF_Myeloid=subset(BALF,idents=c(0,2,11))
BALF_Myeloid <- NormalizeData(BALF_Myeloid)
BALF_Myeloid <- FindVariableFeatures(BALF_Myeloid, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BALF_Myeloid)
BALF_Myeloid <- ScaleData(BALF_Myeloid, features = all.genes)
BALF_Myeloid <- RunPCA(BALF_Myeloid, features = VariableFeatures(object = BALF_Myeloid))

BALF_Myeloid[["RNA"]] <- split(BALF_Myeloid[["RNA"]], f = BALF_Myeloid$sample)
set.seed(100)
BALF_Myeloid <- IntegrateLayers(
  object = BALF_Myeloid, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

BALF_Myeloid <- FindNeighbors(BALF_Myeloid,  reduction = "harmony",dims = 1:30)
BALF_Myeloid <- RunUMAP(BALF_Myeloid, dims = 1:30, return.model = T,reduction = "harmony", verbose = F)
BALF_Myeloid <- FindClusters(BALF_Myeloid, resolution = 0.3)

png(file="lung.mac.dimplot.dim30.res03.png",width=900,height=900)
DimPlot(BALF_Myeloid,raster=FALSE,label = T,reduction = "umap",cols= mypal_4)
dev.off()

png(file="lung.mac.dimplot.dim30.res03.samplesplit.png",width=2100,height=1600)
DimPlot(BALF_Myeloid,raster=FALSE,label = T,reduction = "umap",cols= mypal_4,ncol=5,split.by = "sample")
dev.off()

png(file="lung.mac.dimplot.dim30.res03.groupsplit.png",width=600,height=1200)
DimPlot(BALF_Myeloid,raster=FALSE,label = T,reduction = "umap",cols= mypal_4,ncol=1,split.by = "group")
dev.off()

#Dot plot and Feature plot
png(file="lung.mac.dotplot.dim30.res03.png",width=1800,height=1200)
DotPlot(BALF_Myeloid, features = c("S100A8", "S100A9", "S100A12", "VCAN", "FCN1", "LILRA5", "IL1B", "CCL3", 
                                   "CCL4", "CCL20", "CXCL3", "CXCL8", "EREG", "VEGFA", "IL6", "SPP1", "LGMN", 
                                   "CHI3L1", "MERTK", "CHIT1", "FABP4", "PPARG", "SCD", "TGM2", "RBP4", "CES1", 
                                   "TFRC", "MARCO", "PCOLCE2","HBB","HBA1","HBA2","CD3E","CD3D"))
dev.off()

png(file="lung.mac.featureplot.dim30.res03.png",width=1500,height=1800)
FeaturePlot(BALF_Myeloid, features = c("S100A8", "S100A9", "S100A12", "VCAN", "FCN1", "LILRA5", "IL1B", "CCL3", 
                                       "CCL4", "CCL20", "CXCL3", "CXCL8", "EREG", "VEGFA", "IL6", "SPP1", "LGMN", 
                                       "CHI3L1", "MERTK", "CHIT1", "FABP4", "PPARG", "SCD", "TGM2", "RBP4", "CES1", 
                                       "TFRC", "MARCO", "PCOLCE2","HBB","HBA1","HBA2","CD3E","CD3D"),raster=FALSE,min.cutoff =0,ncol=5,reduction = "umap",cols= c("gray90", "red"))
dev.off()

#Find Markers
BALF_Myeloid <- JoinLayers(BALF_Myeloid)
BALF_Myeloid.markers <- FindAllMarkers(BALF_Myeloid, only.pos = TRUE)
write.csv(BALF_Myeloid.markers,file = './BALF_Myeloid.markers.dim30.res0.3.harmony.csv')

#Quality plot
pdf(file="./lung.myeloid.qc.vlnplot.pdf",width=16,height=6)
VlnPlot(BALF_Myeloid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()



#filter cluter of T cells and Epi
BALF_Myeloid=subset(BALF_Myeloid,idents=c("5","8"),invert=T)

#Re-run the standard process
BALF_Myeloid <- NormalizeData(BALF_Myeloid)
BALF_Myeloid <- FindVariableFeatures(BALF_Myeloid, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BALF_Myeloid)
BALF_Myeloid <- ScaleData(BALF_Myeloid, features = all.genes)
BALF_Myeloid <- RunPCA(BALF_Myeloid, features = VariableFeatures(object = BALF_Myeloid))

BALF_Myeloid[["RNA"]] <- split(BALF_Myeloid[["RNA"]], f = BALF_Myeloid$sample)
set.seed(100)
BALF_Myeloid <- IntegrateLayers(
  object = BALF_Myeloid, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

BALF_Myeloid <- FindNeighbors(BALF_Myeloid,  reduction = "harmony",dims = 1:30)
BALF_Myeloid <- RunUMAP(BALF_Myeloid, dims = 1:30, return.model = T,reduction = "harmony", verbose = F)
BALF_Myeloid <- FindClusters(BALF_Myeloid, resolution = 0.3)


#Dimplot
png(file="lung.mac.dimplot.dim30.res03.filter.png",width=900,height=900)
DimPlot(BALF_Myeloid,raster=FALSE,label = T,reduction = "umap",cols= mypal_4)
dev.off()

png(file="lung.mac.dimplot.dim30.res03.samplesplit.png",width=2100,height=1600)
DimPlot(BALF_Myeloid,raster=FALSE,label = T,reduction = "umap",cols= mypal_4,ncol=5,split.by = "sample")
dev.off()

png(file="lung.mac.dimplot.dim30.res03.groupsplit.png",width=600,height=1200)
DimPlot(BALF_Myeloid,raster=FALSE,label = T,reduction = "umap",cols= mypal_4,ncol=1,split.by = "group")
dev.off()




#Find amrkers
BALF_Myeloid <- JoinLayers(BALF_Myeloid)
BALF_Myeloid.markers <- FindAllMarkers(BALF_Myeloid, only.pos = TRUE)
write.csv(BALF_Myeloid.markers,file = './BALF_Myeloid.markers.dim30.res0.3.harmony.filter.csv')

#Quality plot
pdf(file="./lung.myeloid.qc.vlnplot.afterfilter..pdf",width=16,height=6)
VlnPlot(BALF_Myeloid, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
dev.off()



BALF_Myeloid <- RenameIdents(
  BALF_Myeloid, '0' = 'Macrophage_GPD1','1' = 'Macrophage_FABP3','2' = 'Macrophage_Alveolar_PLA2G7','3' = 'Monocyte_VCAN',
  '4' = 'Macrophage_Alveolar_CCL4','5' = 'Macrophage_Alveolar_MT1X', '6' = 'Macrophage_Alveolar_MALAT1','7'='Macrophage_Alveolar_UCHL1')

BALF_Myeloid$celltype=BALF_Myeloid@active.ident
BALF_Myeloid$celllineage=BALF_Myeloid$celltype
BALF_Myeloid$celllineage<-gsub("Macrophage_Alveolar_PDLIM1", "Macrophage", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Macrophage_Alveolar_FABP3", "Macrophage", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Macrophage_Alveolar_PLA2G7", "Macrophage", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Monocyte_VCAN", "Monocyte", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Macrophage_Alveolar_CCL4", "Macrophage", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Macrophage_Alveolar_MT1X", "Macrophage", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Macrophage_Alveolar_MALAT1", "Macrophage", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Macrophage_Alveolar_UCHL1", "Macrophage", BALF_Myeloid$celllineage)
table(BALF_Myeloid$celllineage)



####remove cluster 7,which is Epi
BALF_Myeloid=subset(BALF_Myeloid,idents='Macrophage_Alveolar_UCHL1',invert=T)

#re run the standard process
BALF_Myeloid <- NormalizeData(BALF_Myeloid)
BALF_Myeloid <- FindVariableFeatures(BALF_Myeloid, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(BALF_Myeloid)
BALF_Myeloid <- ScaleData(BALF_Myeloid, features = all.genes)
BALF_Myeloid <- RunPCA(BALF_Myeloid, features = VariableFeatures(object = BALF_Myeloid))

BALF_Myeloid[["RNA"]] <- split(BALF_Myeloid[["RNA"]], f = BALF_Myeloid$sample)
set.seed(100)
BALF_Myeloid <- IntegrateLayers(
  object = BALF_Myeloid, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

BALF_Myeloid <- FindNeighbors(BALF_Myeloid,  reduction = "harmony",dims = 1:30)
BALF_Myeloid <- RunUMAP(BALF_Myeloid, dims = 1:30, return.model = T,reduction = "harmony", verbose = F)
BALF_Myeloid <- FindClusters(BALF_Myeloid, resolution = 0.3)

###########Cell type annotation
BALF_Myeloid <- RenameIdents(
  BALF_Myeloid, '0' = 'Macrophage_GPD1','1' = 'Macrophage_FABP3','2' = 'Monocyte_CCL2','3' = 'Monocyte_VCAN',
  '4' = 'Macrophage_CCL4','5' = 'Macrophage_MT1X', '6' = 'Macrophage_MALAT1')

BALF_Myeloid$celltype=BALF_Myeloid@active.ident
BALF_Myeloid$celllineage=BALF_Myeloid$celltype
BALF_Myeloid$celllineage<-gsub("Macrophage_GPD1", "Macrophage", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Macrophage_FABP3", "Macrophage", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Monocyte_CCL2", "Monocyte", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Monocyte_VCAN", "Monocyte", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Macrophage_CCL4", "Macrophage", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Macrophage_MT1X", "Macrophage", BALF_Myeloid$celllineage)
BALF_Myeloid$celllineage<-gsub("Macrophage_MALAT1", "Macrophage", BALF_Myeloid$celllineage)
table(BALF_Myeloid$celllineage)

library(SCP)
pdf(file="lung.myeloid.CellDimPlot.celltype.pdf",width=8,height=7)
CellDimPlot(BALF_Myeloid, group.by = "celltype",palcolor = mypal_3, reduction = "UMAP")
dev.off()

#Marker Dotplot
#Markers
target_genes <- c("C1QA","C1QB","FABP3","FABP4", "PPARG", "SCD", "TGM2", "RBP4", "CES1", 
                  "TFRC", "MARCO", "PCOLCE2","S100A8", "S100A9", "S100A12", "VCAN", "FCN1",
                  "IL1B", "CCL3", "CCL4", "CCL20", "CXCL3", "CXCL8", "EREG", "IL6",
                  "MT1X","MT1G","MT1F","MT1E","MRC1","CD163")

png(file="lung.myeloid.dotplot.dim30.res03.png",width=1800,height=1200)
DotPlot(BALF_Myeloid, features = target_genes)
dev.off()


#Markers Heatmap 
avg_expr <- AverageExpression(
  BALF_Myeloid,
  features = target_genes,
  group.by = "seurat_clusters",
  assays = "RNA",       
  slot = "data"         
)$RNA
scaled_expr <- t(scale(t(avg_expr)))
library(ComplexHeatmap)
library(circlize)
heatmap_colors <- colorRamp2(
  c(min(scaled_expr), 0, max(scaled_expr)), 
  c("#0da9ce", "white", "#e74a32")  
)


cluster_colors <- setNames(
  mypal_4[1:nlevels(Idents(BALF_Myeloid))],
  levels(Idents(BALF_Myeloid))
)
top_anno <- HeatmapAnnotation(
  Cluster = levels(Idents(BALF_Myeloid)),
  col = list(Cluster = cluster_colors),
  show_legend = TRUE,
  annotation_name_side = "left"
)

#Heatmap
scaled_expr=scaled_expr[match(target_genes,row.names(scaled_expr)),]
pdf(file="myeloid.heatmap.pdf",width=8,height=16)
Heatmap(
  scaled_expr,
  name = "Z-score", 
  col = heatmap_colors,
  cluster_rows = FALSE,        
  cluster_columns = FALSE,     
  column_names_rot = 45,       
  row_names_gp = gpar(fontface = "italic"),
  top_annotation = top_anno,
  rect_gp = gpar(col = "grey90", lwd = 0.5)
)
dev.off()


###Save the Myeloid data for Pyscenic

library(SeuratDisk)
SaveLoom(
  object = BALF_Myeloid,
  filename = "./BALF_Myeloid.loom",
  overwrite = TRUE
)


#extract BALF_Macrophage
BALF_Macrophage=subset(BALF_Myeloid,idents="Macrophage")




###########################################################################
# AUCell to evaluate the AMPK activity
###########################################################################

BALF <- RenameIdents(
  BALF, '0' = 'Myeloid','1' = 'T_cell','2' = 'Myeloid','3' = 'T_cell',
  '4' = 'Prolif_Myeloid','5' = 'Epithelial','6' = 'T_cell','7' = 'Epithelial',
  '8' = 'DC','9' = 'NK','10' = 'DC','11' = 'Myeloid','12' = 'DC','13' = 'B_cell','14' = 'Mast','15'= 'Epithelial')
BALF$celllineage=BALF@active.ident

library(SCP)
pdf(file="BALF.CellDimPlot.celllineage.pdf",width=8,height=7)
CellDimPlot(BALF, group.by = "celllineage",palcolor = mypal_3, reduction = "UMAP")
dev.off()



#AMPK all
library(AUCell)
library(GSEABase)
library(msigdbr)
library(KEGGREST)
expr_matrix <- as.matrix(GetAssayData(BALF, assay = "RNA", slot = "data"))
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
pdf(file="AUCell_buildRankings.all.pdf",width=9,height = 6)
cells_rankings <- AUCell_buildRankings(
  expr_matrix, 
  plotStats = T, 
  nCores = 16)
dev.off()

set.seed(100)
cells_AUC <- AUCell_calcAUC(
  geneSets, 
  cells_rankings,
  aucMaxRank = 628
)

ampk_scores <- getAUC(cells_AUC)["hsa04152_AMPK_SIGNALING_PATHWAY", ]
BALF$AMPK_AUC <- ampk_scores


###compare the AMPK in all cells
AMPK_CiP_all <- BALF@meta.data[which(BALF@meta.data$group == "CiP"),]$AMPK_AUC
AMPK_Control_all <- BALF@meta.data[which(BALF@meta.data$group == "Control"),]$AMPK_AUC
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_CiP_all, AMPK_Control_all)
# AMPK_HFD_KO vs. AMPK_HFD_WT
mean(AMPK_CiP_all)
mean(AMPK_Control_all)
wilcox.test(AMPK_CiP_all, AMPK_Control_all)$p.value


###compare the AMPK in macrophage
AMPK_CiP_all.mac <- BALF@meta.data[which(BALF@meta.data$group == "CiP" & (colnames(BALF) %in% colnames(BALF_Macrophage))),]$AMPK_AUC
AMPK_Control_all.mac <- BALF@meta.data[which(BALF@meta.data$group == "Control"& (colnames(BALF) %in% colnames(BALF_Macrophage))),]$AMPK_AUC
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_CiP_all.mac, AMPK_Control_all.mac)
# AMPK_HFD_KO vs. AMPK_HFD_WT
mean(AMPK_CiP_all.mac)
mean(AMPK_Control_all.mac)
wilcox.test(AMPK_CiP_all.mac, AMPK_Control_all.mac)$p.value



###########################################################################
# AMPK activity and ETS2 comparison 
###########################################################################
library(ggpubr)
library(ggrepel)
my_comp=list(c("CiP","Control"))

BALF_meta=BALF@meta.data
BALF$group=factor(BALF$group, levels=c("Control","CiP"))
pdf(file="BALF.aucell.boxplot.group.Macrophage.pdf",width=3,height=6)
VlnPlot(subset(BALF,cells=colnames(BALF)[colnames(BALF) %in% colnames(BALF_Macrophage)]),group.by="group",features=c("AMPK_AUC"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()


ets2_expression <- GetAssayData(BALF_Myeloid, assay = "RNA", slot = "data")["ETS2", ]
# ets2_expression add to meta.data
BALF_Myeloid$ETS2_expression <- ets2_expression
my_comp=list(c("CiP","Control"))
pdf(file="BALF.ETS2.boxplot.group.Macrophage.pdf",width=4,height=6)
VlnPlot(BALF_Myeloid,idents="Macrophage",group.by="group",features=c("ETS2_expression"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

BALF_Macrophage
ETS2_CiP_Macrophage <- BALF_Myeloid@meta.data[which(BALF_Myeloid@meta.data$group == "CiP"&BALF_Myeloid@meta.data$celllineage =="Macrophage"),]$ETS2_expression
ETS2_Control_Macrophage <- BALF_Myeloid@meta.data[which(BALF_Myeloid@meta.data$group == "Control"&BALF_Myeloid@meta.data$celllineage =="Macrophage"),]$ETS2_expression
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(ETS2_CiP_Macrophage, ETS2_Control_Macrophage) #2.63e-06
# AMPK_HFD_KO vs. AMPK_HFD_WT
mean(ETS2_CiP_Macrophage)
mean(ETS2_Control_Macrophage)

BALF_Macrophage
pdf(file="BALF.ETS2.dotplot.group.Macrophage.pdf",width=4,height=6)
DotPlot(BALF_Macrophage,group.by="group",features="ETS2",scale = F,scale.min = 0)
dev.off()

library(ggpubr)
library(ggrepel)
my_comp=list(c("CiP","Control"))
pdf(file="BALF.aucell.boxplot.group.Macrophage.pdf",width=4,height=6)
VlnPlot(BALF_Macrophage,group.by="group",features=c("AMPK_AUC"),cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()




##########################################################################################
#Pyscenic visualize
##########################################################################################
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

#regulons_incidMat
regulons_targetgene.list=matrix(,ncol = 2,nrow = nrow(regulons_incidMat))
regulons_targetgene.list[,1]=row.names(regulons_incidMat)
for (i in 1:nrow(regulons_incidMat)) {
  regulons_targetgene.list[i,2]=paste(colnames(regulons_incidMat)[regulons_incidMat[i,]!=0],collapse = ",")
}
write.csv(regulons_targetgene.list,file="balf.myeloid.regulons_targetgene.list.csv")

########################
sub_regulonAUC <- regulonAUC[,match(colnames(BALF_Myeloid),colnames(regulonAUC))]
dim(sub_regulonAUC)
BALF_Myeloid
identical(colnames(sub_regulonAUC), colnames(BALF_Myeloid))

cellClusters <- data.frame(row.names = colnames(BALF_Myeloid), 
                           seurat_clusters = as.character(BALF_Myeloid$celllineage))
cellTypes <- data.frame(row.names = colnames(BALF_Myeloid), 
                        celltype = BALF_Myeloid$celllineage)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 


# pyscenic aucell.loom to auc matrix
auc_mat <- assay(sub_regulonAUC)
regulon_assay <- CreateAssayObject(counts = auc_mat)
BALF_Myeloid[["RegulonAUC"]] <- regulon_assay
DefaultAssay(BALF_Myeloid) <- "RegulonAUC"


#compare the TF activity between Control and Cip 
results_list <- list()
meta_data <- BALF_Myeloid@meta.data
all_celltypes <- unique(meta_data$celllineage)
regulons.to.plot=rownames(auc_mat)


for (ct in all_celltypes) {
  cells_ct <- rownames(meta_data[meta_data$celllineage == ct, ])
  
  df <- FetchData(BALF_Myeloid, vars = c(regulons.to.plot, "group")) %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::filter(cell %in% cells_ct, group %in% c("Control", "CiP"))
  
  df$group <- factor(df$group, levels = c("Control", "CiP"))
  
  for (r in regulons.to.plot) {
    if (length(unique(df$group)) == 2) {
      avg_by_group <- df %>%
        group_by(group) %>%
        summarise(mean_auc = mean(.data[[r]], na.rm = TRUE)) %>%
        pivot_wider(names_from = group, values_from = mean_auc, names_prefix = "mean_")
      delta <- avg_by_group$mean_CiP - avg_by_group$mean_Control
      direction <- ifelse(delta > 0, "up", "down")
      test <- wilcox.test(df[[r]] ~ df$group)
      results_list[[paste(ct, r, sep = "_")]] <- data.frame(
        celltype = ct,
        regulon = r,
        p.value = test$p.value,
        mean_Control = avg_by_group$mean_Control,
        mean_CiP = avg_by_group$mean_CiP,
        delta = delta,
        direction = direction
      )
    }
  }
}

results_ct.cip.control <- bind_rows(results_list) %>%
  group_by(celltype) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

#save results
print(head(results_ct.cip.control, 10))
write.csv(results_ct.cip.control,file="tf.myeloid.results_ct.cip.control.csv")


library(ggpubr)
library(ggrepel)

#ETS2 activity between groups
my_comp=list(c("CiP","Control"))
pdf(file="BALF.ETS2.TF.boxplot.group.Macrophage.pdf",width=4,height=6)
VlnPlot(BALF_Macrophage,group.by="group",features="ETS2(+)",cols=c("#66B99F", "#EC8254"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

BALF_Myeloid@active.assay="RegulonAUC"
BALF_Myeloid@active.ident=BALF_Myeloid$celllineage
pdf(file="BALF.ETS2.TF.dotplot.group.Macrophage.pdf",width=4,height=6)
DotPlot(BALF_Macrophage,group.by="group",features="ETS2(+)",scale = F,scale.min = 0)
dev.off()

#TF activity values extract and compare
ets2tf_activity <- GetAssayData(BALF_Myeloid, assay = "RegulonAUC")["ETS2(+)", ]
BALF_Myeloid$ETS2TF_activity <- ets2tf_activity

ETS2TF_CiP_Macrophage <- BALF_Myeloid@meta.data[which(BALF_Myeloid@meta.data$group == "CiP"&BALF_Myeloid@meta.data$celllineage =="Macrophage"),]$ETS2TF_activity
ETS2TF_Control_Macrophage <- BALF_Myeloid@meta.data[which(BALF_Myeloid@meta.data$group == "Control"&BALF_Myeloid@meta.data$celllineage =="Macrophage"),]$ETS2TF_activity

wilcox.test(ETS2TF_CiP_Macrophage, ETS2TF_Control_Macrophage)
mean(ETS2TF_CiP_Macrophage)
mean(ETS2TF_Control_Macrophage)
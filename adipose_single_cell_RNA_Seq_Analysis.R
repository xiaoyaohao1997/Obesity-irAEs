library(DropletUtils)
library(Seurat)
library(ggsci)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(BiocParallel)
library(harmony)
library(DoubletFinder)
library(clustree)
library(patchwork)
library(tidyverse)
library(paletteer)
library(hdf5r)
library(presto)
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', "#b20000",'#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
mypal_1 <- union(pal_npg("nrc", alpha = 0.7)(10),my36colors)
#red #B53E2B
#blue #1F78B4
mypal_2<-c('#A6CEE3','#1F78B4','#B2DF8A','#33A02C','#FB9A99','#B53E2B','#FDBF6F','#FF7F00','#CAB2D6','#6A3D9A','#FFFF99','#B15928'
           ,"#BBBDDC","#C1E6F3","#FEE1D2","#7DADC6","#F7A073","#E6754B","#AD7F7B")
mypal_3<-union(mypal_2,mypal_1)
mypal_4<-unlist(mapply(brewer.pal, brewer.pal.info[brewer.pal.info$category == 'qual',]$maxcolors, rownames(brewer.pal.info[brewer.pal.info$category == 'qual',])))[-c(6,17,18,19)][-c(7)]
mypal_5<-c("#466791","#60bf37","#953ada","#4fbe6c","#ce49d3","#a7b43d","#5a51dc","#d49f36","#552095","#507f2d","#db37aa","#84b67c","#a06fda","#df462a","#5b83db","#c76c2d","#4f49a3","#82702d","#dd6bbb","#334c22","#d83979","#55baad","#dc4555","#62aad3","#8c3025","#417d61","#862977","#bba672","#403367","#da8a6d","#a79cd4","#71482c","#c689d0","#6b2940","#d593a7","#895c8b","#bd5975")
mypal_6<-c('#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000')

counts_fat12<- Read10X_h5("./filtered_feature_bc_matrix.fat.h5")

##split sample id index
sample_ids <- sapply(colnames(counts_fat12), function(bc) {
  parts <- strsplit(bc, "-")[[1]]
  as.numeric(parts[length(parts)]) 
})
head(sample_ids)

sample_map <- c(
  "1" = "fat1-1", "2" = "fat1-2", "3" = "fat1-3",
  "4" = "fat2-1", "5" = "fat2-2", "6" = "fat2-3",
  "7" = "fat4-1", "8" = "fat4-2", "9" = "fat4-3",
  "10" = "fat3-1", "11" = "fat3-2", "12" = "fat3-3"
)


#####split count into sample level counts
split_counts_fat <- list()
for (id in 1:12) {
  idx <- which(sample_ids == id)
  if (length(idx) > 0) {
    sub_matrix <- counts_fat12[, idx, drop = FALSE]
    split_counts_fat[[sample_map[as.character(id)]]] <- Seurat::CreateSeuratObject(
      counts = sub_matrix,
      project = sample_map[as.character(id)]
    )
  } else {
    warning(paste("", id, ""))
  }
}


#####QC
split_counts_fat <- lapply(split_counts_fat, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  return(obj)
})

#####QC visualize
for (i in 1:12) {
  pdf(file=paste0(split_counts_fat[[i]]@project.name,".vlnplot.counts.features.mt.pdf"),width=12,height=6)
  p=VlnPlot(split_counts_fat[[i]], features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
  print(p)
  dev.off()
}


split_fat_high <- lapply(split_counts_fat, function(obj) {
  obj<-subset(obj,subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)  
  return(obj)
})

######################################################################################
### doubletfinder
######################################################################################
split_fat_singlet <- lapply(split_fat_high, function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- RunUMAP(obj, dims = 1:10)
  ## 2. pK Identification (no ground-truth)
  sweep.res.list_fat12 <- paramSweep(obj, PCs = 1:10, sct = FALSE)
  sweep.stats_fat12 <- summarizeSweep(sweep.res.list_fat12, GT = FALSE)
  png(paste0(obj@project.name,".paramSweep_plot.png"), width = 800, height = 600)
  bcmvn_fat12 <- find.pK(sweep.stats_fat12)
  dev.off()
  mpK_fat12<-as.numeric(as.vector(bcmvn_fat12$pK[which.max(bcmvn_fat12$BCmetric)]))
  ## 3. Run DoubletFinder
  DoubletRate_fat12 = ncol(obj)*8*1e-6 ## ~1000 cells, and a multiplet rate of ~0.8% according to 10XGenomics
  nExp_poi_fat12 <- round(DoubletRate_fat12*ncol(obj)) 
  obj <- doubletFinder(obj, PCs = 1:10, pN = 0.25, pK = mpK_fat12, nExp = nExp_poi_fat12, reuse.pANN = NULL, sct = F)
  ## 4. Remove doublets
  table(obj@meta.data[,paste0("DF.classifications_0.25_",mpK_fat12,"_",nExp_poi_fat12)])
  obj$doubletfinder<-obj@meta.data[,paste0("DF.classifications_0.25_",mpK_fat12,"_",nExp_poi_fat12)]
  seurat_obj_singlet<-subset(obj,doubletfinder=="Singlet")#8284   
  return(seurat_obj_singlet)
})


######################################################################################
###merge sample level data
######################################################################################
sample_names_fat <- sapply(split_fat_singlet, function(obj) obj@project.name)
merged_seurat_fat <- merge(
  x = split_fat_singlet[[1]], 
  y = split_fat_singlet[-1], 
  add.cell.ids = sample_names_fat,  
  merge.data = TRUE             
)
merged_seurat_fat <- JoinLayers(merged_seurat_fat)
Layers(merged_seurat_fat, assay = "RNA")  

pdf(file="./mt15/fat.after.qc.vlnplot.pdf",width=16,height=6)
VlnPlot(merged_seurat_fat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0,group.by = "orig.ident")
dev.off()

######################################################################################
##standard single cell process (with out integration)
######################################################################################
merged_seurat_fat <- NormalizeData(merged_seurat_fat)
merged_seurat_fat <- FindVariableFeatures(merged_seurat_fat, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(merged_seurat_fat)
merged_seurat_fat <- ScaleData(merged_seurat_fat, features = all.genes)
merged_seurat_fat <- RunPCA(merged_seurat_fat, features = VariableFeatures(object = merged_seurat_fat))
ElbowPlot(merged_seurat_fat)
merged_seurat_fat <- FindNeighbors(merged_seurat_fat, dims = 1:30)
merged_seurat_fat <- FindClusters(merged_seurat_fat, resolution = 0.5)
merged_seurat_fat <- RunUMAP(merged_seurat_fat, dims = 1:30)

######################################################################################
##visualize (with out batch effect removing)
######################################################################################
pdf(file="./dimplot.umap.fat.nomerge.pdf",width = 12,height = 12)
DimPlot(merged_seurat_fat, reduction = "umap",label = T,cols = my36colors)
dev.off()
pdf(file="./dimplot.umap.fat.nomerge.group.by.origsampleid.pdf",width = 12,height = 12)
DimPlot(merged_seurat_int, reduction = "umap",label = T,cols = my36colors,group.by = "orig.ident")
dev.off()
pdf(file="./dimplot.umap.fat.nomerge.split.by.origsampleid.pdf",width = 16,height = 12)
DimPlot(merged_seurat_fat, reduction = "umap",split.by = 'orig.ident',cols = my36colors,ncol = 4)
dev.off()


saveRDS(merged_seurat_fat,file="./merged_seurat_fat.RDS")

###############################################################
#######  Harmony integrateion
##############################################################
merged_seurat_fat=readRDS('./merged_seurat_fat.RDS')

options(future.globals.maxSize = 300 * 1024^3)
merged_seurat_fat[["RNA"]] <- split(merged_seurat_fat[["RNA"]], f = merged_seurat_fat$orig.ident)

set.seed(100)
merged_seurat_fat <- IntegrateLayers(
  object = merged_seurat_fat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

merged_seurat_fat <- FindNeighbors(merged_seurat_fat, reduction = "harmony", dims = 1:30)
merged_seurat_fat <- RunUMAP(merged_seurat_fat, reduction = "harmony", dims = 1:30)
merged_seurat_fat <- FindClusters(merged_seurat_fat, resolution = 0.6)
png(file="./dimplot.umap.cluster.label.fat.harmony.merge.dim30.res0.6.validation1.png",width = 1700,height = 800)
p1=DimPlot(merged_seurat_fat, reduction = "umap", label = T,cols = mypal_3,raster = F)
p2=DimPlot(merged_seurat_fat, group.by = "annotation",reduction = "umap", label = T,cols = mypal_3,raster = F)
print(p1+p2)
dev.off()

png(file="fat.feature.plot.validation.png",width=1500,height=900)
FeaturePlot(merged_seurat_fat, features = c("Ptprc","Cd3e","Cd8a","Cd4",
                                            "Trdc","Nkg7","Cd79a","Jchain",
                                            "C1qa","Adgre1","Cd86","Cd163","S100a8",
                                            "Flt3","Kit","Gata3"),raster=FALSE,min.cutoff =0,ncol=5,reduction = "umap",cols= c("gray90", "red"))
dev.off()
png(file="fat.dotplot.validation.0.6.png",width=1500,height=1500)
DotPlot(merged_seurat_fat, features = c("Ptprc","Cd3d","Cd3e","Cd8a","Cd4",
                                        "Trdc","Nkg7","Cd79a","Jchain",
                                        "C1qa","Adgre1","Cd86","Cd163","S100a8",
                                        "Flt3","Kit","Gata3","Mki67","Top2a","Cdh5","Pecam1","Col1a1","Col1a2","Rgs5","Notch3"))
dev.off()



################removing non-immune cell clusters and clusters with mutually exclusive or multiple marker expressions
merged_seurat_fat <- subset(merged_seurat_fat, idents = c(10,11,15,17,19,20), invert = TRUE)#1:30 0.6 seed 100

############### re-run harmony integration
#mapping origin id to samples id
merged_seurat_fat$sample<-merged_seurat_fat$orig.ident
merged_seurat_fat$sample<-gsub("fat1-1", "ND_WT_1", merged_seurat_fat$sample)
merged_seurat_fat$sample<-gsub("fat1-2", "ND_WT_2", merged_seurat_fat$sample)
merged_seurat_fat$sample<-gsub("fat1-3", "ND_WT_3", merged_seurat_fat$sample)
merged_seurat_fat$sample<-gsub("fat2-1", "ND_KO_1", merged_seurat_fat$sample)
merged_seurat_fat$sample<-gsub("fat2-2", "ND_KO_2", merged_seurat_fat$sample)
merged_seurat_fat$sample<-gsub("fat2-3", "ND_KO_3", merged_seurat_fat$sample)
merged_seurat_fat$sample<-gsub("fat4-1", "HFD_WT_1", merged_seurat_fat$sample)
merged_seurat_fat$sample<-gsub("fat4-2", "HFD_WT_2", merged_seurat_fat$sample)
merged_seurat_fat$sample<-gsub("fat4-3", "HFD_WT_3", merged_seurat_fat$sample)
merged_seurat_fat$sample<-gsub("fat3-1", "HFD_KO_1", merged_seurat_fat$sample)
merged_seurat_fat$sample<-gsub("fat3-2", "HFD_KO_2", merged_seurat_fat$sample)
merged_seurat_fat$sample<-gsub("fat3-3", "HFD_KO_3", merged_seurat_fat$sample)

#### re-run standard scRNA process
merged_seurat_fat <- JoinLayers(merged_seurat_fat)
DefaultAssay(merged_seurat_fat) <- "RNA"
merged_seurat_fat <- NormalizeData(merged_seurat_fat)
merged_seurat_fat <- FindVariableFeatures(merged_seurat_fat, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(merged_seurat_fat)
merged_seurat_fat <- ScaleData(merged_seurat_fat, features = all.genes)
merged_seurat_fat <- RunPCA(merged_seurat_fat, features = VariableFeatures(object = merged_seurat_fat))


merged_seurat_fat[["RNA"]] <- split(merged_seurat_fat[["RNA"]], f = merged_seurat_fat$sample)
set.seed(100)
merged_seurat_fat <- IntegrateLayers(
  object = merged_seurat_fat, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)

merged_seurat_fat <- FindNeighbors(merged_seurat_fat, reduction = "harmony", dims = 1:30)
merged_seurat_fat <- RunUMAP(merged_seurat_fat, reduction = "harmony", dims = 1:30, return.model = T, verbose = F)
merged_seurat_fat <- FindClusters(merged_seurat_fat, resolution = 0.2)


png(file="dimplot.umap.cluster.label.fat.harmony.merge.dim30.res060.then.res0.2.afterfilter.validation.rescale.png",width=900,height=800)
DimPlot(merged_seurat_fat, reduction = "umap",cols = mypal_3,label = TRUE,raster=FALSE)
dev.off()

png(file="fat.feature.plot.afterfilter.rescale.reharmony.png",width=1500,height=900)
FeaturePlot(merged_seurat_fat, features = c("Ptprc","Cd3e","Cd8a","Cd4",
                                            "Trdc","Nkg7","Cd79a","Jchain",
                                            "C1qa","Adgre1","Cd86","Cd163","S100a8",
                                            "Flt3","Kit","Gata3"),raster=FALSE,min.cutoff =0,ncol=5,reduction = "umap",cols= c("gray90", "red"))
dev.off()

merged_seurat_fat <- JoinLayers(merged_seurat_fat)



################################################################################################
############Find marker genes and annotation cell types
################################################################################################

fat.merged.harmony.markers <- FindAllMarkers(merged_seurat_fat, only.pos = TRUE)
write.csv(fat.merged.harmony.markers,file = './fat.merged.harmony.markers.afterfilter.csv')

####annotation cell types
merged_seurat_fat <- RenameIdents(
  merged_seurat_fat, '0' = 'Macrophage_Ednrb','1' = 'B_cell','2' = 'Macrophage_Vcam1','3' = 'Macrophage_Ear2',
  '4' = 'CD8_T','5' = 'CD4_T','6' = 'cDC2','7' = 'NK',
  '8' = 'gd_T','9' = 'cDC1','10' = 'Mast','11' = 'DC_Fscn1',
  '12' = 'ILC2','13' = 'Neutrophil','14' = 'Macrophage_Gngt2','15' = 'Plasma')

###celltypes mapping to cell lineage
merged_seurat_fat <- RenameIdents(
  merged_seurat_fat, '0' = 'Macrophage_Gas6','1' = 'B_cell','2' = 'Macrophage_Vcam1','3' = 'Macrophage_Ear2',
  '4' = 'CD8_T','5' = 'CD4_T','6' = 'cDC2','7' = 'NK',
  '8' = 'cDC1','9' = 'gd_T','10' = 'Mast','11' = 'DC_Fscn1',
  '12' = 'ILC2','13' = 'Neutrophil','14' = 'Plasma')


###mapping cell type to cell lineage
table(Idents(merged_seurat_fat))
merged_seurat_fat$celltype<-Idents(merged_seurat_fat)
merged_seurat_fat$celllineage<-merged_seurat_fat$celltype
merged_seurat_fat$celllineage<-gsub("Macrophage_Gas6", "Macrophage", merged_seurat_fat$celllineage)
merged_seurat_fat$celllineage<-gsub("Macrophage_Vcam1", "Macrophage", merged_seurat_fat$celllineage)
merged_seurat_fat$celllineage<-gsub("Macrophage_Ear2", "Macrophage", merged_seurat_fat$celllineage)
merged_seurat_fat$celllineage<-gsub("cDC1", "DC", merged_seurat_fat$celllineage)
merged_seurat_fat$celllineage<-gsub("cDC2", "DC", merged_seurat_fat$celllineage)
merged_seurat_fat$celllineage<-gsub("DC_Fscn1", "DC", merged_seurat_fat$celllineage)
merged_seurat_fat$celllineage<-gsub("CD8_T", "T_cell", merged_seurat_fat$celllineage)
merged_seurat_fat$celllineage<-gsub("CD4_T", "T_cell", merged_seurat_fat$celllineage)
merged_seurat_fat$celllineage<-gsub("gd_T", "T_cell", merged_seurat_fat$celllineage)
table(merged_seurat_fat$celllineage)

##mapping sample group 
merged_seurat_fat$group <- sub("_\\d+$", "", merged_seurat_fat$sample)
merged_seurat_fat@meta.data$group=factor(merged_seurat_fat@meta.data$group,levels = c("HFD_KO","HFD_WT","ND_KO","ND_WT"))
unique(merged_seurat_fat@meta.data$group)
table(merged_seurat_fat@meta.data$group)
library(dplyr)
library(future)


####################################################################################################
### Visualize 
####################################################################################################
DefaultAssay(merged_seurat_fat) <- "RNA"
library(SCP)
pdf(file="fat.CellDimPlot.celllineage.pdf",width=8,height=7)
CellDimPlot(merged_seurat_fat, group.by = "celllineage",palcolor = mypal_3, reduction = "UMAP")
dev.off()

pdf(file="fat.CellDimPlot.group.pdf",width=8,height=7)
CellDimPlot(merged_seurat_fat, group.by = "group",palcolor = c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA")
            , reduction = "UMAP")
dev.off()

pdf(file="fat.CellDimPlot.sample.pdf",width=8,height=7)
CellDimPlot(merged_seurat_fat, group.by = "sample",palcolor = my36colors, reduction = "UMAP")
dev.off()

pdf(file="fat.CellDimPlot.celltype.pdf",width=8,height=7)
CellDimPlot(merged_seurat_fat, group.by = "celltype",palcolor = mypal_4, reduction = "UMAP")
dev.off()


##order celllineage by cell number
table(merged_seurat_fat$celllineage)[order(table(merged_seurat_fat$celllineage),decreasing=T)]
merged_seurat_fat$celllineage<-factor(merged_seurat_fat$celllineage,levels=names(table(merged_seurat_fat$celllineage)[order(table(merged_seurat_fat$celllineage),decreasing=T)]))


##############################################################################################################################
##############################################################################################################################
# cell proportion plots
##############################################################################################################################
##############################################################################################################################

merged_seurat_fat$object<-rep("adipose",ncol(merged_seurat_fat))

cell.prop.sample<-as.data.frame(prop.table(table(merged_seurat_fat$celllineage,merged_seurat_fat$object)))
colnames(cell.prop.sample)<-c("cluster","object","proportion")
cell.prop.sample$object="adipose"

#celllineage proportion plots allsamples
pdf(file = "fat.celllineage.prop.allsamples.pdf",width=3,height=7)
ggplot(cell.prop.sample,aes(object,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.text.x = element_text(angle = 45))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=mypal_3)
dev.off()

#celllineage proportion plots every single sample
pdf(file = "fat.celllineage.prop.singlesample.pdf",width=10,height=7)
cell.prop.sample<-as.data.frame(prop.table(table(merged_seurat_fat$celllineage,merged_seurat_fat$sample)))
colnames(cell.prop.sample)<-c("cluster","sample","proportion")
ggplot(cell.prop.sample,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.text.x = element_text(angle = 45))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=mypal_3)
dev.off()

#cell type proportion plots every single sample
pdf(file = "fat.cell.prop.sample.celltype.pdf",width=10,height=7)
cell.prop.sample<-as.data.frame(prop.table(table(merged_seurat_fat$celltype,merged_seurat_fat$sample)))
colnames(cell.prop.sample)<-c("cluster","sample","proportion")
ggplot(cell.prop.sample,aes(sample,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.5,'cm'),axis.text.x = element_text(angle = 45))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values=mypal_4)
dev.off()

#celllineage proportion box plots across groups
prop_celllineage <- merged_seurat_fat@meta.data %>%
  group_by(sample, group, celllineage) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()


pdf(file = "fat.cell.prop.sample.celllineage.boxplot.pdf",width=14,height=5)
ggplot(prop_celllineage, aes(x = celllineage, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +  
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             size = 1, alpha = 0.6) +  
  labs(x = "Cell Type", y = "Cell Proportion", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "grey90", size = 0.2)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)))+
  scale_fill_manual(values=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA")) 
dev.off()


#cell type proportion box plots across groups
prop_celltype <- merged_seurat_fat@meta.data %>%
  group_by(sample, group, celltype) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
pdf(file = "fat.cell.prop.sample.celltype.boxplot.pdf",width=14,height=5)

ggplot(prop_celltype, aes(x = celltype, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             size = 1, alpha = 0.6) + 
  labs(x = "Cell Type", y = "Cell Proportion", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "grey90", size = 0.2)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)))+
  scale_fill_manual(values=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA")) 
dev.off()






##############################################################################################################################
# Markers for cell type and cell lineage
##############################################################################################################################
library(future)
plan(multisession, workers=8)
markers_merged_seurat_fat_celllineage <- FindAllMarkers(merged_seurat_fat,group.by="celllineage",only.pos = TRUE)
write.csv(markers_merged_seurat_fat_celllineage,"./markers_merged_seurat_fat_celllineage.csv")
plan(sequential)


plan(multisession, workers=8)
markers_merged_seurat_fat_celltype <- FindAllMarkers(merged_seurat_fat,group.by="celltype",only.pos = TRUE)
write.csv(markers_merged_seurat_fat_celltype,"./markers_merged_seurat_fat_celltype.csv")
plan(sequential)

mk <- markers_merged_seurat_fat_celllineage %>%
  group_by(cluster) %>%  # Group by the cell lineage column
  slice_head(n = 5) %>%       # Take the first 5 rows in each group
  ungroup()       
write.csv(mk,"./markers_merged_seurat_fat_celllineage.top5.csv")

##############################################################################################################################
# cell proportion of each macrophage subtype
##############################################################################################################################
fat.mac=subset(merged_seurat_fat,idents=c("Macrophage_Gas6","Macrophage_Vcam1","Macrophage_Ear2"))
fat.mac@active.ident=fat.mac$celltype
prop_celltype.mac <- fat.mac@meta.data %>%
  group_by(sample, group, celltype) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

pdf(file = "fat.mac.cell.prop.sample.celltype.boxplot.pdf",width=10,height=5)
ggplot(prop_celltype.mac, aes(x = celltype, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             size = 1, alpha = 0.6) +
  labs(x = "Cell Type", y = "Cell Proportion", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "grey90", size = 0.2)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)))+
  scale_fill_manual(values=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA")) 
dev.off()

pdf(file = "fat.mac.cell.prop.sample.celltype.boxplot.pdf",width=8,height=5)
ggplot(prop_celltype.mac[prop_celltype.mac$group != "ND_WT",], aes(x = celltype, y = proportion, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2), 
             size = 1, alpha = 0.6) + 
  labs(x = "Cell Type", y = "Cell Proportion", fill = "Group") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(color = "grey90", size = 0.2)) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1)))+
  scale_fill_manual(values=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA")) 
dev.off()



##############################################################################################################################
# cell lineage marker violin plot
##############################################################################################################################
library(ggSCvis)
# Set seed for reproducibility
set.seed(123)
Idents(merged_seurat_fat)<-merged_seurat_fat$celllineage
merged_seurat_fat_downsampled <- subset(merged_seurat_fat, downsample = 500)

# Violin
pdf(file = "fat.pancelltype.marker.vlnplot.pdf",width=14,height=5)
ggscplot(object=merged_seurat_fat_downsampled,features=mk$gene,
         featuresAnno=mk$cluster,
         mapping=aes(x=gene_name,y=celllineage,
                     fill=mean_exp,exp=value))+
  geom_scViolin()+
  facet_new(facet_col="featureAnno",scales="free_x",
            strip.col=mypal_3,
            x.angle=90,x.label.hjust=1)
dev.off()

# Heatmap plot of celllineage markers
pdf(file = "fat.pancelltype.marker.heatmap.pdf",width=7,height=10)
ggscplot(object=merged_seurat_fat_downsampled,features=mk$gene,
         featuresAnno=mk$cluster,
         mapping=aes(y=gene_name,x=celllineage,
                     fill=mean_exp))+
  geom_scTile()+
  facet_new(facet_row="featureAnno",scales="free_y",
            strip.col=mypal_3,x.angle=90,
            switch="y")+scale_fill_gradient(low = "white",high = "red") 
dev.off()


##############################################################################################################################
##############################################################################################################################
# differentially expressed genes analysis 
##############################################################################################################################
##############################################################################################################################

#################3.1. Identify valid cell lineages for downstream analysis
plan(multisession, workers=8)
# Generate a table of cell counts per group and lineage
cell_counts <- merged_seurat_fat@meta.data %>%
  group_by(celllineage, group) %>%
  summarise(n_cells = n()) %>%
  tidyr::pivot_wider(names_from = group, values_from = n_cells, values_fill = 0)
print(cell_counts)
# Identify valid lineages (where both groups have at least 3 cells)
valid_lineages <- cell_counts %>%
  filter(HFD_KO >= 3 & HFD_WT >= 3 & ND_KO >= 3 & ND_WT >= 3) %>%
  pull(celllineage)
print(valid_lineages)

#################3.2. Identify DEGs between HFD_KO and HFD_WT
deg_results_HFD_KOvsHFD_WT <- list()
for (lineage in valid_lineages) {
  subset_cells <- subset(merged_seurat_fat, celllineage == lineage)
  degs <- FindMarkers(
    subset_cells,
    ident.1 = "HFD_KO",
    ident.2 = "HFD_WT",
    group.by = "group",
    test.use = "wilcox",  
    min.pct = 0.1,
    logfc.threshold = 0.25)
  degs$celllineage <- lineage
  degs$genesymbol<-row.names(degs)
  deg_results_HFD_KOvsHFD_WT[[lineage]] <- degs
}
degs_HFD_KOvsHFD_WT <- bind_rows(deg_results_HFD_KOvsHFD_WT, .id = "celllineage")
degs_HFD_KOvsHFD_WT$comparison<-rep("HFD_KOvsHFD_WT",nrow(degs_HFD_KOvsHFD_WT))
degssig_HFD_KOvsHFD_WT<-degs_HFD_KOvsHFD_WT[which(degs_HFD_KOvsHFD_WT$p_val_adj<0.05&abs(degs_HFD_KOvsHFD_WT$avg_log2FC)>0.25),]
degssig_HFD_KOvsHFD_WT <- degssig_HFD_KOvsHFD_WT %>%
  mutate(direction = case_when(
    avg_log2FC > 0~ "up", 
    avg_log2FC < 0~ "down"))
degssig_HFD_KOvsHFD_WT <- degssig_HFD_KOvsHFD_WT %>%
  mutate(direction = factor(direction, levels = c("up", "down"))) %>%
  arrange(celllineage, direction)
head(degssig_HFD_KOvsHFD_WT)
table(degssig_HFD_KOvsHFD_WT$celllineage)


#################3.3. Identify DEGs between HFD_KO and ND_KO
deg_results_HFD_KOvsND_KO <- list()
for (lineage in valid_lineages) {
  subset_cells <- subset(merged_seurat_fat, celllineage == lineage)
  degs <- FindMarkers(
    subset_cells,
    ident.1 = "HFD_KO",
    ident.2 = "ND_KO",
    group.by = "group",
    test.use = "wilcox",  
    min.pct = 0.1,
    logfc.threshold = 0.25)
  degs$genesymbol<-row.names(degs)
  degs$celllineage <- lineage
  deg_results_HFD_KOvsND_KO[[lineage]] <- degs
}
degs_HFD_KOvsND_KO <- bind_rows(deg_results_HFD_KOvsND_KO, .id = "celllineage")
degs_HFD_KOvsND_KO$comparison<-rep("HFD_KOvsND_KO",nrow(degs_HFD_KOvsND_KO))
degssig_HFD_KOvsND_KO<-degs_HFD_KOvsND_KO[which(degs_HFD_KOvsND_KO$p_val_adj<0.05&abs(degs_HFD_KOvsND_KO$avg_log2FC)>0.25),]
degssig_HFD_KOvsND_KO <- degssig_HFD_KOvsND_KO %>%
  mutate(direction = case_when(
    avg_log2FC > 0~ "up", 
    avg_log2FC < 0~ "down"))
degssig_HFD_KOvsND_KO <- degssig_HFD_KOvsND_KO %>%
  mutate(direction = factor(direction, levels = c("up", "down"))) %>%
  arrange(celllineage, direction)
head(degssig_HFD_KOvsND_KO)
table(degssig_HFD_KOvsND_KO$celllineage)

#################3.4. Identify DEGs between HFD_KO and ND_WT
deg_results_HFD_KOvsND_WT <- list()
for (lineage in valid_lineages) {
  subset_cells <- subset(merged_seurat_fat, celllineage == lineage)
  degs <- FindMarkers(
    subset_cells,
    ident.1 = "HFD_KO",
    ident.2 = "ND_WT",
    group.by = "group",
    test.use = "wilcox",  
    min.pct = 0.1,
    logfc.threshold = 0.25)
  degs$genesymbol<-row.names(degs)
  degs$celllineage <- lineage
  deg_results_HFD_KOvsND_WT[[lineage]] <- degs
}
degs_HFD_KOvsND_WT <- bind_rows(deg_results_HFD_KOvsND_WT, .id = "celllineage")
degs_HFD_KOvsND_WT$comparison<-rep("HFD_KOvsND_WT",nrow(degs_HFD_KOvsND_WT))

degssig_HFD_KOvsND_WT<-degs_HFD_KOvsND_WT[which(degs_HFD_KOvsND_WT$p_val_adj<0.05&abs(degs_HFD_KOvsND_WT$avg_log2FC)>0.25),]
degssig_HFD_KOvsND_WT <- degssig_HFD_KOvsND_WT %>%
  mutate(direction = case_when(
    avg_log2FC > 0~ "up", 
    avg_log2FC < 0~ "down"))
degssig_HFD_KOvsND_WT <- degssig_HFD_KOvsND_WT %>%
  mutate(direction = factor(direction, levels = c("up", "down"))) %>%
  arrange(celllineage, direction)
head(degssig_HFD_KOvsND_WT)
table(degssig_HFD_KOvsND_WT$celllineage)

dir.create("DEGs")
write.csv(degssig_HFD_KOvsHFD_WT,"./DEGs/fat.degssig_HFD_KOvsHFD_WT.logfc0.25.csv")
write.csv(degssig_HFD_KOvsND_KO,"./DEGs/fat.degssig_HFD_KOvsND_KO.logfc0.25.csv")

#################3.5. DEG visualization
# 1. degssig_HFD_KOvsHFD_WT
degssig_HFD_KOvsHFD_WT$celllineage<-factor(degssig_HFD_KOvsHFD_WT$celllineage,levels=names(table(merged_seurat_fat$celllineage)[order(table(merged_seurat_fat$celllineage),decreasing=T)]))
degssig_HFD_KOvsHFD_WT <- degssig_HFD_KOvsHFD_WT %>%
  mutate(log10p_val_adj = -log10(replace(p_val_adj, p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm = TRUE)))
  )
library(tidyverse)
library(ggrepel)
top10 <- degssig_HFD_KOvsHFD_WT %>% 
  group_by(celllineage) %>%
  top_n(10, abs(avg_log2FC)) %>%
  ungroup()

pdf(file="fat.volcano.pancelltype.hfdko.hfdwt.pdf",width=20,height=5)
p<-ggplot() +
  geom_point(data = degssig_HFD_KOvsHFD_WT, aes(x = avg_log2FC, y = log10p_val_adj),size = 0.8, color ='grey') +
  coord_flip() +
  facet_grid(. ~ celllineage,scales ="free") +
  geom_point(data = top10, aes(x = avg_log2FC, y = log10p_val_adj,color = celllineage)) +
  geom_vline(xfatercept = c(-0.25, 0.25), size = 0.5, color ="grey50", lty ='dashed')+
  scale_color_manual(values = mypal_3) +
  xlab(label ="avg_log2FC(HFD_KOvsHFD)") +
  ylab(label ="") +
  theme_bw()+
  theme( legend.position ='none',
         panel.grid = element_blank(),
         axis.text = element_text(size = 10),
         axis.text.x = element_text(angle = 45, vjust = 0.8),
         strip.text.x = element_text(size = 10, face ='bold')
  )
p + geom_text_repel(
  data = top10,
  aes(x = avg_log2FC, y = log10p_val_adj, label = genesymbol, color = celllineage),
  fontface = 'italic',
  seed = 233,
  size = 3,
  min.segment.length = 0,  
  force = 12,              
  force_pull = 2,         
  box.padding = 0.1,       
  max.overlaps = Inf,      
  segment.linetype = 3,    
  segment.alpha = 0.5,     
  direction = "y",         
  hjust = 0                
)

dev.off()

# 2. degssig_HFD_KOvsND_KO
degssig_HFD_KOvsND_KO$celllineage<-factor(degssig_HFD_KOvsND_KO$celllineage,levels=names(table(merged_seurat_fat$celllineage)[order(table(merged_seurat_fat$celllineage),decreasing=T)]))
degssig_HFD_KOvsND_KO <- degssig_HFD_KOvsND_KO %>%
  mutate(log10p_val_adj = -log10(replace(p_val_adj, p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm = TRUE)))
  )

top10 <- degssig_HFD_KOvsND_KO %>% 
  group_by(celllineage) %>%
  top_n(10, abs(avg_log2FC)) %>%
  ungroup()
pdf(file="fat.volcano.pancelltype.hfdko.ndko.pdf",width=20,height=5)
p<-ggplot() +
  geom_point(data = degssig_HFD_KOvsND_KO, aes(x = avg_log2FC, y = log10p_val_adj),size = 0.8, color ='grey') +
  coord_flip() +
  facet_grid(. ~ celllineage,scales ="free") +
  geom_point(data = top10, aes(x = avg_log2FC, y = log10p_val_adj,color = celllineage)) +
  geom_vline(xfatercept = c(-0.25, 0.25), size = 0.5, color ="grey50", lty ='dashed')+
  scale_color_manual(values = mypal_3) +
  xlab(label ="avg_log2FC(HFD_KOvsND_KO)") +
  ylab(label ="") +
  theme_bw()+
  theme( legend.position ='none',
         panel.grid = element_blank(),
         axis.text = element_text(size = 10),
         axis.text.x = element_text(angle = 45, vjust = 0.8),
         strip.text.x = element_text(size = 10, face ='bold')
  )
p + geom_text_repel(
  data = top10,
  aes(x = avg_log2FC, y = log10p_val_adj, label = genesymbol, color = celllineage),
  fontface = 'italic',
  seed = 233,
  size = 3,
  min.segment.length = 0,  
  force = 12,             
  force_pull = 2,          
  box.padding = 0.1,       
  max.overlaps = Inf,      
  segment.linetype = 3,   
  segment.alpha = 0.5,     
  direction = "y",         
  hjust = 0               
)

dev.off()

#################3.6. Identify overlaped DEGs between different comparisons in the same celllineage
library(dplyr)
# Rename columns in each dataframe to include full comparison names
deg1 <- degssig_HFD_KOvsHFD_WT %>%
  rename_with(~ paste0(., "_HFD_KOvsHFD_WT"), -c(genesymbol, direction, celllineage))
deg2 <- degssig_HFD_KOvsND_KO %>%
  rename_with(~ paste0(., "_HFD_KOvsND_KO"), -c(genesymbol, direction, celllineage))
# Merge all three dataframes, keeping only genes with consistent direction
degssig_overlap <- deg1 %>%
  inner_join(deg2, by = c("genesymbol", "direction", "celllineage"))
# Check the column names to verify
colnames(degssig_overlap)
degssig_overlap <- degssig_overlap %>%
  select(celllineage,direction,genesymbol, everything())
degssig_overlap <- degssig_overlap %>%
  arrange(celllineage, direction)  
degssig_overlap <- degssig_overlap %>%
  mutate(direction = factor(direction, levels = c("up", "down"))) %>%
  arrange(celllineage, direction)
table(degssig_overlap$celllineage)

head(degssig_overlap[which(degssig_overlap$celllineage=="Macrophage"),])
write.csv(degssig_overlap,"./DEGs/degssig_overlap.logfc0.25.csv")

degssig_overlap_up<-degssig_overlap[which(degssig_overlap$direction=="up"),]
degssig_overlap_down<-degssig_overlap[which(degssig_overlap$direction=="down"),]


#################3.7. Overlap degs visualization
# 1. Calculate gene counts and total per lineage for ranking
count_data <- degssig_overlap %>%
  group_by(celllineage, direction) %>%
  summarise(count = n_distinct(genesymbol), .groups = "drop") %>%
  mutate(
    celllineage = reorder(celllineage, count, FUN = sum)  # Rank by total genes (descending)
  )

# 2. Plot with flipped axes and sorted lineages
pdf(file="fat.pancelltype.overlap.degs.number.counts.pdf",width=6,height=6)
ggplot(count_data, aes(y = celllineage, x = count, fill = direction)) +
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
    labels = c("Upregulated", "Downregulated")
  ) +
  labs(
    y = "Cell Lineage", 
    x = "Number of Genes", 
    title = "Differentially Expressed Genes by Cell Lineage",
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


#################3.8. Overlap degs visualization
#################3.8.1 degssig_HFD_KOvsHFD_WT
degssig_HFD_KOvsHFD_WT <- degssig_HFD_KOvsHFD_WT %>%
  mutate(log10p_val_adj = -log10(replace(p_val_adj, p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm = TRUE))))
degssig_HFD_KOvsHFD_WT_mac<-degssig_HFD_KOvsHFD_WT[which(degssig_HFD_KOvsHFD_WT$celllineage=="Macrophage"),]
# ?޸ĺ??Ľű?

pdf(file="fat.mac.pancelltype.overlap.degs.hfdko.hfdwt.volcano.pdf",width=5,height=4.5)
ggplot(degssig_HFD_KOvsHFD_WT_mac, 
       aes(avg_log2FC, log10p_val_adj, color = direction)) +
  geom_point(size = 1) + 
  theme_bw() +
  ggtitle("degssig_HFD_KOvsHFD_WT") +
  xlim(c(-4, 4)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dotted") +
  scale_color_manual(values = c("#FF3333", "#0066CC")) +
  geom_point(data = degssig_HFD_KOvsHFD_WT_mac %>% 
               filter(genesymbol %in% c("Il1b","Cxcl1", "Cxcl2", "Tnf", "Nfkb2")),
             size = 3, color = "black") +  
  geom_point(data = degssig_HFD_KOvsHFD_WT_mac %>% 
               filter(genesymbol %in% c("Il1b", "Cxcl1","Cxcl2", "Tnf", "Nfkb2")),
             size = 2) +  
  geom_text_repel(data = degssig_HFD_KOvsHFD_WT_mac %>% 
                    filter(genesymbol %in% c("Il1b", "Cxcl1","Cxcl2", "Tnf", "Nfkb2")),
                  aes(label = genesymbol),
                  size = 4,
                  col = 'black',
                  box.padding = 0.5,  
                  max.overlaps = Inf)  
dev.off()


#################3.8.2 degssig_HFD_KOvsND_KO
degssig_HFD_KOvsND_KO <- degssig_HFD_KOvsND_KO %>%
  mutate(log10p_val_adj = -log10(replace(p_val_adj, p_val_adj == 0, min(p_val_adj[p_val_adj > 0], na.rm = TRUE))))
degssig_HFD_KOvsND_KO_mac<-degssig_HFD_KOvsND_KO[which(degssig_HFD_KOvsND_KO$celllineage=="Macrophage"),]

pdf(file="fat.mac.pancelltype.overlap.degs.hfdko.ndko.volcano.pdf",width=5,height=4.5)
ggplot(degssig_HFD_KOvsND_KO_mac, 
       aes(avg_log2FC, log10p_val_adj, color = direction)) +
  geom_point(size = 1) +  
  theme_bw() +
  ggtitle("degssig_HFD_KOvsND_KO") +
  xlim(c(-4, 4)) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted") +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dotted") +
  scale_color_manual(values = c("#FF3333", "#0066CC")) +
  geom_point(data = degssig_HFD_KOvsND_KO_mac %>% 
               filter(genesymbol %in% c("Il1b", "Cxcl1","Cxcl2", "Tnf", "Nfkb2")),
             size = 3, color = "black") + 
?
  geom_point(data = degssig_HFD_KOvsND_KO_mac %>% 
               filter(genesymbol %in% c("Il1b","Cxcl1","Cxcl2", "Tnf", "Nfkb2")),
             size = 2) +  
  # ???ӱ?ǩ
  geom_text_repel(data = degssig_HFD_KOvsND_KO_mac %>% 
                    filter(genesymbol %in% c("Il1b", "Cxcl1","Cxcl2", "Tnf", "Nfkb2")),
                  aes(label = genesymbol),
                  size = 4,
                  col = 'black',
                  box.padding = 0.5,  
                  max.overlaps = Inf)

dev.off()


library(ggpubr)
#################3.8.3 Violin plot
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
pdf(file="fat.mac.pancelltype.group.stack.vlnplot.pdf",width=3,height=8.5)
StackedVlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),ident="Macrophage", group.by="group",c("Il1b","Cxcl2","Cxcl1","Tnf","Nfkb2"), pt.size=0, cols=c("#FF8080", "#FDCABF","#9B8EB6"))
dev.off()


#################3.9 comparison of cytokines
cytokines_expression <- GetAssayData(merged_seurat_fat, assay = "RNA",slot="data")[c("Il1b","Cxcl2","Cxcl1","Tnf","Nfkb2"), ]
merged_seurat_fat=AddMetaData(merged_seurat_fat,t(cytokines_expression))

Il1b_HFD_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Il1b
Il1b_HFD_WT_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Il1b
Il1b_ND_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Il1b
wilcox.test(Il1b_HFD_KO_mac, Il1b_HFD_WT_mac) #1.793917e-182
wilcox.test(Il1b_HFD_KO_mac, Il1b_ND_KO_mac) #1.276347e-218


Cxcl2_HFD_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Cxcl2
Cxcl2_HFD_WT_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Cxcl2
Cxcl2_ND_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Cxcl2
wilcox.test(Cxcl2_HFD_KO_mac, Cxcl2_HFD_WT_mac) #1.343017e-185
wilcox.test(Cxcl2_HFD_KO_mac, Cxcl2_ND_KO_mac) #< 2.225e-308

####
Tnf_HFD_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Tnf
Tnf_HFD_WT_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Tnf
Tnf_ND_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Tnf
wilcox.test(Tnf_HFD_KO_mac, Tnf_HFD_WT_mac) #1.097085e-188
wilcox.test(Tnf_HFD_KO_mac, Tnf_ND_KO_mac) #1.206267e-81

####
Nfkb2_HFD_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Nfkb2
Nfkb2_HFD_WT_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Nfkb2
Nfkb2_ND_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Nfkb2
wilcox.test(Nfkb2_HFD_KO_mac, Nfkb2_HFD_WT_mac) #6.968859e-185
wilcox.test(Nfkb2_HFD_KO_mac, Nfkb2_ND_KO_mac) #4.399398e-183




##############################################################################################################################
##############################################################################################################################
# 4. AUCELL to evaluate KEGG pathway
##############################################################################################################################
##############################################################################################################################
library(AUCell)
library(GSEABase)
library(msigdbr)
library(KEGGREST)

#################5.1 AMPK signaling pathway
expr_matrix <- as.matrix(GetAssayData(merged_seurat_fat, assay = "RNA", slot = "data"))
pathway_id <- "mmu04152"
pathway_info <- keggGet(pathway_id)
genes <- pathway_info[[1]]$GENE
gene_symbols <- genes[seq(2, length(genes), by = 2)]  
gene_symbols <- gsub(";.*", "", gene_symbols)  
gene_symbols <- trimws(gene_symbols)           
print(gene_symbols)
ampk_pathway <- list("mmu04152_AMPK_SIGNALING_PATHWAY" = gene_symbols)
ampk_geneset <- GeneSet(
  geneIds = ampk_pathway[[1]], 
  geneIdType = SymbolIdentifier(), 
  setName = names(ampk_pathway)[1]
)
geneSets <- GeneSetCollection(ampk_geneset)
set.seed(100)
pdf(file="AUCell_buildRankings.pdf",width=8,height=6)
cells_rankings <- AUCell_buildRankings(
  expr_matrix, 
  plotStats = T)
dev.off()

# Quantiles for the number of genes detected by cell: 
#   (Non-detected genes are shuffled at the end of the ranking. Keep it in mind when choosing the threshold for calculating the AUC).
# min   1%   5%  10%  50% 100% 
# 212  414  677  950 1900 5997 

cells_AUC <- AUCell_calcAUC(
  geneSets, 
  cells_rankings,
  aucMaxRank = 950
)
ampk_scores <- getAUC(cells_AUC)["mmu04152_AMPK_SIGNALING_PATHWAY", ]
merged_seurat_fat$AMPK_AUC <- ampk_scores

#compare the AUCell score of AMPK in All the cell types
AMPK_HFD_KO <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"),]$AMPK_AUC
AMPK_HFD_WT <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"),]$AMPK_AUC
AMPK_ND_KO <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"),]$AMPK_AUC
AMPK_ND_WT <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"),]$AMPK_AUC
wilcox.test(AMPK_HFD_KO, AMPK_HFD_WT)
wilcox.test(AMPK_HFD_KO, AMPK_ND_KO) 



#compare the AUCell score of AMPK inMacrophage
AMPK_HFD_KO_Macrophage <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$AMPK_AUC
AMPK_HFD_WT_Macrophage <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$AMPK_AUC
AMPK_ND_KO_Macrophage <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$AMPK_AUC
AMPK_ND_WT_Macrophage <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$AMPK_AUC
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_Macrophage, AMPK_HFD_WT_Macrophage)
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_Macrophage, AMPK_ND_KO_Macrophage) 

library(ggpubr)
target.group=c("HFD_KO","HFD_WT","ND_KO")
my_comp=list(c("HFD_KO","ND_KO"),c("HFD_KO","HFD_WT"))



pdf(file="fat.aucell.boxplot.group.inone.pdf",width=4,height=6)
VlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),group.by="group",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()



pdf(file="fat.aucell.boxplot.group.macrophage.pdf",width=5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="Macrophage",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()


pdf(file="fat.aucell.boxplot.group.macrophage.pdf",width=4,height=6)
VlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),group.by="group",ident="Macrophage",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()


fat.mac=subset(merged_seurat_fat,idents=c("Macrophage"))
fat.mac@active.ident=fat.mac$celltype
prop_celltype.mac <- fat.mac@meta.data %>%
  group_by(sample, group, celltype) %>%
  summarise(n = n()) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

pdf(file="fat.Cxcl16.expr.vlnplot.pdf",width=12,height=4)
VlnPlot(fat.mac,features="Cxcl16",split.by="group",cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)
dev.off()

pdf(file="fat.aucell.boxplot.group.macrophage_Vcam1.pdf",width=4,height=6)
VlnPlot(subset(fat.mac,subset = group %in% target.group),group.by="group",ident="Macrophage_Vcam1",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()




pdf(file="fat.aucell.boxplot.group.macrophage_Gas6.pdf",width=5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="Macrophage_Gas6",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()

pdf(file="fat.aucell.boxplot.group.macrophage_Gas6.pdf",width=4,height=6)
VlnPlot(subset(fat.mac,subset = group %in% target.group),group.by="group",ident="Macrophage_Gas6",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

###################fat mac subcluster
fat.mac@meta.data$combined_group <- paste(
  fat.mac@meta.data$celltype, 
  fat.mac@meta.data$group, 
  sep = "_"
)
pdf(file="fat.vcam1.mac.features.dotplot.pdf",width=8,height=6)
DotPlot(fat.mac,features=c("Cxcl16","Cd93","Ms4a7","Cd300ld"),group.by="combined_group")
dev.off()


pdf(file="fat.vcam1.mac.features.group.vlnplot.cxcl16.pdf",width=4,height=6)
VlnPlot(subset(fat.mac,subset = group %in% target.group),group.by="group",ident="Macrophage_Vcam1",c("Cxcl16"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

pdf(file="fat.vcam1.mac.features.group.vlnplot.Ms4a7.pdf",width=4,height=6)
VlnPlot(subset(fat.mac,subset = group %in% target.group),group.by="group",ident="Macrophage_Vcam1",c("Ms4a7"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

pdf(file="fat.vcam1.mac.features.group.vlnplot.cd93.pdf",width=4,height=6)
VlnPlot(subset(fat.mac,subset = group %in% target.group),group.by="group",ident="Macrophage_Vcam1",c("Cd93"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

pdf(file="fat.vcam1.mac.features.group.vlnplot.cd300ld.pdf",width=4,height=6)
VlnPlot(subset(fat.mac,subset = group %in% target.group),group.by="group",ident="Macrophage_Vcam1",c("Cd300ld"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

pdf(file="fat.aucell.boxplot.group.macrophage_Ear2.pdf",width=4,height=6)
VlnPlot(subset(fat.mac,subset = group %in% target.group),group.by="group",ident="Macrophage_Ear2",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

pdf(file="fat.aucell.boxplot.group.Macrophage_Ace.pdf",width=4,height=6)
VlnPlot(subset(fat.mac,subset = group %in% target.group),group.by="group",ident="Macrophage_Ace",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()


pdf(file="fat.aucell.boxplot.group.DC.pdf",width=5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="DC",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()

pdf(file="fat.aucell.boxplot.group.DC.pdf",width=4,height=6)
VlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),ident="DC",group.by="group",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

pdf(file="fat.aucell.boxplot.group.Tcells.pdf",width=5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="T_cell",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()

pdf(file="fat.aucell.boxplot.group.Tcells.pdf",width=4,height=6)
VlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),ident="T_cell",group.by="group",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

pdf(file="fat.aucell.boxplot.group.Bcells.pdf",width=5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="B_cell",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()

pdf(file="fat.aucell.boxplot.group.Bcells.pdf",width=4,height=6)
VlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),ident="B_cell",group.by="group",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()


pdf(file="fat.aucell.boxplot.group.NK.pdf",width=5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="NK",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()


pdf(file="fat.aucell.boxplot.group.Mast.pdf",width=5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="Mast",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()

pdf(file="fat.aucell.boxplot.group.ILC2.pdf",width=5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="ILC2",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()

pdf(file="fat.aucell.boxplot.group.Neutrophil.pdf",width=5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="Neutrophil",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()

pdf(file="fat.aucell.boxplot.group.Plasma.pdf",width=5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="Plasma",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()





#################Macrophage
AMPK_HFD_KO_Macrophage_Vcam1 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celltype =="Macrophage_Vcam1"),]$AMPK_AUC
AMPK_HFD_WT_Macrophage_Vcam1 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celltype =="Macrophage_Vcam1"),]$AMPK_AUC
AMPK_ND_KO_Macrophage_Vcam1 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celltype =="Macrophage_Vcam1"),]$AMPK_AUC
AMPK_ND_WT_Macrophage_Vcam1 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celltype =="Macrophage_Vcam1"),]$AMPK_AUC
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_Macrophage_Vcam1, AMPK_HFD_WT_Macrophage_Vcam1)
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_Macrophage_Vcam1, AMPK_ND_KO_Macrophage_Vcam1) 
# AMPK_HFD_KO vs. AMPK_ND_WT
wilcox.test(AMPK_HFD_KO_Macrophage_Vcam1, AMPK_ND_WT_Macrophage_Vcam1)  



tmp=fat.mac@meta.data[,c("celltype","AMPK_AUC","group")]
df_ampk_mac_subtype=tmp %>%
  group_by(celltype,group) %>%
  summarise(mean_auc = mean(AMPK_AUC))
df_ampk_mac_subtype
# celltype         group  mean_auc
# <fct>            <fct>     <dbl>
#   1 Macrophage_Gas6  HFD_KO   0.0331
# 2 Macrophage_Gas6  HFD_WT   0.0336
# 3 Macrophage_Gas6  ND_KO    0.0340
# 4 Macrophage_Gas6  ND_WT    0.0359
# 5 Macrophage_Vcam1 HFD_KO   0.0327
# 6 Macrophage_Vcam1 HFD_WT   0.0319
# 7 Macrophage_Vcam1 ND_KO    0.0317
# 8 Macrophage_Vcam1 ND_WT    0.0338
# 9 Macrophage_Ear2  HFD_KO   0.0266
# 10 Macrophage_Ear2  HFD_WT   0.0270
# 11 Macrophage_Ear2  ND_KO    0.0275
# 12 Macrophage_Ear2  ND_WT    0.0284

#################Macrophage
AMPK_HFD_KO_Macrophage_Ear2 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celltype =="Macrophage_Ear2"),]$AMPK_AUC
AMPK_HFD_WT_Macrophage_Ear2 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celltype =="Macrophage_Ear2"),]$AMPK_AUC
AMPK_ND_KO_Macrophage_Ear2 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celltype =="Macrophage_Ear2"),]$AMPK_AUC
AMPK_ND_WT_Macrophage_Ear2 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celltype =="Macrophage_Ear2"),]$AMPK_AUC
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_Macrophage_Ear2, AMPK_HFD_WT_Macrophage_Ear2)
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_Macrophage_Ear2, AMPK_ND_KO_Macrophage_Ear2) 
# AMPK_HFD_KO vs. AMPK_ND_WT
wilcox.test(AMPK_HFD_KO_Macrophage_Ear2, AMPK_ND_WT_Macrophage_Ear2)  

#################Macrophage
AMPK_HFD_KO_Macrophage_Gas6 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celltype =="Macrophage_Gas6"),]$AMPK_AUC
AMPK_HFD_WT_Macrophage_Gas6 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celltype =="Macrophage_Gas6"),]$AMPK_AUC
AMPK_ND_KO_Macrophage_Gas6 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celltype =="Macrophage_Gas6"),]$AMPK_AUC
AMPK_ND_WT_Macrophage_Gas6 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celltype =="Macrophage_Gas6"),]$AMPK_AUC
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_Macrophage_Gas6, AMPK_HFD_WT_Macrophage_Gas6)
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_Macrophage_Gas6, AMPK_ND_KO_Macrophage_Gas6) 
# AMPK_HFD_KO vs. AMPK_ND_WT
wilcox.test(AMPK_HFD_KO_Macrophage_Gas6, AMPK_ND_WT_Macrophage_Gas6)  


AMPK_Macrophage_Gas6 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$celltype =="Macrophage_Gas6"),]$AMPK_AUC
AMPK_Macrophage_Ear2 <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$celltype =="Macrophage_Ear2"),]$AMPK_AUC
wilcox.test(AMPK_Macrophage_Gas6, AMPK_Macrophage_Ear2) 
mean(AMPK_Macrophage_Gas6)
mean(AMPK_Macrophage_Ear2)

#################DC
AMPK_HFD_KO_DC <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="DC"),]$AMPK_AUC
AMPK_HFD_WT_DC <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="DC"),]$AMPK_AUC
AMPK_ND_KO_DC <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="DC"),]$AMPK_AUC
AMPK_ND_WT_DC <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celllineage =="DC"),]$AMPK_AUC
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_DC, AMPK_HFD_WT_DC)
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_DC, AMPK_ND_KO_DC) 
# AMPK_HFD_KO vs. AMPK_ND_WT
wilcox.test(AMPK_HFD_KO_DC, AMPK_ND_WT_DC) 

#################T_cell
AMPK_HFD_KO_T_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="T_cell"),]$AMPK_AUC
AMPK_HFD_WT_T_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="T_cell"),]$AMPK_AUC
AMPK_ND_KO_T_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="T_cell"),]$AMPK_AUC
AMPK_ND_WT_T_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celllineage =="T_cell"),]$AMPK_AUC
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_T_cell, AMPK_HFD_WT_T_cell)
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_T_cell, AMPK_ND_KO_T_cell) 
# AMPK_HFD_KO vs. AMPK_ND_WT
wilcox.test(AMPK_HFD_KO_T_cell, AMPK_ND_WT_T_cell) 
# > wilcox.test(AMPK_HFD_KO_T_cell, AMPK_HFD_WT_T_cell)
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  AMPK_HFD_KO_T_cell and AMPK_HFD_WT_T_cell
# W = 6079709, p-value = 2.973e-05
# alternative hypothesis: true location shift is not equal to 0
# 
# > # AMPK_HFD_KO vs. AMPK_HFD_WT
#   > wilcox.test(AMPK_HFD_KO_T_cell, AMPK_ND_KO_T_cell) 
# 
# Wilcoxon rank sum test with continuity correction
# 
# data:  AMPK_HFD_KO_T_cell and AMPK_ND_KO_T_cell
# W = 7244793, p-value = 6.99e-05
# alternative hypothesis: true location shift is not equal to 0


#################B_cell
AMPK_HFD_KO_B_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="B_cell"),]$AMPK_AUC
AMPK_HFD_WT_B_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="B_cell"),]$AMPK_AUC
AMPK_ND_KO_B_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="B_cell"),]$AMPK_AUC
AMPK_ND_WT_B_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celllineage =="B_cell"),]$AMPK_AUC
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_B_cell, AMPK_HFD_WT_B_cell)
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_B_cell, AMPK_ND_KO_B_cell) 
# AMPK_HFD_KO vs. AMPK_ND_WT
wilcox.test(AMPK_HFD_KO_B_cell, AMPK_ND_WT_B_cell)  

#################Heatmap visualization
celllineage<-names(table(merged_seurat_fat@meta.data$celllineage))
data_AMPK<-data.frame()
for (i in celllineage){
  AMPK_HFD_KO_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage ==i),]$AMPK_AUC
  AMPK_HFD_WT_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage ==i),]$AMPK_AUC
  AMPK_ND_KO_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage ==i),]$AMPK_AUC
  AMPK_ND_WT_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celllineage ==i),]$AMPK_AUC
  data_AMPK1<-data.frame(AMPK_HFD_KO=mean(AMPK_HFD_KO_cell),
                         AMPK_HFD_WT=mean(AMPK_HFD_WT_cell),
                         AMPK_ND_KO=mean(AMPK_ND_KO_cell),
                         AMPK_ND_WT=mean(AMPK_ND_WT_cell))
  data_AMPK<-rbind(data_AMPK,data_AMPK1)
}
row.names(data_AMPK)<-celllineage

pdf(file="fat.ampk.pancelltype.heatmap.pdf",width=3.2,height=6)
pheatmap::pheatmap(data_AMPK[,1:3],show_rownames=T,show_colnames=T,cluster_col = FALSE,cluster_row = FALSE,scale="row",
                   color=rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)))
dev.off()



#################5.2 KEGG signalling transduction

pathways_id <- c("mmu04010","mmu04012","mmu04014","mmu04015","mmu04310","mmu04330",
                 "mmu04340","mmu04350","mmu04390","mmu04370","mmu04371","mmu04630",
                 "mmu04064","mmu04668","mmu04066","mmu04068","mmu04020","mmu04070",
                 "mmu04072","mmu04071","mmu04024","mmu04022","mmu04151","mmu04152",
                 "mmu04150")

# ??ʼ??һ?????б?��?洢 GeneSet ????
geneSets_list <- list()
# ʹ?? for ѭ??????????ÿ??ͨ·
for (pid in pathways_id) {
  # ??ȡ KEGG ͨ·??Ϣ
  pathway_info <- keggGet(pid)
  # ??ȡ???????ţ?KEGG ???صĸ?ʽ?ǽ????Ļ???ID?ͷ??ţ?
  genes <- pathway_info[[1]]$GENE
  gene_symbols <- genes[seq(2, length(genes), by = 2)]  # ??ȡż??λ?õĻ???????
  # ???????????ţ?ȥ?????ܵ?ע?ͣ?
  gene_symbols <- gsub(";.*", "", gene_symbols)  # ?Ƴ??ֺź???????
  gene_symbols <- trimws(gene_symbols)           # ȥ??ǰ???ո?
  # ??ȡͨ·???Ʋ?????
  pathway_name <- pathway_info[[1]]$NAME
  pathway_name <- gsub(" - Mus musculus \\$house mouse\\$", "", pathway_name)  # ?Ƴ???????Ϣ
  pathway_name <- gsub(" ", "_", pathway_name)  # ?滻?ո?Ϊ?»???
  # ???? GeneSet ????
  gs <- GeneSet(
    geneIds = gene_symbols,
    geneIdType = SymbolIdentifier(),
    setName = paste0(pid, "_", pathway_name)  # ??????ʽ: "mmuXXXXX_Pathway_Name"
  )
  # ???ӵ? geneSets_list
  geneSets_list[[length(geneSets_list) + 1]] <- gs
}
# ???б?ת??Ϊ GeneSetCollection
geneSets <- GeneSetCollection(geneSets_list)

# ??ÿ??ϸ???Ļ?????????????????????֮ǰû??????
if (!exists("cells_rankings")) {
  cells_rankings <- AUCell_buildRankings(
    expr_matrix,
    plotStats = FALSE,
    nCores = 8
  )
}


set.seed(100)
cells_AUC <- AUCell_calcAUC(
  geneSets,
  cells_rankings,
  aucMaxRank = 950
)

# ?????????ӵ? Seurat ????
auc_matrix <- t(getAUC(cells_AUC))  # ת?þ???ʹϸ??Ϊ?У?ͨ·Ϊ??
merged_seurat_fat <- AddMetaData(merged_seurat_fat, auc_matrix)

print(tail(colnames(merged_seurat_fat@meta.data), 5))  # ?鿴????5??ͨ·?ķ?????


# Calculate significance differences
celllineage<-names(table(merged_seurat_fat@meta.data$celllineage))
pathway<-colnames(auc_matrix)

aucell_pathway<-data.frame()
for (i in pathway){
  HFD_KO_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"),][,i]
  HFD_WT_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"),][,i]
  ND_KO_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"),][,i]
  ND_WT_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"),][,i]
  
  p_HFD_KOvsHFD_WT1<-wilcox.test(HFD_KO_cell,HFD_WT_cell)$p.value
  d_HFD_KOvsHFD_WT1<-mean(HFD_KO_cell)-mean(HFD_WT_cell)
  
  p_HFD_KOvsND_KO1<-wilcox.test(HFD_KO_cell,ND_KO_cell)$p.value
  d_HFD_KOvsND_KO1<-mean(HFD_KO_cell)-mean(ND_KO_cell)
  
  p_HFD_KOvsND_WT1<-wilcox.test(HFD_KO_cell,ND_WT_cell)$p.value
  d_HFD_KOvsND_WT1<-mean(HFD_KO_cell)-mean(ND_WT_cell)
  
  aucell_pathway1<-data.frame(p_HFD_KOvsHFD_WT=p_HFD_KOvsHFD_WT1,
                              d_HFD_KOvsHFD_WT=d_HFD_KOvsHFD_WT1,
                              p_HFD_KOvsND_KO=p_HFD_KOvsND_KO1,
                              d_HFD_KOvsND_KO=d_HFD_KOvsND_KO1,
                              p_HFD_KOvsND_WT=p_HFD_KOvsND_WT1,
                              d_HFD_KOvsND_WT=d_HFD_KOvsND_WT1,
                              mean_HFD_KO=mean(HFD_KO_cell),
                              mean_HFD_WT=mean(HFD_WT_cell),
                              mean_ND_KO=mean(ND_KO_cell),
                              mean_ND_WT=mean(ND_WT_cell))
  aucell_pathway<-rbind(aucell_pathway,aucell_pathway1)
}
row.names(aucell_pathway)<-gsub("_-_Mus_musculus_(house_mouse)","",pathway,fixed = TRUE)

aucell_pathway_up<-aucell_pathway[which(aucell_pathway$p_HFD_KOvsHFD_WT<0.05&
                                          aucell_pathway$d_HFD_KOvsHFD_WT>0&
                                          aucell_pathway$p_HFD_KOvsND_KO<0.05&
                                          aucell_pathway$d_HFD_KOvsND_KO>0#&
                                        #aucell_pathway$p_HFD_KOvsND_WT<0.05&
                                        #aucell_pathway$d_HFD_KOvsND_WT>0
),]
sum(aucell_pathway$p_HFD_KOvsHFD_WT<0.05&aucell_pathway$d_HFD_KOvsHFD_WT>0) #4
sum(aucell_pathway$p_HFD_KOvsND_KO<0.05&aucell_pathway$d_HFD_KOvsND_KO>0) #11
aucell_pathway_down<-aucell_pathway[which(aucell_pathway$p_HFD_KOvsHFD_WT<0.05&
                                            aucell_pathway$d_HFD_KOvsHFD_WT<0&
                                            aucell_pathway$p_HFD_KOvsND_KO<0.05&
                                            aucell_pathway$d_HFD_KOvsND_KO<0#&
                                          #aucell_pathway$p_HFD_KOvsND_WT<0.05&
                                          #aucell_pathway$d_HFD_KOvsND_WT<0
),]
sum(aucell_pathway$p_HFD_KOvsHFD_WT<0.05&aucell_pathway$d_HFD_KOvsHFD_WT<0) #17
sum(aucell_pathway$p_HFD_KOvsND_KO<0.05&aucell_pathway$d_HFD_KOvsND_KO<0) #11

nrow(aucell_pathway_up) #4
nrow(aucell_pathway_down)#9


pdf(file="fat.aucell.pathways.compare.heatmap.4groups.pdf",width = 6,height=6)
pheatmap::pheatmap(aucell_pathway[c(row.names(aucell_pathway_up),row.names(aucell_pathway_down)),c("mean_HFD_KO","mean_HFD_WT","mean_ND_KO","mean_ND_WT")],show_rownames=T,show_colnames=T,cluster_col = FALSE,cluster_row = FALSE,scale="row",
                   color=rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)))
dev.off()
pdf(file="fat.aucell.pathways.compare.heatmap.pdf",width = 5,height=6)
pheatmap::pheatmap(aucell_pathway[c(row.names(aucell_pathway_up),row.names(aucell_pathway_down)),c("mean_HFD_KO","mean_HFD_WT","mean_ND_KO")],show_rownames=T,show_colnames=T,cluster_col = FALSE,cluster_row = FALSE,scale="row",
                   color=rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)))
dev.off()

celllineage<-names(table(merged_seurat_fat@meta.data$celllineage))
data_AMPK<-data.frame()
for (i in celllineage){
  AMPK_HFD_KO_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage ==i),][,"mmu04152_AMPK_signaling_pathway_-_Mus_musculus_(house_mouse)"]
  AMPK_HFD_WT_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage ==i),][,"mmu04152_AMPK_signaling_pathway_-_Mus_musculus_(house_mouse)"]
  AMPK_ND_KO_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage ==i),][,"mmu04152_AMPK_signaling_pathway_-_Mus_musculus_(house_mouse)"]
  AMPK_ND_WT_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celllineage ==i),][,"mmu04152_AMPK_signaling_pathway_-_Mus_musculus_(house_mouse)"]
  data_AMPK1<-data.frame(AMPK_HFD_KO=mean(AMPK_HFD_KO_cell),
                         AMPK_HFD_WT=mean(AMPK_HFD_WT_cell),
                         AMPK_ND_KO=mean(AMPK_ND_KO_cell),
                         AMPK_ND_WT=mean(AMPK_ND_WT_cell))
  data_AMPK<-rbind(data_AMPK,data_AMPK1)
}
row.names(data_AMPK)<-celllineage

pdf(file="fat.aucell.AMPKpathways.compare.heatmap.4groups.pdf",width = 4,height=6)
pheatmap::pheatmap(data_AMPK,show_rownames=T,show_colnames=T,cluster_col = FALSE,cluster_row = FALSE,scale="row",
                   color=rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)))
dev.off()
pdf(file="fat.aucell.AMPKpathways.compare.heatmap.3groups.pdf",width = 3.5,height=6)
pheatmap::pheatmap(data_AMPK[,c("AMPK_HFD_KO","AMPK_HFD_WT","AMPK_ND_KO")],show_rownames=T,show_colnames=T,cluster_col = FALSE,cluster_row = FALSE,scale="row",
                   color=rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)))
dev.off()


###################################################Pathways in Macrophage
aucell_pathway_mac <- data.frame()
for (i in pathway){
  HFD_KO_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),][,i]
  HFD_WT_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),][,i]
  ND_KO_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),][,i]
  ND_WT_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),][,i]
  
  p_HFD_KOvsHFD_WT1<-wilcox.test(HFD_KO_cell,HFD_WT_cell)$p.value
  d_HFD_KOvsHFD_WT1<-mean(HFD_KO_cell)-mean(HFD_WT_cell)
  
  p_HFD_KOvsND_KO1<-wilcox.test(HFD_KO_cell,ND_KO_cell)$p.value
  d_HFD_KOvsND_KO1<-mean(HFD_KO_cell)-mean(ND_KO_cell)
  
  p_HFD_KOvsND_WT1<-wilcox.test(HFD_KO_cell,ND_WT_cell)$p.value
  d_HFD_KOvsND_WT1<-mean(HFD_KO_cell)-mean(ND_WT_cell)
  
  aucell_pathway1<-data.frame(p_HFD_KOvsHFD_WT=p_HFD_KOvsHFD_WT1,
                              d_HFD_KOvsHFD_WT=d_HFD_KOvsHFD_WT1,
                              p_HFD_KOvsND_KO=p_HFD_KOvsND_KO1,
                              d_HFD_KOvsND_KO=d_HFD_KOvsND_KO1,
                              p_HFD_KOvsND_WT=p_HFD_KOvsND_WT1,
                              d_HFD_KOvsND_WT=d_HFD_KOvsND_WT1,
                              mean_HFD_KO=mean(HFD_KO_cell),
                              mean_HFD_WT=mean(HFD_WT_cell),
                              mean_ND_KO=mean(ND_KO_cell),
                              mean_ND_WT=mean(ND_WT_cell))
  aucell_pathway_mac<-rbind(aucell_pathway_mac,aucell_pathway1)
}
row.names(aucell_pathway_mac)<-gsub("_-_Mus_musculus_(house_mouse)","",pathway,fixed = TRUE)

aucell_pathway_mac_up<-aucell_pathway_mac[which(aucell_pathway_mac$p_HFD_KOvsHFD_WT<0.05&
                                                  aucell_pathway_mac$d_HFD_KOvsHFD_WT>0&
                                                  aucell_pathway_mac$p_HFD_KOvsND_KO<0.05&
                                                  aucell_pathway_mac$d_HFD_KOvsND_KO>0#&
                                                #aucell_pathway$p_HFD_KOvsND_WT<0.05&
                                                #aucell_pathway$d_HFD_KOvsND_WT>0
),]
aucell_pathway_mac_down<-aucell_pathway_mac[which(aucell_pathway_mac$p_HFD_KOvsHFD_WT<0.05&
                                                    aucell_pathway_mac$d_HFD_KOvsHFD_WT<0&
                                                    aucell_pathway_mac$p_HFD_KOvsND_KO<0.05&
                                                    aucell_pathway_mac$d_HFD_KOvsND_KO<0#&
                                                  #aucell_pathway$p_HFD_KOvsND_WT<0.05&
                                                  #aucell_pathway$d_HFD_KOvsND_WT<0
),]
aucell_pathway_mac_up
aucell_pathway_mac_down


pdf(file="fat.aucell.pathways.compare.heatmap.4groups.pdf",width = 6,height=6)
pheatmap::pheatmap(aucell_pathway_mac[c(row.names(aucell_pathway_mac_up),row.names(aucell_pathway_mac_down)),c("mean_HFD_KO","mean_HFD_WT","mean_ND_KO","mean_ND_WT")],show_rownames=T,show_colnames=T,cluster_col = FALSE,cluster_row = FALSE,scale="row",
                   color=rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)))
dev.off()
pdf(file="fat.aucell.pathways.compare.heatmap.pdf",width = 5,height=6)
pheatmap::pheatmap(aucell_pathway_mac[c(row.names(aucell_pathway_mac_up),row.names(aucell_pathway_mac_down)),c("mean_HFD_KO","mean_HFD_WT","mean_ND_KO")],show_rownames=T,show_colnames=T,cluster_col = FALSE,cluster_row = FALSE,scale="row",
                   color=rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)))
dev.off()









celllineage<-names(table(merged_seurat_fat@meta.data$celllineage))
data_NFKB<-data.frame()
for (i in celllineage){
  NFKB_HFD_KO_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage ==i),][,"mmu04064_NF-kappa_B_signaling_pathway_-_Mus_musculus_(house_mouse)"]
  NFKB_HFD_WT_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage ==i),][,"mmu04064_NF-kappa_B_signaling_pathway_-_Mus_musculus_(house_mouse)"]
  NFKB_ND_KO_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage ==i),][,"mmu04064_NF-kappa_B_signaling_pathway_-_Mus_musculus_(house_mouse)"]
  NFKB_ND_WT_cell <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celllineage ==i),][,"mmu04064_NF-kappa_B_signaling_pathway_-_Mus_musculus_(house_mouse)"]
  data_NFKB1<-data.frame(NFKB_HFD_KO=mean(NFKB_HFD_KO_cell),
                         NFKB_HFD_WT=mean(NFKB_HFD_WT_cell),
                         NFKB_ND_KO=mean(NFKB_ND_KO_cell),
                         NFKB_ND_WT=mean(NFKB_ND_WT_cell))
  data_NFKB<-rbind(data_NFKB,data_NFKB1)
}
row.names(data_NFKB)<-celllineage

pdf(file="fat.aucell.NFKBpathways.compare.heatmap.4groups.pdf",width = 5,height=6)
pheatmap::pheatmap(data_NFKB,show_rownames=T,show_colnames=T,cluster_col = FALSE,cluster_row = FALSE,scale="row",
                   color=rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)))
dev.off()



pdf(file="fat.aucell.AMPKpathways.vlnplot.pdf",width = 5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",c("mmu04152_AMPK_signaling_pathway_-_Mus_musculus_(house_mouse)"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()

pdf(file="fat.aucell.AMPKpathways.macrophage.vlnplot.pdf",width = 5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="Macrophage",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()

pdf(file="fat.aucell.AMPKpathways.DC.vlnplot.pdf",width = 5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="DC",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()

pdf(file="fat.aucell.AMPKpathways.Tcells.vlnplot.pdf",width = 5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="T_cell",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()


pdf(file="fat.aucell.AMPKpathways.Bcells.vlnplot.pdf",width = 5,height=6)
VlnPlot(merged_seurat_fat,group.by="group",ident="B_cell",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+  
  NoLegend()
dev.off()


#################5.3 DEGs in AMPK pathway
AMPK_genes<-unlist(geneIds(geneSets[[24]]), use.names = FALSE)
degssig_overlap_up_AMPK<-degssig_overlap_up[which(degssig_overlap_up$genesymbol %in% AMPK_genes),]
degssig_overlap_up_AMPK
degssig_overlap_down_AMPK<-degssig_overlap_down[which(degssig_overlap_down$genesymbol %in% AMPK_genes),]
degssig_overlap_down_AMPK
write.csv(degssig_overlap_up_AMPK,file="degssig_overlap_up_AMPK.csv")
write.csv(degssig_overlap_down_AMPK,file="degssig_overlap_down_AMPK.csv")


###################################################################################pyscenic
###################################################################################pyscenic
###################################################################################pyscenic
###################################################################################pyscenic
###################################################################################pyscenic
###################################################################################pyscenic
###################################################################################pyscenic
###################################################################################pyscenic
loom <- open_loom("./pyscenic/aucell.loom")
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

######################## Vis
sub_regulonAUC <- regulonAUC[,match(colnames(merged_seurat_fat),colnames(regulonAUC))]
dim(sub_regulonAUC)
merged_seurat_fat
#ȷ???Ƿ?һ??
identical(colnames(sub_regulonAUC), colnames(merged_seurat_fat))

cellClusters <- data.frame(row.names = colnames(merged_seurat_fat), 
                           seurat_clusters = as.character(merged_seurat_fat$celllineage))
cellTypes <- data.frame(row.names = colnames(merged_seurat_fat), 
                        celltype = merged_seurat_fat$celllineage)
head(cellTypes)
head(cellClusters)
sub_regulonAUC[1:4,1:4] 

auc_mat <- assay(sub_regulonAUC)
regulon_assay <- CreateAssayObject(counts = auc_mat)
merged_seurat_fat[["RegulonAUC"]] <- regulon_assay
DefaultAssay(merged_seurat_fat) <- "RegulonAUC"

regulons_targetgene.list=matrix(,ncol = 2,nrow = nrow(regulons_incidMat))
regulons_targetgene.list[,1]=row.names(regulons_incidMat)
for (i in 1:nrow(regulons_incidMat)) {
  regulons_targetgene.list[i,2]=paste(colnames(regulons_incidMat)[regulons_incidMat[i,]!=0],collapse = ",")
}
write.csv(regulons_targetgene.list,file="fat.regulons_targetgene.list.csv")

##########################################################test 
results_list <- list()
meta_data <- merged_seurat_fat@meta.data
all_celltypes <- unique(meta_data$celllineage)
regulons.to.plot=rownames(auc_mat)


for (ct in all_celltypes) {
  cells_ct <- rownames(meta_data[meta_data$celllineage == ct, ])
  
  df <- FetchData(merged_seurat_fat, vars = c(regulons.to.plot, "group")) %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::filter(cell %in% cells_ct, group %in% c("ND_KO", "HFD_KO"))
  
  df$group <- factor(df$group, levels = c("ND_KO", "HFD_KO"))  # ?̶?˳?򷽱??Ƚ?
  
  for (r in regulons.to.plot) {
    if (length(unique(df$group)) == 2) {
      # ??ֵ????
      avg_by_group <- df %>%
        group_by(group) %>%
        summarise(mean_auc = mean(.data[[r]], na.rm = TRUE)) %>%
        pivot_wider(names_from = group, values_from = mean_auc, names_prefix = "mean_")
      
      # ???췽??
      delta <- avg_by_group$mean_HFD_KO - avg_by_group$mean_ND_KO
      direction <- ifelse(delta > 0, "up", "down")
      
      # Wilcoxon ????
      test <- wilcox.test(df[[r]] ~ df$group)
      
      # ????????
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

# ???Ͻ???????У?? p ֵ
results_ct.hfdko.ndko <- bind_rows(results_list) %>%
  group_by(celltype) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

# ??ӡǰ???п?һ??
print(head(results_ct.hfdko.ndko, 10))
write.csv(results_ct.hfdko.ndko,file="fat.results_tf.hfdko.ndko.csv")

results_ct.hfdko.ndko=read.csv("fat.results_tf.hfdko.ndko.csv")
results_ct.hfdko.ndko$sig="not_sig"
results_ct.hfdko.ndko$sig[results_ct.hfdko.ndko$p.adj<0.05]="sig"
head(results_ct.hfdko.ndko)
results_ct.hfdko.ndko=results_ct.hfdko.ndko[results_ct.hfdko.ndko$p.adj<0.05,]
results_ct.hfdko.ndko$index=paste(results_ct.hfdko.ndko$celltype,results_ct.hfdko.ndko$regulon,results_ct.hfdko.ndko$direction,sep="_")



########### hfd_ko hfd_wt
results_list <- list()
meta_data <- merged_seurat_fat@meta.data
all_celltypes <- unique(meta_data$celllineage)
regulons.to.plot=rownames(auc_mat)


for (ct in all_celltypes) {
  cells_ct <- rownames(meta_data[meta_data$celllineage == ct, ])
  
  df <- FetchData(merged_seurat_fat, vars = c(regulons.to.plot, "group")) %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::filter(cell %in% cells_ct, group %in% c("HFD_WT", "HFD_KO"))
  
  df$group <- factor(df$group, levels = c("HFD_WT", "HFD_KO"))  # ?̶?˳?򷽱??Ƚ?
  
  for (r in regulons.to.plot) {
    if (length(unique(df$group)) == 2) {
      # ??ֵ????
      avg_by_group <- df %>%
        group_by(group) %>%
        summarise(mean_auc = mean(.data[[r]], na.rm = TRUE)) %>%
        pivot_wider(names_from = group, values_from = mean_auc, names_prefix = "mean_")
      
      # ???췽??
      delta <- avg_by_group$mean_HFD_KO - avg_by_group$mean_HFD_WT
      direction <- ifelse(delta > 0, "up", "down")
      
      # Wilcoxon ????
      test <- wilcox.test(df[[r]] ~ df$group)
      
      # ????????
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

# ???Ͻ???????У?? p ֵ
results_ct.hfdko.hfdwt <- bind_rows(results_list) %>%
  group_by(celltype) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()

# ??ӡǰ???п?һ??
print(head(results_ct.hfdko.hfdwt, 10))

write.csv(results_ct.hfdko.hfdwt,file="fat.results_tf.hfdko.hfdwt.csv")


#########common
results_ct.hfdko.hfdwt=read.csv("fat.results_tf.hfdko.hfdwt.csv")
results_ct.hfdko.hfdwt$sig="not_sig"
results_ct.hfdko.hfdwt$sig[results_ct.hfdko.hfdwt$p.adj<0.05]="sig"
results_ct.hfdko.hfdwt=results_ct.hfdko.hfdwt[results_ct.hfdko.hfdwt$p.adj<0.05,]
head(results_ct.hfdko.hfdwt)
results_ct.hfdko.hfdwt$index=paste(results_ct.hfdko.hfdwt$celltype,results_ct.hfdko.hfdwt$regulon,results_ct.hfdko.hfdwt$direction,sep="_")

commonids=intersect(results_ct.hfdko.ndko$index,results_ct.hfdko.hfdwt$index)
results_ct.hfdko.ndko.common=results_ct.hfdko.ndko[match(commonids, results_ct.hfdko.ndko$index),]
results_ct.hfdko.hfdwt.common=results_ct.hfdko.hfdwt[match(commonids, results_ct.hfdko.hfdwt$index),]

results_ct.hfdko.ndko.and.hfdko.hfdwt.common=cbind(results_ct.hfdko.ndko.common,results_ct.hfdko.hfdwt.common)
write.csv(results_ct.hfdko.ndko.and.hfdko.hfdwt.common,file="fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.csv")


#####################################################################################
#####################################################################################
#####################################################################################
fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common=read.csv("../fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.csv")
library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)
regulons_targetgene.list=read.csv("./fat.regulons_targetgene.list.csv")
regulons_targetgene.list=regulons_targetgene.list[,-1]
# head(regulons_targetgene.list,1)
# [,1]       [,2]                                                         
# [1,] "Ascl2(+)" "Lrrc7,Kdm4c,2210417A02Rik,Lmod3,Mfap3l,Nup107,Rad51ap2,Lmf1"
fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$target=regulons_targetgene.list[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$regulon, regulons_targetgene.list[,1]),2]


fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common=na.omit(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common)
fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up=fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="Macrophage" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up",]

fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2=gsub("\\(\\+\\)","",fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$regulon)

tfs=c("Ets2(+)","Nfe2l2(+)","Spi1(+)","Pparg(+)","Irf5(+)")
fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up.tfs=fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up$regulon %in% tfs,]

###dif genes Ets2(+), Nfe2l2(+)
target.list=list()
for (i in 1:nrow(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up.tfs)) {
  target.list[[i]]=strsplit(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common.mac.up.tfs$target[i],split = ",")
}


###tf all, each cell each tf and it's target in a data.frame
target.list.all=list()
for (i in 1:nrow(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common)) {
  target.list.all[[i]]=strsplit(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$target[i],split = ",")
}
tmp=length(unlist(target.list.all))

target.df.all=fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common[rep(1, tmp),]
target.df.all$target.gene=""
target.df.all$target.gene=unlist(target.list.all)

####replace with true data
index=1
for (i in 1:nrow(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common)) {
  target.df.all[index:(index+length(target.list.all[[i]][[1]])-1),1:24]=fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common[i,1:24]
  index=index+length(target.list.all[[i]][[1]])
}

####visualize
head(target.df.all)
target.df.all$source.gene=target.df.all$regulon
target.df.all$source.gene=gsub("\\(\\+\\)","",target.df.all$source.gene)

target.df.all$relation="regulate"

###add regulon and mac
tmp=fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common
tmp$target.gene=fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$regulon
tmp$target.gene=gsub("\\(\\+\\)","",tmp$target.gene)
tmp$source.gene=tmp$celltype
tmp$relation="Direct"
target.df.all=rbind(tmp,target.df.all)
#target.df.all[,]

length(unique(unlist(target.list.all)))

tmp1=target.df.all[target.df.all$regulon %in% c("Ets2(+)","Nfe2l2(+)"),]

#########to node and edege
edges <- target.df.all %>% 
  rename(from = source.gene, to = target.gene, relation = relation)

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
write.csv(tmp,file="fat.mac.edge.csv",fileEncoding = "UTF-8",row.names = F)
nodes.mac <- data.frame(
  id = unique(c(edges.mac$from, edges.mac$to)),
  label = unique(c(edges.mac$from, edges.mac$to))
)
nodes.mac$type="genes"
nodes.mac$type[1]="Celltype"

###match up tfs in nodes.mac, then update this data
nodes.mac$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="Macrophage" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.mac[,1])]="up"
nodes.mac$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="Macrophage" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.mac[,1])]="down"
#nodes.mac$type[2:(1+length(unique(edges.mac.up$regulon)))]="TFs"
write.csv(nodes.mac,file="fat.mac.node.csv",fileEncoding = "UTF-8",row.names = F)

#####bcell
edges.bcell=edges[edges$celltype=="B_cell",]
which(colnames(edges.bcell)=="direction")
tmp=edges.bcell[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
tmp[tmp[,3]=="regulate",4]=""
write.csv(tmp,file="fat.bcell.edge.csv",fileEncoding = "UTF-8",row.names = F)
nodes.bcell <- data.frame(
  id = unique(c(edges.bcell$from, edges.bcell$to)),
  label = unique(c(edges.bcell$from, edges.bcell$to))
)
nodes.bcell$type="genes"
nodes.bcell$type[1]="Celltype"

nodes.bcell$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="B_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.bcell[,1])]="up"
nodes.bcell$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="B_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.bcell[,1])]="down"
#nodes.bcell$type[2:(1+length(unique(edges.bcell.up$regulon)))]="TFs"
write.csv(nodes.bcell,file="fat.bcell.node.csv",fileEncoding = "UTF-8",row.names = F)

################################################## tcell
edges.tcell=edges[edges$celltype=="T_cell",]
which(colnames(edges.tcell)=="direction")
tmp=edges.tcell[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
tmp[tmp[,3]=="regulate",4]=""
write.csv(tmp,file="fat.tcell.edge.csv",fileEncoding = "UTF-8",row.names = F)
nodes.tcell <- data.frame(
  id = unique(c(edges.tcell$from, edges.tcell$to)),
  label = unique(c(edges.tcell$from, edges.tcell$to))
)
nodes.tcell$type="genes"
nodes.tcell$type[1]="Celltype"

nodes.tcell$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="T_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.tcell[,1])]="up"
nodes.tcell$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="T_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.tcell[,1])]="down"
#nodes.tcell$type[2:(1+length(unique(edges.tcell.up$regulon)))]="TFs"
write.csv(nodes.tcell,file="fat.tcell.node.csv",fileEncoding = "UTF-8",row.names = F)

###############################################################
edges.dc=edges[edges$celltype=="DC",]
which(colnames(edges.dc)=="direction")
tmp=edges.dc[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
tmp[tmp[,3]=="regulate",4]=""
write.csv(tmp,file="fat.dc.edge.csv",fileEncoding = "UTF-8",row.names = F)
nodes.dc <- data.frame(
  id = unique(c(edges.dc$from, edges.dc$to)),
  label = unique(c(edges.dc$from, edges.dc$to))
)
nodes.dc$type="genes"
nodes.dc$type[1]="Celltype"

nodes.dc$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="DC" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.dc[,1])]="up"
nodes.dc$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="DC" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.dc[,1])]="down"
#nodes.dc$type[2:(1+length(unique(edges.dc.up$regulon)))]="TFs"
write.csv(nodes.dc,file="fat.dc.node.csv",fileEncoding = "UTF-8",row.names = F)






#####################
tmp=edges.mac.up[,25:27]
colnames(tmp)
colnames(tmp)=c("Source","Target","Type")
write.csv(tmp,file="edges.mac.up.edge.csv",fileEncoding = "UTF-8",row.names = F)

nodes.mac.up <- data.frame(
  id = unique(c(edges.mac.up$from, edges.mac.up$to)),
  label = unique(c(edges.mac.up$from, edges.mac.up$to))
)

tmp1=nodes.mac.up
tmp1$type="genes"
tmp1$type[1]="Celltype"
tmp1$type[2:(1+length(unique(edges.mac.up$regulon)))]="TFs"
write.csv(tmp1,file="nodes.mac.up.csv",fileEncoding = "UTF-8",row.names = F)


##################down
edges.mac.down=edges[edges$celltype=="Macrophage" & edges$direction=="down",]
tmp=edges.mac.down[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
write.csv(tmp,file="fat.mac.edge.down.csv",fileEncoding = "UTF-8",row.names = F)

nodes.mac.down <- data.frame(
  id = unique(c(edges.mac.down$from, edges.mac.down$to)),
  label = unique(c(edges.mac.down$from, edges.mac.down$to))
)
nodes.mac.down$type="genes"
nodes.mac.down$type[1]="Celltype"

nodes.mac.down$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="Macrophage" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.mac.down[,1])]="down"
nodes.mac.down$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="Macrophage" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.mac.down[,1])]="down"
#nodes.mac$type[2:(1+length(unique(edges.mac.down$regulon)))]="TFs"
write.csv(nodes.mac.down,file="fat.mac.node.down.csv",fileEncoding = "UTF-8",row.names = F)


############################b cell down
edges.bcell.down=edges[edges$celltype=="B_cell" & edges$direction=="down",]
tmp=edges.bcell.down[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
write.csv(tmp,file="fat.bcell.edge.down.csv",fileEncoding = "UTF-8",row.names = F)

nodes.bcell.down <- data.frame(
  id = unique(c(edges.bcell.down$from, edges.bcell.down$to)),
  label = unique(c(edges.bcell.down$from, edges.bcell.down$to))
)
nodes.bcell.down$type="genes"
nodes.bcell.down$type[1]="Celltype"

nodes.bcell.down$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="B_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.bcell.down[,1])]="down"
nodes.bcell.down$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="B_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.bcell.down[,1])]="down"
#nodes.bcell$type[2:(1+length(unique(edges.bcell.down$regulon)))]="TFs"
write.csv(nodes.bcell.down,file="fat.bcell.node.down.csv",fileEncoding = "UTF-8",row.names = F)


######################b cell up
edges.bcell.up=edges[edges$celltype=="B_cell" & edges$direction=="up",]
tmp=edges.bcell.up[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
write.csv(tmp,file="fat.bcell.edge.up.csv",fileEncoding = "UTF-8",row.names = F)

nodes.bcell.up <- data.frame(
  id = unique(c(edges.bcell.up$from, edges.bcell.up$to)),
  label = unique(c(edges.bcell.up$from, edges.bcell.up$to))
)
nodes.bcell.up$type="genes"
nodes.bcell.up$type[1]="Celltype"

nodes.bcell.up$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="B_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.bcell.up[,1])]="up"
nodes.bcell.up$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="B_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.bcell.up[,1])]="up"
#nodes.bcell$type[2:(1+length(unique(edges.bcell.up$regulon)))]="TFs"
write.csv(nodes.bcell.up,file="fat.bcell.node.up.csv",fileEncoding = "UTF-8",row.names = F)

###################t cell up
edges.tcell.up=edges[edges$celltype=="T_cell" & edges$direction=="up",]
tmp=edges.tcell.up[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
write.csv(tmp,file="fat.tcell.edge.up.csv",fileEncoding = "UTF-8",row.names = F)

nodes.tcell.up <- data.frame(
  id = unique(c(edges.tcell.up$from, edges.tcell.up$to)),
  label = unique(c(edges.tcell.up$from, edges.tcell.up$to))
)
nodes.tcell.up$type="genes"
nodes.tcell.up$type[1]="Celltype"

nodes.tcell.up$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="T_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.tcell.up[,1])]="up"
nodes.tcell.up$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="T_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.tcell.up[,1])]="up"
#nodes.tcell$type[2:(1+length(unique(edges.tcell.up$regulon)))]="TFs"
write.csv(nodes.tcell.up,file="fat.tcell.node.up.csv",fileEncoding = "UTF-8",row.names = F)


#####################t cell down
edges.tcell.down=edges[edges$celltype=="T_cell" & edges$direction=="down",]
tmp=edges.tcell.down[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
write.csv(tmp,file="fat.tcell.edge.down.csv",fileEncoding = "UTF-8",row.names = F)

nodes.tcell.down <- data.frame(
  id = unique(c(edges.tcell.down$from, edges.tcell.down$to)),
  label = unique(c(edges.tcell.down$from, edges.tcell.down$to))
)
nodes.tcell.down$type="genes"
nodes.tcell.down$type[1]="Celltype"

nodes.tcell.down$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="T_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.tcell.down[,1])]="down"
nodes.tcell.down$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="T_cell" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.tcell.down[,1])]="down"
#nodes.tcell$type[2:(1+length(unique(edges.tcell.down$regulon)))]="TFs"
write.csv(nodes.tcell.down,file="fat.tcell.node.down.csv",fileEncoding = "UTF-8",row.names = F)


#############DC down
edges.dc.down=edges[edges$celltype=="DC" & edges$direction=="down",]
tmp=edges.dc.down[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
write.csv(tmp,file="fat.dc.edge.down.csv",fileEncoding = "UTF-8",row.names = F)

nodes.dc.down <- data.frame(
  id = unique(c(edges.dc.down$from, edges.dc.down$to)),
  label = unique(c(edges.dc.down$from, edges.dc.down$to))
)
nodes.dc.down$type="genes"
nodes.dc.down$type[1]="Celltype"

nodes.dc.down$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="DC" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.dc.down[,1])]="down"
nodes.dc.down$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="DC" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="down"], nodes.dc.down[,1])]="down"
#nodes.dc$type[2:(1+length(unique(edges.dc.down$regulon)))]="TFs"
write.csv(nodes.dc.down,file="fat.dc.node.down.csv",fileEncoding = "UTF-8",row.names = F)


##################DC up
edges.dc.up=edges[edges$celltype=="DC" & edges$direction=="up",]
tmp=edges.dc.up[,c(25:27,9)]
colnames(tmp)=c("Source","Target","Type","Direction")
write.csv(tmp,file="fat.dc.edge.up.csv",fileEncoding = "UTF-8",row.names = F)

nodes.dc.up <- data.frame(
  id = unique(c(edges.dc.up$from, edges.dc.up$to)),
  label = unique(c(edges.dc.up$from, edges.dc.up$to))
)
nodes.dc.up$type="genes"
nodes.dc.up$type[1]="Celltype"

nodes.dc.up$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="DC" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.dc.up[,1])]="up"
nodes.dc.up$type[match(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$X.2[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype=="DC" & fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$direction=="up"], nodes.dc.up[,1])]="up"
#nodes.dc$type[2:(1+length(unique(edges.dc.up$regulon)))]="TFs"
write.csv(nodes.dc.up,file="fat.dc.node.up.csv",fileEncoding = "UTF-8",row.names = F)


###############normal node color gray [204,204,204] CCCCCC
#up AA0000
#down node color 0066FF
#label color 0,0,0
#sig tf label 204 0 0, size 30
#cell type size 20 labelsize 0.9
#tf size 15 labelsize 0.8
#target sig label 255 0 0 labelsize 0.9
#Levels: Macrophage B_cell T_cell DC NK Mast ILC2 Neutrophil Plasma

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
pdf(file="fat.regulon.mac.vlnplot.pdf",width=3,height=6)
StackedVlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),ident="Macrophage", group.by="group",c("Ets2","Nfe2l2","Irf5"), pt.size=0, cols=c("#FF8080", "#FDCABF","#9B8EB6"))
dev.off()

pdf(file="fat.regulon.mac.vlnplot.raw.pdf",width=9,height=6)
VlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),ident="Macrophage", group.by="group",c("Ets2","Nfe2l2","Irf5"), pt.size=0, cols=c("#FF8080", "#FDCABF","#9B8EB6"))
dev.off()


pdf(file="fat.aucell.boxplot.group.Bcells.pdf",width=4,height=6)
VlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),ident="B_cell",group.by="group",c("AMPK_AUC"),cols=c("#FF8080", "#FDCABF","#9B8EB6"), pt.size=0)+ 
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()

merged_seurat_fat@active.assay="RegulonAUC"

pdf(file="fat.regulon.mac.vlnplot.activaty.pdf",width=4,height=6)
p=StackedVlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),ident="Macrophage",features=c("Ets2(+)","Nfe2l2(+)","Irf5(+)"),pt.size=0, cols=c("#FF8080", "#FDCABF","#9B8EB6"),group.by="group")+
  ylab("Activaty Level")
print(p)
dev.off()

merged_seurat_fat@active.assay="RegulonAUC"
pdf(file="fat.regulon.mac.vlnplot.activaty.raw.pdf",width=9,height=6)
VlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),ident="Macrophage", group.by="group",features=c("Ets2(+)","Nfe2l2(+)","Irf5(+)"), pt.size=0, cols=c("#FF8080", "#FDCABF","#9B8EB6"))
dev.off()

merged_seurat_fat@active.assay="RegulonAUC"
target.group=c("HFD_KO","HFD_WT","ND_KO")
my_comp=list(c("HFD_KO","ND_KO"),c("HFD_KO","HFD_WT"))

pdf(file="fat.regulon.mac.vlnplot.activaty.raw.ets2.pdf",width=4,height=6)
VlnPlot(subset(merged_seurat_fat,subset = group %in% target.group),ident="Macrophage", group.by="group",features=c("Ets2(+)"), pt.size=0, cols=c("#FF8080", "#FDCABF","#9B8EB6"))+
  geom_boxplot(width=.2,col="black",fill="white")+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())+
  NoLegend()
dev.off()


ETS2_expression <- GetAssayData(merged_seurat_fat, assay = "RNA",slot="data")[c("Ets2"), ,drop=F]
merged_seurat_fat=AddMetaData(merged_seurat_fat,t(ETS2_expression))

Ets2_HFD_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Ets2
Ets2_HFD_WT_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Ets2
Ets2_ND_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$Ets2
wilcox.test(Ets2_HFD_KO_mac, Ets2_HFD_WT_mac)
wilcox.test(Ets2_HFD_KO_mac, Ets2_ND_KO_mac) 
mean(Ets2_HFD_KO_mac)
mean(Ets2_HFD_WT_mac)
mean(Ets2_ND_KO_mac)
wilcox.test(Ets2_HFD_KO_mac, Ets2_HFD_WT_mac)$p.value #4.67659e-105
wilcox.test(Ets2_HFD_KO_mac, Ets2_ND_KO_mac)$p.value #1.485113e-118

#####ets2 activaty
ETS2_activaty <- GetAssayData(merged_seurat_fat, assay = "RegulonAUC")[c("Ets2(+)"), ,drop=F]
merged_seurat_fat=AddMetaData(merged_seurat_fat,t(ETS2_activaty))

Ets2_activaty_HFD_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$`Ets2(+)`
Ets2_activaty_HFD_WT_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$`Ets2(+)`
Ets2_activaty_ND_KO_mac <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$`Ets2(+)`
wilcox.test(Ets2_activaty_HFD_KO_mac, Ets2_activaty_HFD_WT_mac)
wilcox.test(Ets2_activaty_HFD_KO_mac, Ets2_activaty_ND_KO_mac) 
mean(Ets2_activaty_HFD_KO_mac)
mean(Ets2_activaty_HFD_WT_mac)
mean(Ets2_activaty_ND_KO_mac)
wilcox.test(Ets2_activaty_HFD_KO_mac, Ets2_activaty_HFD_WT_mac)$p.value #3.058177e-264
wilcox.test(Ets2_activaty_HFD_KO_mac, Ets2_activaty_ND_KO_mac)$p.value #2.347772e-158

###################################
################################## dotplot

for (i in unique(fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype)) {
  eval(parse(text = paste0("fat.tf.common.",i,"=fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common[fat.results_ct.hfdko.ndko.and.hfdko.hfdwt.common$celltype==\"",i,"\",]")))
  eval(parse(text = paste0("fat.tf.common.",i,"$delta.mean=(fat.tf.common.",i,"$delta + fat.tf.common.",i,"$delta.1) / 2")))
  eval(parse(text = paste0("fat.tf.common.",i,"$logfoldchange1=log2(fat.tf.common.",i,"$mean_HFD_KO /fat.tf.common.",i,"$mean_ND_KO)")))
  eval(parse(text = paste0("fat.tf.common.",i,"$logfoldchange2=log2(fat.tf.common.",i,"$mean_HFD_KO / fat.tf.common.",i,"$mean_HFD_WT)")))
  eval(parse(text = paste0("fat.tf.common.",i,"$logfc.mean=(fat.tf.common.",i,"$logfoldchange1 + fat.tf.common.",i,"$logfoldchange2) / 2")))
  
}
#######################################################delta
library(ggplot2)
library(ggrepel)  # 用于智能标签防重叠
library(dplyr)    # 用于数据处理

df <- data.frame(
  Regulon = fat.tf.common.T_cell$regulon,
  Delta_KO_ND = fat.tf.common.T_cell$delta,       # 横坐标
  Delta_KO_WT = fat.tf.common.T_cell$delta.1,     # 纵坐标
  Direction = factor(fat.tf.common.T_cell$direction, levels = c("up", "down")),
  Delta_Mean = fat.tf.common.T_cell$delta.mean,     # 新增delta.mean列
  logfc.mean = fat.tf.common.T_cell$logfc.mean
)

data_range <- range(c(df$Delta_KO_ND, df$Delta_KO_WT))
max_abs <- max(abs(data_range)) * 1.05  # 扩展20%边界
limit <- c(-max_abs, max_abs)

# 1. 颜色映射设置
color_mapping <- c("up" = "#AA0000", "down" = "#0066FF")

# 2. 标记delta.mean极值点
# 找出delta.mean最高和最低的各5个基因
top_high <- df %>%
  arrange(desc(Delta_Mean)) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "High")

top_low <- df %>%
  arrange(Delta_Mean) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "Low")

extreme_points <- bind_rows(top_high, top_low)

# 3. 创建基础散点图
p <- ggplot(df, aes(x = Delta_KO_ND, y = Delta_KO_WT)) +
  # 根据方向添加彩色点
  geom_point(aes(color = Direction), size = 3, alpha = 0.7) +
  
  # 应用自定义颜色
  scale_color_manual(values = color_mapping) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # 添加坐标轴线（交汇于0点）[3](@ref)
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  coord_equal(ratio = 1, xlim = limit, ylim = limit) +
  
  # 坐标轴标签
  labs(
    x = "Delta in HFD.KO vs ND.KO",
    y = "Delta in HFD.KO vs HFD.WT",
    color = "Regulation Direction"
  ) +
  
  # 4. 移除网格线和背景
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # 移除主要网格线[6,7](@ref)
    panel.grid.minor = element_blank(),   # 移除次要网格线[6,7](@ref)
    panel.background = element_blank(),   # 移除背景色
    axis.line = element_line(color = "black")  # 保留坐标轴线[7](@ref)
  )

# 5. 添加极值点标签（使用防重叠算法）
p_labeled <- p +
  geom_point(data = extreme_points, aes(color = Direction), size = 4, shape = 1, stroke = 1.5) +  # 空心圆标记极值点
  ggrepel::geom_text_repel(
    data = extreme_points,
    aes(label = paste0(Regulon, "\n(", round(logfc.mean, 2), ")")),  # 基因名+delta.mean值
    size = 3.5,
    color = "black",
    box.padding = 0.6,  # 标签与点间距
    max.overlaps = 20,
    segment.color = "gray50"  # 连接线颜色
  )

# 显示图形
pdf(file="T_cell_tf.Delta.dotplot.pdf", width = 8, height = 8)
print(p_labeled)
dev.off()
# 保存高质量图像（可选）
ggsave("T_cell_tf.Delta.dotplot.pdf", p_labeled, width = 8, height = 8)


library(ggplot2)
library(ggrepel)  # 用于智能标签防重叠
library(dplyr)    # 用于数据处理

df <- data.frame(
  Regulon = fat.tf.common.Macrophage$regulon,
  Delta_KO_ND = fat.tf.common.Macrophage$delta,       # 横坐标
  Delta_KO_WT = fat.tf.common.Macrophage$delta.1,     # 纵坐标
  Direction = factor(fat.tf.common.Macrophage$direction, levels = c("up", "down")),
  Delta_Mean = fat.tf.common.Macrophage$delta.mean,     # 新增delta.mean列
  logfc.mean = fat.tf.common.Macrophage$logfc.mean
)

data_range <- range(c(df$Delta_KO_ND, df$Delta_KO_WT))
max_abs <- max(abs(data_range)) * 1.05  # 扩展20%边界
limit <- c(-max_abs, max_abs)

# 1. 颜色映射设置
color_mapping <- c("up" = "#AA0000", "down" = "#0066FF")

# 2. 标记delta.mean极值点
# 找出delta.mean最高和最低的各5个基因
top_high <- df %>%
  arrange(desc(Delta_Mean)) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "High")

top_low <- df %>%
  arrange(Delta_Mean) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "Low")

extreme_points <- bind_rows(top_high, top_low)

# 3. 创建基础散点图
p <- ggplot(df, aes(x = Delta_KO_ND, y = Delta_KO_WT)) +
  # 根据方向添加彩色点
  geom_point(aes(color = Direction), size = 3, alpha = 0.7) +
  
  # 应用自定义颜色
  scale_color_manual(values = color_mapping) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # 添加坐标轴线（交汇于0点）[3](@ref)
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  coord_equal(ratio = 1, xlim = limit, ylim = limit) +
  
  # 坐标轴标签
  labs(
    x = "Delta in HFD.KO vs ND.KO",
    y = "Delta in HFD.KO vs HFD.WT",
    color = "Regulation Direction"
  ) +
  
  # 4. 移除网格线和背景
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # 移除主要网格线[6,7](@ref)
    panel.grid.minor = element_blank(),   # 移除次要网格线[6,7](@ref)
    panel.background = element_blank(),   # 移除背景色
    axis.line = element_line(color = "black")  # 保留坐标轴线[7](@ref)
  )

# 5. 添加极值点标签（使用防重叠算法）
p_labeled <- p +
  geom_point(data = extreme_points, aes(color = Direction), size = 4, shape = 1, stroke = 1.5) +  # 空心圆标记极值点
  ggrepel::geom_text_repel(
    data = extreme_points,
    aes(label = paste0(Regulon, "\n(", round(logfc.mean, 2), ")")),  # 基因名+delta.mean值
    size = 3.5,
    color = "black",
    box.padding = 0.6,  # 标签与点间距
    max.overlaps = 20,
    segment.color = "gray50"  # 连接线颜色
  )

# 显示图形
pdf(file="Macrophage_tf.Delta.dotplot.pdf", width = 8, height = 8)
print(p_labeled)
dev.off()
# 保存高质量图像（可选）
ggsave("Macrophage_tf.Delta.dotplot.pdf", p_labeled, width = 8, height = 8)


library(ggplot2)
library(ggrepel)  # 用于智能标签防重叠
library(dplyr)    # 用于数据处理

df <- data.frame(
  Regulon = fat.tf.common.DC$regulon,
  Delta_KO_ND = fat.tf.common.DC$delta,       # 横坐标
  Delta_KO_WT = fat.tf.common.DC$delta.1,     # 纵坐标
  Direction = factor(fat.tf.common.DC$direction, levels = c("up", "down")),
  Delta_Mean = fat.tf.common.DC$delta.mean,     # 新增delta.mean列
  logfc.mean = fat.tf.common.DC$logfc.mean
)

data_range <- range(c(df$Delta_KO_ND, df$Delta_KO_WT))
max_abs <- max(abs(data_range)) * 1.05  # 扩展20%边界
limit <- c(-max_abs, max_abs)

# 1. 颜色映射设置
color_mapping <- c("up" = "#AA0000", "down" = "#0066FF")

# 2. 标记delta.mean极值点
# 找出delta.mean最高和最低的各5个基因
top_high <- df %>%
  arrange(desc(Delta_Mean)) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "High")

top_low <- df %>%
  arrange(Delta_Mean) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "Low")

extreme_points <- bind_rows(top_high, top_low)

# 3. 创建基础散点图
p <- ggplot(df, aes(x = Delta_KO_ND, y = Delta_KO_WT)) +
  # 根据方向添加彩色点
  geom_point(aes(color = Direction), size = 3, alpha = 0.7) +
  
  # 应用自定义颜色
  scale_color_manual(values = color_mapping) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # 添加坐标轴线（交汇于0点）[3](@ref)
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  coord_equal(ratio = 1, xlim = limit, ylim = limit) +
  
  # 坐标轴标签
  labs(
    x = "Delta in HFD.KO vs ND.KO",
    y = "Delta in HFD.KO vs HFD.WT",
    color = "Regulation Direction"
  ) +
  
  # 4. 移除网格线和背景
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # 移除主要网格线[6,7](@ref)
    panel.grid.minor = element_blank(),   # 移除次要网格线[6,7](@ref)
    panel.background = element_blank(),   # 移除背景色
    axis.line = element_line(color = "black")  # 保留坐标轴线[7](@ref)
  )

# 5. 添加极值点标签（使用防重叠算法）
p_labeled <- p +
  geom_point(data = extreme_points, aes(color = Direction), size = 4, shape = 1, stroke = 1.5) +  # 空心圆标记极值点
  ggrepel::geom_text_repel(
    data = extreme_points,
    aes(label = paste0(Regulon, "\n(", round(logfc.mean, 2), ")")),  # 基因名+delta.mean值
    size = 3.5,
    color = "black",
    box.padding = 0.6,  # 标签与点间距
    max.overlaps = 20,
    segment.color = "gray50"  # 连接线颜色
  )

# 显示图形
pdf(file="DC_tf.Delta.dotplot.pdf", width = 8, height = 8)
print(p_labeled)
dev.off()
# 保存高质量图像（可选）
ggsave("DC_tf.Delta.dotplot.pdf", p_labeled, width = 8, height = 8)

library(ggplot2)
library(ggrepel)  # 用于智能标签防重叠
library(dplyr)    # 用于数据处理

df <- data.frame(
  Regulon = fat.tf.common.B_cell$regulon,
  Delta_KO_ND = fat.tf.common.B_cell$delta,       # 横坐标
  Delta_KO_WT = fat.tf.common.B_cell$delta.1,     # 纵坐标
  Direction = factor(fat.tf.common.B_cell$direction, levels = c("up", "down")),
  Delta_Mean = fat.tf.common.B_cell$delta.mean,     # 新增delta.mean列
  logfc.mean = fat.tf.common.B_cell$logfc.mean
)

data_range <- range(c(df$Delta_KO_ND, df$Delta_KO_WT))
max_abs <- max(abs(data_range)) * 1.05  # 扩展20%边界
limit <- c(-max_abs, max_abs)

# 1. 颜色映射设置
color_mapping <- c("up" = "#AA0000", "down" = "#0066FF")

# 2. 标记delta.mean极值点
# 找出delta.mean最高和最低的各5个基因
top_high <- df %>%
  arrange(desc(Delta_Mean)) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "High")

top_low <- df %>%
  arrange(Delta_Mean) %>%
  slice_head(n = 5) %>%
  mutate(Label_Type = "Low")

extreme_points <- bind_rows(top_high, top_low)

# 3. 创建基础散点图
p <- ggplot(df, aes(x = Delta_KO_ND, y = Delta_KO_WT)) +
  # 根据方向添加彩色点
  geom_point(aes(color = Direction), size = 3, alpha = 0.7) +
  
  # 应用自定义颜色
  scale_color_manual(values = color_mapping) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # 添加坐标轴线（交汇于0点）[3](@ref)
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  coord_equal(ratio = 1, xlim = limit, ylim = limit) +
  
  # 坐标轴标签
  labs(
    x = "Delta in HFD.KO vs ND.KO",
    y = "Delta in HFD.KO vs HFD.WT",
    color = "Regulation Direction"
  ) +
  
  # 4. 移除网格线和背景
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),   # 移除主要网格线[6,7](@ref)
    panel.grid.minor = element_blank(),   # 移除次要网格线[6,7](@ref)
    panel.background = element_blank(),   # 移除背景色
    axis.line = element_line(color = "black")  # 保留坐标轴线[7](@ref)
  )

# 5. 添加极值点标签（使用防重叠算法）
p_labeled <- p +
  geom_point(data = extreme_points, aes(color = Direction), size = 4, shape = 1, stroke = 1.5) +  # 空心圆标记极值点
  ggrepel::geom_text_repel(
    data = extreme_points,
    aes(label = paste0(Regulon, "\n(", round(logfc.mean, 2), ")")),  # 基因名+delta.mean值
    size = 3.5,
    color = "black",
    box.padding = 0.6,  # 标签与点间距
    max.overlaps = 20,
    segment.color = "gray50"  # 连接线颜色
  )

# 显示图形
pdf(file="B_cell_tf.Delta.dotplot.pdf", width = 8, height = 8)
print(p_labeled)
dev.off()
# 保存高质量图像（可选）
ggsave("B_cell_tf.Delta.dotplot.pdf", p_labeled, width = 8, height = 8)
##############################################################################################################
##############################################################################################################
##############################################################################################################
##############################################################################################################
#####number plot
library(dplyr)
library(tidyr)
tmp=results_ct.hfdko.ndko.and.hfdko.hfdwt.common[,1:10]
tmp=tmp[!is.na(tmp[,1]),]

summary_df.tf.common <- tmp %>%
  group_by(celltype, direction) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(
    celltype = reorder(celltype, count, FUN = sum)  # Rank by total genes (descending)
  )
summary_df.tf.common

pdf(file="fat.pancelltype.overlap.detfs.number.counts.pdf",width=6,height=6)
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


















###
results_ct.hfdko.ndko.and.hfdko.hfdwt.common[,1]=results_ct.hfdko.ndko.and.hfdko.hfdwt.common[,3]
results_ct.hfdko.ndko.and.hfdko.hfdwt.common[,1]=gsub("\\(\\+\\)","",results_ct.hfdko.ndko.and.hfdko.hfdwt.common[,1])

colnames(results_ct.hfdko.ndko.and.hfdko.hfdwt.common)[1] <- "genesymbol"
colnames(results_ct.hfdko.ndko.and.hfdko.hfdwt.common)[colnames(results_ct.hfdko.ndko.and.hfdko.hfdwt.common) == "celltype"] <- "celllineage"
head(degssig_overlap)  # 确认含genesymbol和celllineage列
head(results_ct.hfdko.ndko.and.hfdko.hfdwt.common)


colnames(results_ct.hfdko.ndko.and.hfdko.hfdwt.common)[13]="celllineage.1"
colnames(results_ct.hfdko.ndko.and.hfdko.hfdwt.common)[19]="direction.1"
merged_data <- merge(
  x = degssig_overlap,
  y = results_ct.hfdko.ndko.and.hfdko.hfdwt.common,
  by = c("genesymbol", "celllineage","direction"),  # 按基因和细胞类型匹配
  all = TRUE                           # 保留所有行（类似全外连接）
)

merged_data=merged_data[!is.na(merged_data$target),]
merged_data=merged_data[!is.na(merged_data$p_val_HFD_KOvsHFD_WT),]


###tf indegs plot
###########################overlap mac up only go and kegg
library(ggplot2)
tmp=readxl::read_xlsx("./mac_up/metascape_result.xlsx",sheet = 2)
tmp=as.data.frame(tmp)
tmp=tmp[grep("Summary",tmp[,1]),]
tmp$Description=gsub(" - Mus musculus \\(house mouse\\)","",tmp$Description)
tmp$Pathways=paste(tmp$Term,tmp$Description,sep=": ")
tmp$neglogP=-tmp$LogP
tmp$neglogQ=-tmp$`Log(q-value)`
tmp=head(tmp,20)
head(tmp)

tmp$Pathways[12]
tmp$Pathways[12]="GO:1901224: positive regulation of \n non-canonical NF-kappaB signal transduction"

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
dev.copy2pdf(file="fat.mac.metascape.up.kegg.go.top20.pdf",width=9,height=0.3*(nrow(tmp)+1))







###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
#####fat m1 m2
pdf(file="fat.m2.Egr2.featueplot.pdf",width=9,height=9)
FeaturePlot(merged_seurat_fat,features="Egr2",raster=FALSE,min.cutoff =0,reduction = "umap",cols= c("gray90", "red"))
dev.off()

pdf(file="fat.m2.Fcgr2b.featueplot.pdf",width=9,height=9)
FeaturePlot(merged_seurat_fat,features="Fcgr2b",raster=FALSE,min.cutoff =0,reduction = "umap",cols= c("gray90", "red"))
dev.off()

# m2.markers=c("Arg1","Arg2","Il10","Fcgr2b","Cd163", "Fcer2a","Cd200r1","Pdcd1lg2",
#              "Cd274", "Marco", "Csf1r", "Mrc1","Il1rn", "Il1r2", "Il4ra", "Ccl4", 
#              "Ccl20", "Ccl17", "Ccl22", "Ccl24", "Lyve1", "Vegfa", "Vegfb", "Vegfc",
#              "Vegfd", "Egf", "Ctsa","Ctsb","Cst3","Ctsd","Tgfb1", "Tgfb2",
#              "Tgfb3", "Mmp14", "Mmp19", "Mmp9","Clec7a","Wnt7b", "Fasl","Tnfsf12" ,
#              "Tnfsf8","Cd276", "Vtcn1", "Msr1","Fn1", "Irf4")

m2.markers=c("Il4ra","Ccl4","Ccl20","Ccl17","Ccl22","Ccl24","Lyve1","Vegfa","Vegfb",
             "Vegfc","Vegfd","Egf","Ctsa","Ctsb","Ctsc","Ctsd","Tgfb1","Tgfb2","Tgfb3",
             "Mmp14","Mmp19","Mmp9","Clec7a","Wnt7b","Fasl","Tnfsf12","Tnfsf8","Cd276",
             "Vtcn1","Msr1","Fn1","Irf4","Cd163","Egr2")

png(file="fat.m2.markers1.featueplot.png",width=1600,height=1600)
FeaturePlot(merged_seurat_fat,features=m2.markers[1:16],raster=FALSE,min.cutoff =0,reduction = "umap",cols= c("gray90", "red"),ncol=4)
dev.off()

png(file="fat.m2.markers2.featueplot.png",width=1600,height=1600)
FeaturePlot(merged_seurat_fat,features=m2.markers[17:32],raster=FALSE,min.cutoff =0,reduction = "umap",cols= c("gray90", "red"),ncol=4)
dev.off()
# png(file="fat.m2.markers3.featueplot.png",width=1600,height=1600)
# FeaturePlot(merged_seurat_fat,features=m2.markers[33:46],raster=FALSE,min.cutoff =0,reduction = "umap",cols= c("gray90", "red"),ncol=4)
# dev.off()

# m1.markers=c("Nos2","Il12a","Il12b","Fcgr1","Cd80","Il23a","Il12b","Cxcl9",
#              "Cxcl10","Cxcl11","Cd86","Il1a","Il1b","Il6", "Tnf", "H2???Aa","H2???Ab1", "H2???Eb1",
#              "Ccl5","Irf5","Irf1","Cd40","Ido1","Kynu","Ccr7")
m1.markers=c("Il23a","Tnf","Cxcl9","Cxcl10","Cxcl11","Cd86","Il1a","Il1b","Il6","Ccl5","Irf5","Irf1",
             "Cd40","Ido1","Kynu","Ccr7")

png(file="fat.m1.markers1.featueplot.png",width=1600,height=1200)
FeaturePlot(merged_seurat_fat,features=m1.markers[1:16],raster=FALSE,min.cutoff =0,reduction = "umap",cols= c("gray90", "red"),ncol=4)
dev.off()

# png(file="fat.m1.markers2.featueplot.png",width=1600,height=800)
# FeaturePlot(merged_seurat_fat,features=m1.markers[17:27],raster=FALSE,min.cutoff =0,reduction = "umap",cols= c("gray90", "red"),ncol=4)
# dev.off()

########################################dotplot
########################################
fat.mac@active.ident=fat.mac$celltype
pdf(file="fat.mac.m2.markers1.dotplot.pdf",width=16,height=6)
DotPlot(fat.mac,features=m2.markers[1:16],cols= c("gray90", "red"))
dev.off()
pdf(file="fat.mac.m2.markers2.dotplot.pdf",width=16,height=6)
DotPlot(fat.mac,features=m2.markers[17:34],cols= c("gray90", "red"))
dev.off()
pdf(file="fat.mac.m1.markers1.dotplot.pdf",width=16,height=6)
DotPlot(fat.mac,features=m1.markers[1:16],cols= c("gray90", "red"))
dev.off()


####################################################
geneSets.m1=list(
  m1.genes = m1.markers  # ?????б?Ԫ??
)
set.seed(100)
cells_AUC <- AUCell_calcAUC(
  geneSets.m1, 
  cells_rankings,
  aucMaxRank = 950,
)

m1_scores <- getAUC(cells_AUC)["m1.genes", ]
merged_seurat_fat$m1_scores <- m1_scores

####m2
geneSets.m2=list(
  m2.genes = m2.markers  # ?????б?Ԫ??
)
set.seed(100)
cells_AUC <- AUCell_calcAUC(
  geneSets.m2, 
  cells_rankings,
  aucMaxRank = 950
)

m2_scores <- getAUC(cells_AUC)["m2.genes", ]
merged_seurat_fat$m2_scores <- m2_scores



fat.mac <- subset(merged_seurat_fat, 
                  subset = celllineage %in% c("Macrophage"))
##################################################m1 m2 score plot
library(ggrepel)
m1m2avg_scores <- fat.mac@meta.data %>%
  group_by(celltype) %>%
  summarise(
    avg_m1 = mean(m1_scores, na.rm = TRUE),
    avg_m2 = mean(m2_scores, na.rm = TRUE)
  )

# ????ɢ??ͼ
p <- ggplot(m1m2avg_scores, aes(x = avg_m1, y = avg_m2, color = celltype)) +
  geom_point(size = 8, alpha = 0.8) +
  geom_text_repel(aes(label = celltype), size = 5, point.padding = 0.8) +
  scale_color_brewer(palette = "Set1") +
  labs(x = "M1 Score", y = "M2 Score", 
       title = "Macrophage Subtype Polarization") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))
pdf(file="m1m2score.dotplot.pdf",width=8,height=8)
print(p)
dev.off()


#######################################################
########## m1 m2 ratio
#######################################################

library(tidyverse)
meta_data <- fat.mac@meta.data %>%
  dplyr::select(sample, group, celltype)

cell_counts <- meta_data %>%
  filter(celltype %in% c("Macrophage_Ear2", "Macrophage_Gas6")) %>%
  group_by(sample, group, celltype) %>%
  tally(name = "cell_num") %>%
  tidyr::spread(key = celltype, value = cell_num, fill = 0)  # ????0????NA

# ????Ear2/Ednrb????
ratio_data <- cell_counts %>%
  mutate(ratio = Macrophage_Ear2 / Macrophage_Gas6) %>%
  dplyr::select(sample, group, ratio)

# ?鿴????
head(ratio_data)

####ratio test 3 groups
cols=c("#FF8080", "#FDCABF","#9B8EB6")
pdf(file="hfdko_hfdwt_ndko.m1m2ratio.pdf",width=6,height=5)
ggplot(ratio_data[1:9,], aes(x = group, y = ratio, fill = group)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("HFD_KO" = "#FF8080", "HFD_WT" = "#FDCABF","ND_KO" = "#9B8EB6")) +
  labs(x = NULL, y = "Macrophage_Ear2 / Macrophage_Gas6 Ratio",
       title = "") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5))
dev.off()

hfd_data <- ratio_data %>% 
  filter(group %in% c("HFD_KO", "HFD_WT"))

wilcox_test <- wilcox.test(ratio ~ group, data = hfd_data)
p_value <- wilcox_test$p.value

# ??????????ֵ????
mean_diff <- hfd_data %>%
  group_by(group) %>%
  summarise(mean_ratio = mean(ratio)) %>%
  spread(key = group, value = mean_ratio) %>%
  mutate(diff_ratio = HFD_KO - HFD_WT)

# ????????
cat("HFD_KO vs. HFD_WT: p =", round(p_value, 4), 
    "| Mean Ratio Diff =", round(mean_diff$diff_ratio, 2))

pdf(file="hfdko_hfdwt.m1m2ratio.pdf",width=6,height=6)
ggplot(hfd_data, aes(x = group, y = ratio, fill = group)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("HFD_KO" = "#FF8080", "HFD_WT" = "#FDCABF")) +
  labs(x = NULL, y = "Macrophage_Ear2 / Macrophage_Gas6 Ratio",
       title = "diff=0.14, p=0.4",
       caption = paste("Mean Diff =", round(mean_diff$diff_ratio, 2))) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5))
dev.off()



# ??ȡHFD??????
hfd_data <- ratio_data %>% 
  filter(group %in% c("HFD_KO", "ND_KO"))

# Wilcoxon?Ⱥͼ??飨?ǲ??????飩
wilcox_test <- wilcox.test(ratio ~ group, data = hfd_data)
p_value <- wilcox_test$p.value

# ??????????ֵ????
mean_diff <- hfd_data %>%
  group_by(group) %>%
  summarise(mean_ratio = mean(ratio)) %>%
  spread(key = group, value = mean_ratio) %>%
  mutate(diff_ratio = HFD_KO - ND_KO)

# ????????
cat("HFD_KO vs. ND_KO: p =", round(p_value, 4), 
    "| Mean Ratio Diff =", round(mean_diff$diff_ratio, 2))

pdf(file="hfdko_ndko.m1m2ratio.pdf",width=6,height=6)
ggplot(hfd_data, aes(x = group, y = ratio, fill = group)) +
  geom_boxplot(width = 0.6, alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.8) +
  scale_fill_manual(values = c("HFD_KO" = "#FF8080", "ND_KO" = "#FDCABF")) +
  labs(x = NULL, y = "Macrophage_Ear2 / Macrophage_Gas6 Ratio",
       title = "diff=0.15, p=0.4",
       caption = paste("Mean Diff =", round(mean_diff$diff_ratio, 2))) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5))
dev.off()



##################
library("scProportionTest")
mac.prop_test_fat <- sc_utils(fat.mac)
mac.prop_test_fat_KO_HFD_vs_ND <- permutation_test(
  mac.prop_test_fat, cluster_identity = "celltype",
  sample_1 = "ND_KO", sample_2 = "HFD_KO",
  sample_identity = "group"
)

mac.prop_test_fat_HFD_KO_vs_WT <- permutation_test(
  mac.prop_test_fat, cluster_identity = "celltype",
  sample_1 = "HFD_WT", sample_2 = "HFD_KO",
  sample_identity = "group"
)

mac.prop_test_fat_WT_HFD_vs_ND <- permutation_test(
  mac.prop_test_fat, cluster_identity = "celltype",
  sample_1 = "ND_WT", sample_2 = "HFD_WT",
  sample_identity = "group"
)

pdf(file="./mac.prop_test_fat_KO_HFD_vs_ND.pdf",width=6,height=9)
permutation_plot(mac.prop_test_fat_KO_HFD_vs_ND,log2FD_threshold=0)
dev.off()

pdf(file="./mac.prop_test_fat_WT_HFD_vs_ND.pdf",width=6,height=9)
permutation_plot(mac.prop_test_fat_WT_HFD_vs_ND,log2FD_threshold=0)
dev.off()

pdf(file="./mac.prop_test_fat_HFD_KO_vs_WT.pdf",width=6,height=9)
permutation_plot(mac.prop_test_fat_HFD_KO_vs_WT,log2FD_threshold=0)
dev.off()

pdf(file="./mac.prop_test_fat_KO_HFD_vs_ND.0.25.pdf",width=6,height=9)
permutation_plot(mac.prop_test_fat_KO_HFD_vs_ND,log2FD_threshold=0.25)
dev.off()

pdf(file="./mac.prop_test_fat_WT_HFD_vs_ND.0.25.pdf",width=6,height=9)
permutation_plot(mac.prop_test_fat_WT_HFD_vs_ND,log2FD_threshold=0.25)
dev.off()

pdf(file="./mac.prop_test_fat_HFD_KO_vs_WT.0.25.pdf",width=6,height=9)
permutation_plot(mac.prop_test_fat_HFD_KO_vs_WT,log2FD_threshold=0.25)
dev.off()


pdf(file="fat.pdcd1.featureplot.pdf",width=16,height=4)
FeaturePlot(merged_seurat_fat,features="Pdcd1",split.by="group",raster=F,ncol=2,min.cutoff =0,reduction = "umap",cols= c("gray90", "red"),label = T)
dev.off()

pdf(file="fat.pdcd1.vlnplot.pdf",width=16,height=4)
VlnPlot(merged_seurat_fat,features="Pdcd1",split.by="group",group.by="celllineage",pt.size=0,cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"))
dev.off()


pdf(file="int.pdcd1.featureplot.pdf",width=16,height=4)
FeaturePlot(merged_seurat_fat,features="Pdcd1",split.by="group",raster=F,ncol=2,min.cutoff =0,reduction = "umap",cols= c("gray90", "red"),label = T)
dev.off()

pdf(file="int.pdcd1.vlnplot.pdf",width=16,height=4)
VlnPlot(merged_seurat_fat,features="Pdcd1",split.by="group",group.by="celllineage",pt.size=0,cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"))
dev.off()

my_comp=list(c("HFD_KO","ND_KO"),c("HFD_KO","HFD_WT"),c("HFD_WT","ND_WT"),c("ND_KO","ND_WT"))

pdf(file="int.pdcd1.vlnplot.pdf",width=16,height=4)
VlnPlot(merged_seurat_fat,features="Pdcd1",split.by="group",group.by="celllineage",pt.size=0,cols=c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA"))+ 
  stat_compare_means(method = "wilcox.test",comparisons = my_comp,label = "p.signif",tip.length = 0,size = 6)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(axis.title.x = element_blank())
dev.off()


df=FetchData(merged_seurat_fat, vars = c("Pdcd1","celllineage","group"))
result <- df %>%
  group_by(group, celllineage) %>%
  summarise(
    Pdcd1_mean = mean(Pdcd1, na.rm = TRUE),
    .groups = "drop"  # ????????????
  )

# ?鿴????
print(result)
write.csv(result,file="fat.pdcd1.meanexpr.csv")



merged_seurat_fat@meta.data$combined_group <- paste(
  merged_seurat_fat@meta.data$celllineage, 
  merged_seurat_fat@meta.data$group, 
  sep = "_"
)

pdf(file="fat.pdcd1.dotplot.pdf",width=6,height=12)
DotPlot(merged_seurat_fat,features="Pdcd1",group.by="combined_group")
dev.off()


#####################################Pdcd1 prop
pdcd1_expr <- GetAssayData(merged_seurat_fat, assay = "RNA", slot = "data")["Pdcd1", ]
merged_seurat_fat$pdcd1_positive <- pdcd1_expr > 0
library(dplyr)
prop_df.pdcd1 <- merged_seurat_fat@meta.data %>%
  group_by(celllineage, group, sample) %>%         # ??ϸ??????+????+????????????
  summarise(
    total_cells = n(),
    positive_cells = sum(pdcd1_positive),
    prop = positive_cells / total_cells,
    .groups = 'drop'
  ) %>%
  ungroup()

my_comp <- list(c("HFD_KO", "ND_KO"), c("HFD_KO", "HFD_WT"), 
                c("HFD_WT", "ND_WT"), c("ND_KO", "ND_WT"))
library(ggplot2)

pdf(file="prop_df.pdcd1.barplot.compare.pdf",width=12,height=8)
ggplot(prop_df.pdcd1, aes(x = group, y = prop, fill = group)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.7) +        # ??????ֵ??????״ͼ
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6, color = "gray20") + # ??????????
  scale_fill_manual(values = c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA")) +   # Ӧ???Զ?????ɫ
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ celllineage, scales = "free_y", ncol = 3) +  # ??ϸ?????ͷ???
  labs(x = "Sample Group", y = "Adipose Pdcd1+ Cells Proportion", fill = "Group") +
  
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )
dev.off()

################################################################################
#################### expr > 1
pdcd1_expr <- GetAssayData(merged_seurat_fat, assay = "RNA", slot = "data")["Pdcd1", ]
merged_seurat_fat$pdcd1_positive1 <- pdcd1_expr > 1
library(dplyr)
prop_df.pdcd1 <- merged_seurat_fat@meta.data %>%
  group_by(celllineage, group, sample) %>%         # ??ϸ??????+????+????????????
  summarise(
    total_cells = n(),
    positive_cells = sum(pdcd1_positive1),
    prop = positive_cells / total_cells,
    .groups = 'drop'
  ) %>%
  ungroup()

my_comp <- list(c("HFD_KO", "ND_KO"), c("HFD_KO", "HFD_WT"), 
                c("HFD_WT", "ND_WT"), c("ND_KO", "ND_WT"))
library(ggplot2)

prop_df.pdcd1=prop_df.pdcd1[prop_df.pdcd1$celllineage != "Neutrophil",]
pdf(file="prop_df.pdcd1.thre1.barplot.compare.pdf",width=12,height=8)
ggplot(prop_df.pdcd1, aes(x = group, y = prop, fill = group)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.7) +        # ??????ֵ??????״ͼ
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.6, color = "gray20") + # ??????????
  scale_fill_manual(values = c("#FF8080", "#FDCABF","#9B8EB6", "#CCD1DA")) +   # Ӧ???Զ?????ɫ
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.1))) +
  
  facet_wrap(~ celllineage, scales = "free_y", ncol = 3) +  # ??ϸ?????ͷ???
  labs(x = "Sample Group", y = "Adipose Pdcd1+ Cells Proportion", fill = "Group") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )
dev.off()





setwd("E:/jiaoda/xianjiaoda/reclusterformac/harmony.dim30.res0.60.afterfilter.then020/")
dir("./metascape/")
fat.all.metascape=dir("./metascape/")
fat.all.metascape
fat.all.metascape=fat.all.metascape[-grep("zip",fat.all.metascape)]
fat.all.metascape.list=list()
for (i in 1:length(fat.all.metascape)) {
  tmp=readxl::read_xlsx(paste0("./metascape/",fat.all.metascape[i],"/metascape_result.xlsx"),sheet = 2)
  tmp=as.data.frame(tmp)
  tmp=tmp[grep("Summary",tmp[,1]),]
  tmp$Description=gsub(" - Mus musculus \\(house mouse\\)","",tmp$Description)
  tmp$Pathways=paste(tmp$Term,tmp$Description,sep=": ")
  tmp$neglogP=-tmp$LogP
  tmp$neglogQ=-tmp$`Log(q-value)`
  tmp=head(tmp,10)
  
  fat.all.metascape.list[[i]]=tmp
  
}
names(fat.all.metascape.list)=fat.all.metascape
#red #B53E2B
#blue #1F78B4

####top5
fat.all.metascape.list.top5=list()
for (i in 1:length(fat.all.metascape)) {
  tmp=readxl::read_xlsx(paste0("./metascape/",fat.all.metascape[i],"/metascape_result.xlsx"),sheet = 2)
  tmp=as.data.frame(tmp)
  tmp=tmp[grep("Summary",tmp[,1]),]
  tmp$Description=gsub(" - Mus musculus \\(house mouse\\)","",tmp$Description)
  tmp$Pathways=paste(tmp$Term,tmp$Description,sep=": ")
  tmp$neglogP=-tmp$LogP
  tmp$neglogQ=-tmp$`Log(q-value)`
  tmp=head(tmp,5)
  
  fat.all.metascape.list.top5[[i]]=tmp
  
}
names(fat.all.metascape.list.top5)=fat.all.metascape

library(ggplot2)
setwd("metascape_plots")
####don't divide by up down
for (i in 1:length(fat.all.metascape.list)) {
  nrows=nrow(fat.all.metascape.list[[i]])
  pdf(file=paste0(fat.all.metascape[i],".pdf"),width = 12,height=(0.6*nrows)+0.6)
  if(length(grep("up",fat.all.metascape[i]))==1){
    p=ggplot(data=fat.all.metascape.list[[i]], aes(x=reorder(Pathways,neglogQ), y=neglogQ)) + #reorder 
      geom_bar(stat="identity", width=0.8,fill="#B53E2B") + coord_flip() + 
      #  scale_x_discrete(labels=labels) +
      xlab("") + ylab("")+
      labs(title = "")+
      theme_bw()+
      theme(axis.text=element_text(face = "bold", color="black",size = 12),
            axis.text.x = element_text(size = 18),  # X???̶?
            axis.text.y = element_text(size = 18),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())}
  else{
    p=ggplot(data=fat.all.metascape.list[[i]], aes(x=reorder(Pathways,neglogQ), y=neglogQ)) + #reorder 
      geom_bar(stat="identity", width=0.8,fill="#1F78B4") + coord_flip() + 
      #  scale_x_discrete(labels=labels) +
      xlab("") + ylab("")+
      labs(title = "")+
      theme_bw()+
      theme(axis.text=element_text(face = "bold", color="black",size = 12),
            axis.text.x = element_text(size = 18),  # X???̶?
            axis.text.y = element_text(size = 18),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())}
  print(p)
  dev.off()
}

fat.mac.metascape.up=fat.all.metascape.list[["fat.mac.up"]]
fat.mac.metascape.down=fat.all.metascape.list[["fat.mac.down"]]


########combine
fat.all.cellname=strsplit(fat.all.metascape,split = "[.]")
fat.all.cellname=sapply(fat.all.cellname,function(x){c(x[[2]])})
fat.all.cellname=unique(fat.all.cellname)
fat.all.cellname
for (i in 1:length(fat.all.metascape)) {
  if(length(grep("up",fat.all.metascape[i]))==1){
    fat.all.metascape.list[[i]]$direction="up"
  }
  else{fat.all.metascape.list[[i]]$direction="down"
  }
}

for (i in fat.all.cellname) {
  eval(parse(text = paste0("fat.",i,".metascape=rbind(fat.all.metascape.list[[\"fat.",i,".up\"]],fat.all.metascape.list[[\"fat.",i,".down\"]])")))
}

fat.all.cellname
fat.all.cellname.full=paste("fat",fat.all.cellname,"metascape",sep=".")
fat.all.cellname.full
for (i in fat.all.cellname.full) {
  pdf(file=paste0("./enrich/",i,".pdf"),width=9,height=0.3*(eval(parse(text = paste0("nrow(",i,")")))+1))
  p=ggplot(data=get(i), aes(x=reorder(Pathways,neglogQ), y=neglogQ,fill=direction)) + #reorder 
    geom_bar(stat="identity", width=0.8) + coord_flip() + 
    scale_fill_manual(values = c(up = "#B53E2B", down = "#1F78B4")) +
    xlab("") + ylab("")+
    labs(title = "")+
    theme_bw()+
    theme(axis.text=element_text(face = "bold", color="black",size = 12),
          axis.text.x = element_text(size = 12),  # X???̶?
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(p)
  dev.off()
  
}
fat.nk.metascape$Pathways
fat.nk.metascape$Pathways[11]
fat.nk.metascape$Pathways[11]="GO:0097011: cellular response to granulocyte macrophage \n colony-stimulating factor stimulus"

#####single
#red #B53E2B
#blue #1F78B4
names(fat.all.metascape.list)
# [1] "Bcells_down" "Bcells_up"   "DC_down"     "DC_up"       "mac_down"    "mac_up"      "NK_down"     "NK_up"       "Tcells_down"
# [10] "Tcells_up"  
fat.all.metascape.list[[3]]
fat.all.metascape.list[[3]]$Pathways
fat.all.metascape.list[[3]]$Pathways[[3]]="GO:1900151: regulation of nuclear-transcribed mRNA \n catabolic process, deadenylation-dependent decay"
fat.all.metascape.list.top5[[3]]
fat.all.metascape.list.top5[[3]]$Pathways
fat.all.metascape.list.top5[[3]]$Pathways[[3]]="GO:1900151: regulation of nuclear-transcribed mRNA \n catabolic process, deadenylation-dependent decay"

library(scales)
p=ggplot(data=fat.all.metascape.list.top5[[3]], aes(x=reorder(Pathways,neglogQ), y=neglogQ)) + #reorder 
  geom_bar(stat="identity", width=0.8,fill="#1F78B4") + coord_flip() + 
  #  scale_x_discrete(labels=labels) +
  xlab("") + ylab("")+
  labs(title = "")+
  theme_bw()+
  theme(axis.text=element_text(face = "bold", color="black",size = 12),
        axis.text.x = element_text(size = 15),  # X???̶?
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),aspect.ratio = 1)+
  scale_y_continuous(labels = label_number(accuracy = 1))
print(p)
dev.copy2pdf(file="DC.down.pdf",width=9,height=0.3*(nrow(fat.all.metascape.list.top5[[3]]))+2)

ggplot(data=fat.all.metascape.list.top5[[6]], aes(x=reorder(Pathways,neglogQ), y=neglogQ)) + #reorder 
  geom_bar(stat="identity", width=0.8,fill="#B53E2B") + coord_flip() + 
  #  scale_x_discrete(labels=labels) +
  xlab("") + ylab("")+
  labs(title = "")+
  theme_bw()+
  theme(axis.text=element_text(face = "bold", color="black",size = 12),
        axis.text.x = element_text(size = 15),  # X???̶?
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ,aspect.ratio = 1)
dev.copy2pdf(file="fat.mac.up.pdf",width=7,height=0.3*(nrow(fat.all.metascape.list.top5[[6]])*(5/nrow(fat.all.metascape.list.top5[[6]])))+2)


for (i in 1:length(fat.all.metascape.list.top5)) {
  fat.all.metascape.list.top5[[i]]=eval(parse(text = paste0("fat.all.metascape.list.top5[[",i,"]][fat.all.metascape.list.top5[[",i,"]]$neglogQ>-log10(0.05),]")))
}

library(scales)
for (i in 1:length(fat.all.metascape.list.top5)) {
  nrows=nrow(fat.all.metascape.list.top5[[i]])
  pdf(file=paste0(fat.all.metascape[i],".pdf"),width = 0.15*max(nchar(fat.all.metascape.list.top5[[i]]$Pathways)),height=0.3*(nrows)*(5/nrows)+2.5)
  if(length(grep("up",fat.all.metascape[i]))==1){
    p=ggplot(data=fat.all.metascape.list.top5[[i]], aes(x=reorder(Pathways,neglogQ), y=neglogQ)) + #reorder 
      geom_bar(stat="identity", width=0.8,fill="#B53E2B") + coord_flip() + 
      #  scale_x_discrete(labels=labels) +
      xlab("") + ylab("")+
      labs(title = "")+
      theme_bw()+
      theme(axis.text=element_text(face = "bold", color="black",size = 12),
            axis.text.x = element_text(size = 15),  # X???̶?
            axis.text.y = element_text(size = 15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),aspect.ratio = 1)+
      scale_y_continuous(labels = label_number(accuracy = 1))}
  else{
    p=ggplot(data=fat.all.metascape.list.top5[[i]], aes(x=reorder(Pathways,neglogQ), y=neglogQ)) + #reorder 
      geom_bar(stat="identity", width=0.8,fill="#1F78B4") + coord_flip() + 
      #  scale_x_discrete(labels=labels) +
      xlab("") + ylab("")+
      labs(title = "")+
      theme_bw()+
      theme(axis.text=element_text(face = "bold", color="black",size = 12),
            axis.text.x = element_text(size = 15),  # X???̶?
            axis.text.y = element_text(size = 15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),aspect.ratio = 1)+
      scale_y_continuous(labels = label_number(accuracy = 1))}
  print(p)
  dev.off()
}
fat.all.metascape.list.top5[[i]]$Pathways
max(nchar(fat.all.metascape.list.top5[[i]]$Pathways))


################################################prop test in fat
prop_test_fat <- sc_utils(merged_seurat_fat)

prop_test_fat_ND_KO_vs_WT <- permutation_test(
  prop_test_fat, cluster_identity = "celltype",
  sample_1 = "ND_WT", sample_2 = "ND_KO",
  sample_identity = "group"
)

prop_test_fat_KO_HFD_vs_ND <- permutation_test(
  prop_test_fat, cluster_identity = "celltype",
  sample_1 = "ND_KO", sample_2 = "HFD_KO",
  sample_identity = "group"
)

prop_test_fat_HFD_KO_vs_WT <- permutation_test(
  prop_test_fat, cluster_identity = "celltype",
  sample_1 = "HFD_WT", sample_2 = "HFD_KO",
  sample_identity = "group"
)

prop_test_fat_WT_HFD_vs_ND <- permutation_test(
  prop_test_fat, cluster_identity = "celltype",
  sample_1 = "ND_WT", sample_2 = "HFD_WT",
  sample_identity = "group"
)


###############################plots
pdf(file="./prop.test/prop_test_fat_ND_KO_vs_WT.0.025.pdf",width=6,height=9)
permutation_plot(prop_test_fat_ND_KO_vs_WT,log2FD_threshold = 0.25)
dev.off()
pdf(file="./prop.test/prop_test_fat_KO_HFD_vs_ND.0.025.pdf",width=6,height=9)
permutation_plot(prop_test_fat_KO_HFD_vs_ND,log2FD_threshold = 0.25)
dev.off()
pdf(file="./prop.test/prop_test_fat_HFD_KO_vs_WT.0.025.pdf",width=6,height=9)
permutation_plot(prop_test_fat_HFD_KO_vs_WT,log2FD_threshold = 0.25)
dev.off()
pdf(file="./prop.test/prop_test_fat_WT_HFD_vs_ND.0.025.pdf",width=6,height=9)
permutation_plot(prop_test_fat_WT_HFD_vs_ND,log2FD_threshold = 0.25)
dev.off()


################################pan
prop_test_fat.pan <- sc_utils(merged_seurat_fat)

prop_test_fat_ND_KO_vs_WT.pan <- permutation_test(
  prop_test_fat.pan, cluster_identity = "celllineage",
  sample_1 = "ND_WT", sample_2 = "ND_KO",
  sample_identity = "group"
)

prop_test_fat_KO_HFD_vs_ND.pan <- permutation_test(
  prop_test_fat.pan, cluster_identity = "celllineage",
  sample_1 = "ND_KO", sample_2 = "HFD_KO",
  sample_identity = "group"
)

prop_test_fat_HFD_KO_vs_WT.pan <- permutation_test(
  prop_test_fat.pan, cluster_identity = "celllineage",
  sample_1 = "HFD_WT", sample_2 = "HFD_KO",
  sample_identity = "group"
)

prop_test_fat_WT_HFD_vs_ND.pan <- permutation_test(
  prop_test_fat.pan, cluster_identity = "celllineage",
  sample_1 = "ND_WT", sample_2 = "HFD_WT",
  sample_identity = "group"
)

pdf(file="./prop.test/prop_test_fat_ND_KO_vs_WT.pan.0.25.pdf",width=6,height=9)
permutation_plot(prop_test_fat_ND_KO_vs_WT.pan,log2FD_threshold = 0.25)
dev.off()

pdf(file="./prop.test/prop_test_fat_KO_HFD_vs_ND.pan.0.25.pdf",width=6,height=9)
permutation_plot(prop_test_fat_KO_HFD_vs_ND.pan,log2FD_threshold = 0.25)
dev.off()

pdf(file="./prop.test/prop_test_fat_HFD_KO_vs_WT.pan.0.25.pdf",width=6,height=9)
permutation_plot(prop_test_fat_HFD_KO_vs_WT.pan,log2FD_threshold = 0.25)
dev.off()

pdf(file="./prop.test/prop_test_fat_WT_HFD_vs_ND.pan.0.25.pdf",width=6,height=9)
permutation_plot(prop_test_fat_WT_HFD_vs_ND.pan,log2FD_threshold = 0.25)
dev.off()



###############################addmodulescore
"mmu04152_AMPK_signaling_pathway_-_Mus_musculus_(house_mouse)"

library(GSEABase)
genes <- geneIds(geneSets[["mmu04152_AMPK_signaling_pathway_-_Mus_musculus_(house_mouse)"]])  # 提取基因符号向量
# 转换为列表格式（AddModuleScore要求输入list）
gene_list <- list(AMPK_addmodule = genes)  # 命名为"AMPK"，后续打分结果列名为"AMPK1"

set.seed(100)
merged_seurat_fat <- AddModuleScore(
  merged_seurat_fat,
  features = gene_list,
  name = "AMPK_addmodule",
  ctrl = 100
) 
AMPK_HFD_KO_Macrophage_addmodule <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$AMPK_addmodule1
AMPK_HFD_WT_Macrophage_addmodule <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "HFD_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$AMPK_addmodule1
AMPK_ND_KO_Macrophage_addmodule <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_KO"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$AMPK_addmodule1
AMPK_ND_WT_Macrophage_addmodule <- merged_seurat_fat@meta.data[which(merged_seurat_fat@meta.data$group == "ND_WT"&merged_seurat_fat@meta.data$celllineage =="Macrophage"),]$AMPK_addmodule1
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_Macrophage_addmodule, AMPK_HFD_WT_Macrophage_addmodule)
# AMPK_HFD_KO vs. AMPK_HFD_WT
wilcox.test(AMPK_HFD_KO_Macrophage_addmodule, AMPK_ND_KO_Macrophage_addmodule)
####dim50, res0.25 sig, 0.01799
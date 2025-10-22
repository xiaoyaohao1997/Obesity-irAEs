library(RColorBrewer)
library(ggrepel)
library(ggplot2)
library(ggsci)
library(dplyr)
library(viridis)
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




################################################################################################
###Code for BMDM Bulk RNA seq data analysis
################################################################################################


#read  bulk RNA data
gene_counts=read.csv("gene_counts.txt",sep = "\t",row.names = 1)
group_list = factor(c(rep("NDWT",3),rep("NDKO",3),rep("HFDWT",3),rep("HFDKO",3), rep("HFDKOMET",3)),levels = c("NDWT","NDKO","HFDWT","HFDKO","HFDKOMET"))
colour_group = rainbow(length(unique(group_list)))
colour = mypal_3[as.numeric(factor(group_list))]

library(edgeR)
dge <- DGEList(counts = gene_counts, group = group_list)
# TMM
dge <- calcNormFactors(dge, method = "TMM")

#normalized counts
tmm_norm_counts <- cpm(dge,log = T) 
tmm_norm_counts_nolog <- cpm(dge,log = F) 

write.table(tmm_norm_counts,file="normalized_counts_tmm.filter.txt",sep="\t",row.names = T,col.names = T,quote = F)

#groups
design <- model.matrix(~0 + group_list)
colnames(design) <- levels(group_list) 
design

# 差异分析
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design)

####set test groups hfdko and hfdwt
cont_hfdko_hfdwt <- makeContrasts(HFDKO - HFDWT, levels = design)
res_hfdko_hfdwt <- glmQLFTest(fit, contrast = cont_hfdko_hfdwt)
resgene_hfdko_hfdwt <- topTags(res_hfdko_hfdwt, n = Inf, adjust.method = "fdr")$table
resgene_hfdko_hfdwt.sig <- subset(resgene_hfdko_hfdwt, FDR < 0.05)

###############hfdko ndko
cont_hfdko_ndko <- makeContrasts(HFDKO - NDKO, levels = design)
res_hfdko_ndko <- glmQLFTest(fit, contrast = cont_hfdko_ndko)
resgene_hfdko_ndko <- topTags(res_hfdko_ndko, n = Inf, adjust.method = "fdr")$table
resgene_hfdko_ndko.sig <- subset(resgene_hfdko_ndko, FDR < 0.05)

###############hfdkomet hfdko
cont_hfdkomet_hfdfko <- makeContrasts(HFDKOMET - HFDKO, levels = design)
res_hfdkomet_hfdfko <- glmQLFTest(fit, contrast = cont_hfdkomet_hfdfko)
resgene_hfdkomet_hfdfko <- topTags(res_hfdkomet_hfdfko, n = Inf, adjust.method = "fdr")$table
resgene_hfdkomet_hfdfko.sig <- subset(resgene_hfdkomet_hfdfko, FDR < 0.05)


head(resgene_hfdko_hfdwt.sig)
resgene_hfdko_hfdwt.sig$direction="down"
resgene_hfdko_hfdwt.sig$direction[resgene_hfdko_hfdwt.sig$logFC>0]="up"
resgene_hfdko_ndko.sig$direction="down"
resgene_hfdko_ndko.sig$direction[resgene_hfdko_ndko.sig$logFC>0]="up"
resgene_hfdkomet_hfdfko.sig$direction="down"
resgene_hfdkomet_hfdfko.sig$direction[resgene_hfdkomet_hfdfko.sig$logFC>0]="up"

#save results
write.csv(resgene_hfdko_hfdwt.sig,file="resgene_hfdko_hfdwt.sig.all.csv")
write.csv(resgene_hfdko_ndko.sig,file="resgene_hfdko_ndko.sig.all.csv")
write.csv(resgene_hfdkomet_hfdfko.sig,file="resgene_hfdkomet_hfdfko.sig.all.csv")
table(resgene_hfdkomet_hfdfko.sig$direction)


#####cytokines heatmap
c("Il1b","Nfkb2","Il12a","Il36a","Tnf","Cxcl10","Ccl4","Cxcl9","Ccl8","Ccl2","Ccl12")
cytokines_df=tmm_norm_counts[match(c("Il1b","Il12a","Il36a","Tnf","Cxcl10","Ccl4","Cxcl9","Ccl8","Ccl2","Ccl12"),row.names(tmm_norm_counts)),]

pheatmap::pheatmap(cytokines_df[,c(10:15)],show_rownames=T,show_colnames=T,cluster_col = FALSE,cluster_row = T,scale="row",treeheight_row = 0,cellwidth = 20,cellheight = 20,
                   color=rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)))
dev.copy2pdf(file="cytokines_heatmap.pdf",width=8,height=8)


####common genes between two groups
table(resgene_hfdko_hfdwt.sig$direction)
table(resgene_hfdko_ndko.sig$direction)

common_hfdko_hfdwt.and.hfdko_ndko.up=
   merge(
    resgene_hfdko_hfdwt.sig, 
    resgene_hfdko_ndko.sig,
    by = "row.names",  
    suffixes = c("_hfdwt", "_ndko")
  )
library(dplyr)
rownames(common_hfdko_hfdwt.and.hfdko_ndko.up) <- common_hfdko_hfdwt.and.hfdko_ndko.up$Row.names
common_hfdko_hfdwt.and.hfdko_ndko.up$Row.names <- NULL
common_hfdko_hfdwt.and.hfdko_ndko.direction <- common_hfdko_hfdwt.and.hfdko_ndko.up %>%
  filter(
    (direction_hfdwt == "up" & direction_ndko == "up") |
      (direction_hfdwt == "down" & direction_ndko == "down")
  ) %>%
  mutate(consensus_direction = ifelse(
    direction_hfdwt == "up", "Co-upregulated", "Co-downregulated"
  ))

write.csv(common_hfdko_hfdwt.and.hfdko_ndko.direction,file="common_hfdko_hfdwt.and.hfdko_ndko.direction.all.csv")
table(common_hfdko_hfdwt.and.hfdko_ndko.direction$consensus_direction)


###########volcano plot for two comparisons
#hfdko vs ndko 
resgene_hfdko_ndko.plot=as.data.frame(resgene_hfdko_ndko.sig)
resgene_hfdko_ndko.plot <- resgene_hfdko_ndko.plot %>%
  mutate(log10p_val_adj = -log10(replace(FDR, FDR == 0, min(FDR[FDR > 0], na.rm = TRUE))))
resgene_hfdko_ndko.plot$genesymbol=row.names(resgene_hfdko_ndko.plot)

table(resgene_hfdko_ndko.plot$direction)
pdf(file="resgene_hfdko_ndko.volcano.all.pdf",width=5,height=4.5)
ggplot(resgene_hfdko_ndko.plot, 
       aes(logFC, log10p_val_adj, color = direction)) +
  geom_point(size = 1) + 
  theme_bw() +
  ggtitle("hfdko vs ndko") +
  xlim(c(-10, 10)) +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dotted") +
  scale_color_manual(values = c("#0066CC","#FF3333")) +

  geom_point(data = resgene_hfdko_ndko.plot %>% 
               filter(genesymbol %in% c("Il1b","Nfkb2","Tnf","Il12a","Il36a","Cxcl9","Cxcl10","Ccl2","Ccl4","Ccl8","Ccl12")),
             size = 3, color = "black") +

  geom_point(data = resgene_hfdko_ndko.plot %>% 
               filter(genesymbol %in% c("Il1b","Nfkb2","Tnf","Il12a","Il36a","Cxcl9","Cxcl10","Ccl2","Ccl4","Ccl8","Ccl12")),
             size = 2) +

  geom_text_repel(data = resgene_hfdko_ndko.plot %>% 
                    filter(genesymbol %in% c("Il1b","Nfkb2","Tnf","Il12a","Il36a","Cxcl9","Cxcl10","Ccl2","Ccl4","Ccl8","Ccl12")),
                  aes(label = genesymbol),
                  size = 4,force = 1,
                  col = 'black',
                  box.padding = 0.5,
                  max.overlaps = Inf) 
dev.off()



#####hfdko hfdwt
resgene_hfdko_hfwt.plot=as.data.frame(resgene_hfdko_hfdwt.sig)
resgene_hfdko_hfwt.plot <- resgene_hfdko_hfwt.plot %>%
  mutate(log10p_val_adj = -log10(replace(FDR, FDR == 0, min(FDR[FDR > 0], na.rm = TRUE))))
resgene_hfdko_hfwt.plot$genesymbol=row.names(resgene_hfdko_hfwt.plot)
resgene_hfdko_hfwt.plot$direction="down"
resgene_hfdko_hfwt.plot$direction[resgene_hfdko_hfwt.plot$logFC>0]="up"

pdf(file="resgene_hfdko_hfdwt.volcano.all.pdf",width=5,height=4.5)
ggplot(resgene_hfdko_hfwt.plot, 
       aes(logFC, log10p_val_adj, color = direction)) +
  geom_point(size = 1) +  
  theme_bw() +
  ggtitle("hfdko vs hfdwt") +
  xlim(c(-10, 10)) +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dotted") +
  scale_color_manual(values = c("#0066CC","#FF3333")) +
  geom_point(data = resgene_hfdko_hfwt.plot %>% 
               filter(genesymbol %in% c("Il1b","Nfkb2","Tnf","Il12a","Il36a","Cxcl9","Cxcl10","Ccl2","Ccl4","Ccl8","Ccl12")),
             size = 3, color = "black") +
  geom_point(data = resgene_hfdko_hfwt.plot %>% 
               filter(genesymbol %in% c("Il1b","Nfkb2","Tnf","Il12a","Il36a","Cxcl9","Cxcl10","Ccl2","Ccl4","Ccl8","Ccl12")),
             size = 2) + 
  geom_text_repel(data = resgene_hfdko_hfwt.plot %>% 
                    filter(genesymbol %in% c("Il1b","Nfkb2","Tnf","Il12a","Il36a","Cxcl9","Cxcl10","Ccl2","Ccl4","Ccl8","Ccl12")),
                  aes(label = genesymbol),
                  size = 4,force = 1.5,
                  col = 'black',
                  box.padding = 0.5,  
                  max.overlaps = Inf) 
dev.off()

#####hfdkomet hfdko
resgene_hfdkomet_hfdfko.plot=as.data.frame(resgene_hfdkomet_hfdfko.sig)
resgene_hfdkomet_hfdfko.plot <- resgene_hfdkomet_hfdfko.plot %>%
  mutate(log10p_val_adj = -log10(replace(FDR, FDR == 0, min(FDR[FDR > 0], na.rm = TRUE))))
resgene_hfdkomet_hfdfko.plot$genesymbol=row.names(resgene_hfdkomet_hfdfko.plot)
resgene_hfdkomet_hfdfko.plot$direction="down"
resgene_hfdkomet_hfdfko.plot$direction[resgene_hfdkomet_hfdfko.plot$logFC>0]="up"

pdf(file="resgene_hfdkomet_hfdko.volcano.all.pdf",width=5,height=4.5)
ggplot(resgene_hfdkomet_hfdfko.plot, 
       aes(logFC, log10p_val_adj, color = direction)) +
  geom_point(size = 1) +  
  theme_bw() +
  ggtitle("hfdko.met vs hfdko") +
  xlim(c(-10, 10)) +
  geom_hline(yintercept = -log(0.05, 10), linetype = "dotted") +
  scale_color_manual(values = c("#0066CC","#FF3333")) +

  geom_point(data = resgene_hfdkomet_hfdfko.plot %>% 
               filter(genesymbol %in% c("Il1b","Nfkb2","Tnf","Il12a","Il36a","Cxcl9","Cxcl10","Ccl2","Ccl4","Ccl8","Ccl12")),
             size = 3, color = "black") +

  geom_point(data = resgene_hfdkomet_hfdfko.plot %>% 
               filter(genesymbol %in% c("Il1b","Nfkb2","Tnf","Il12a","Il36a","Cxcl9","Cxcl10","Ccl2","Ccl4","Ccl8","Ccl12")),
             size = 2) + 

  geom_text_repel(data = resgene_hfdkomet_hfdfko.plot %>% 
                    filter(genesymbol %in% c("Il1b","Nfkb2","Tnf","Il12a","Il36a","Cxcl9","Cxcl10","Ccl2","Ccl4","Ccl8","Ccl12")),
                  aes(label = genesymbol),
                  size = 4,force = 1.5,
                  col = 'black',
                  box.padding = 0.5,
                  max.overlaps = Inf) 
dev.off()



##Ets2 expression plot
library(ggplot2)
plotCounts(dds, gene = "Ets2", intgroup = "condition")
df_ets2=tmm_norm_counts_nolog["Ets2",,drop=F]
df_ets2=t(df_ets2)
df_ets2=as.data.frame(df_ets2)
df_ets2$group=c(rep("NDWT",3),rep("NDKO",3),rep("HFDWT",3),rep("HFDKO",3), rep("HFDKOMET",3))
df_ets2=df_ets2[4:12,]
df_ets2$group=as.factor(df_ets2$group)
levels(df_ets2$group)=c("HFDKO","HFDWT","NDKO")
ggplot(df_ets2, aes(group, Ets2, fill = group)) +
  stat_summary(fun = mean, geom = "bar") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(width = 0.2)+
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black")
    ) +
  ylim(0,300)+
  scale_fill_manual(values = c("#FF8080", "#FDCABF","#9B8EB6"))+
  ylab("CPM of Ets2")

dev.copy2pdf(file="Ets2.bulkRNA.barplot.pdf",width=4,height=5)





################################################################################################
##metascape_result plot
################################################################################################
library(ggplot2)
#extract data
tmp=read.csv("./hfdkomet_hfdfko_down_metascape/Enrichment_GO/_FINAL_GO.csv")
tmp=as.data.frame(tmp)
first_idx <- which(!duplicated(tmp$GROUP_ID))
first_idx
tmp=tmp[first_idx,]

#extract top 40 pathways
tmp$Description=gsub(" - Mus musculus \\(house mouse\\)","",tmp$Description)
tmp$Pathways=paste(tmp$GO,tmp$Description,sep=": ")
tmp$neglogP=-tmp$LogP
tmp$neglogQ=-tmp$Log.q.value.
tmp=head(tmp,40)


#bar plot
ggplot(data=tmp, aes(x=reorder(Pathways,neglogQ), y=neglogQ)) + #reorder 
  geom_bar(stat="identity", width=0.8,fill="#1F78B4") + coord_flip() + 
  xlab("") + ylab("")+
  labs(title = "")+
  theme_bw()+
  theme(axis.text=element_text(face = "bold", color="black",size = 12),
        axis.text.x = element_text(size = 12),  # X???̶?
        axis.text.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.copy2pdf(file="hfdkomet_hfdfko.metascape.down.kegg.go.top40.pdf",width=12,height=0.25*(nrow(tmp)+1))



######################################################################
#metascape_result plot for upregulated genes in both hfdko vs hfdwt and hfdko vs ndko
library(ggplot2)
tmp=readxl::read_xlsx("./common_hfdko_hfdwt.and.hfdko_ndko.direction.up/metascape_result.xlsx",sheet = 2)
tmp=as.data.frame(tmp)
tmp=tmp[grep("Summary",tmp[,1]),]
tmp$Description=gsub(" - Mus musculus \\(house mouse\\)","",tmp$Description)
tmp$Pathways=paste(tmp$Term,tmp$Description,sep=": ")
tmp$neglogP=-tmp$LogP
tmp$neglogQ=-tmp$`Log(q-value)`
tmp=head(tmp,50)
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
dev.copy2pdf(file="common_hfdko_hfdwt.and.hfdko_ndko.metascape.up.kegg.go.top20.pdf",width=9,height=0.3*(nrow(tmp)+1))


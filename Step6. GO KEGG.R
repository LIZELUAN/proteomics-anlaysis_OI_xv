library(tidyverse)
library(openxlsx)
library(org.Hs.eg.db)
library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
library(clusterProfiler)
memory.limit(size=100000000)

library(Rgraphviz)

data_up <- data.frame(Gene=DEG_up)
data_down <- data.frame(Gene=DEG_down)

# 转换基因
# 上调基因处理
keytypes(org.Hs.eg.db)
data_up$ID = mapIds(x=org.Hs.eg.db,
                    keys= data_up$Gene,
                    keytype="SYMBOL",
                    column="ENTREZID")

data_up <- na.omit(data_up)

# 下调基因处理
data_down$ID = mapIds(x=org.Hs.eg.db,
                    keys= data_down$Gene,
                    keytype="SYMBOL",
                    column="ENTREZID")

data_down <- na.omit(data_down)

# Create list for gene list with color in KEGG viewer
data_up_color <- data_up %>%
  select(ID) %>%
  mutate(color="pink")

data_down_color <- data_down %>%
  select(ID) %>%
  mutate(color="cyan")


gene_for_view <- rbind(data_up_color,data_down_color)


write.csv(gene_for_view, file="gene_for_view_uc.csv",row.names = F)



#############################################

#  上调基因GO富集
##############################################

ego_up_all <- enrichGO(gene          = data_up$ID,
                   OrgDb         = 'org.Hs.eg.db',
                   keyType       = 'ENTREZID',
                   ont           = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
head(ego_up_all,100)


p_up <- dotplot(ego_up_all,showCategory = 20,split='ONTOLOGY') + ggtitle("Top20 GO terms for each subclass for up-DEG in OI_WNT")+
  theme(axis.title = element_text(size=7),plot.title=element_text(size=7,face = "bold.italic"),axis.text.y = element_text(size=7),axis.text.x = element_text(size=7))+
  ggforce::facet_col(vars(ONTOLOGY), scales = 'free', space = 'free', 
                     strip.position='right') 

pdf('GO_results2/GO_up_top20.pdf',w=5,h=10)
p_up
dev.off()

###################################################
# 下调基因GO富集
####################################################
ego_down_all <- enrichGO(gene     = data_down$ID,
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'ENTREZID',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
head(ego_down_all)



p_down <- dotplot(ego_down_all,showCategory = 20,split='ONTOLOGY') + ggtitle("Top20 GO terms for each subclass for down-DEG in OI_WNT")+theme(axis.title = element_text(size=7),plot.title=element_text(size=7,face = "bold.italic"),axis.text.y = element_text(size=7),axis.text.x = element_text(size=7))+
  ggforce::facet_col(vars(ONTOLOGY), scales = 'free', space = 'free', 
                     strip.position='right') 


pdf('GO_results2/GO_all_down_10.pdf',w=5,h=11)
p_down
dev.off()

###############################################################
############################K E G G 分 析#########################
###############################################################
# 上调基因


# install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto' )


ekegg_up <- enrichKEGG(gene=data_up$ID,
                   organism = 'hsa',
                   keyType="ncbi-geneid",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  =  0.05,
                   pAdjustMethod = "BH")

head(ekegg_up)

p1 <- dotplot(ekegg_up, showCategory=100,title="KEGG for up-DEGs in OI_WNT") +
  theme(axis.title = element_text(size=15),
        plot.title=element_text(size=10, face = "bold.italic",hjust = 0.5),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=15))

pdf('GO_results2/KEGG_up.pdf',w=6,h=4.3)
p1
dev.off()

# 下调基因
ekegg_down <- enrichKEGG(gene=data_down$ID,
                    organism = 'hsa',
                    keyType="ncbi-geneid",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  =  0.05,
                    pAdjustMethod = "BH")
head(ekegg_down)

p2 <- dotplot(ekegg_down, showCategory=20,title="KEGG for down-DEGs in OI_WNT")+
  theme(axis.title = element_text(size=15),
        plot.title=element_text(size=10, face = "bold.italic",hjust = 0.5),
        axis.text.y = element_text(size=12),axis.text.x = element_text(size=15))

pdf('GO_results2/KEGG_down.pdf',w=5.5,h=4.3)
p2
dev.off()

##########################################################
# 研究KEGG重要通路
# KEGG up
browseKEGG(ekegg_up, 'hsa04216')

browseKEGG(ekegg_up, 'hsa00980')

browseKEGG(ekegg_up, 'hsa00480')

# KEGG down
browseKEGG(ekegg_down, 'hsa04510')

browseKEGG(ekegg_down, 'hsa04810')






# total KEGG关联网络图
options(ggrepel.max.overlaps = Inf)

KEGG_up <- pairwise_termsim(ekegg_up)
pdf('GO_results2/pathway network_up.pdf',w=7,h=7)
enrichplot::emapplot(KEGG_up,showCategory =20, 
                     color = "p.adjust", 
                     layout = "kk")
dev.off()

KEGG_down <- pairwise_termsim(ekegg_down)
pdf('GO_results2/pathway network_down.pdf',w=12,h=9)
enrichplot::emapplot(KEGG_down,showCategory =20, 
                     color = "p.adjust", 
                     layout = "kk")
dev.off()

# deg up
ekegg_up <- setReadable(ekegg_up, 'org.Hs.eg.db', 'ENTREZID')

pdf("GO_results2/KEGG_gene_pathway_up.pdf",w=8,h=5)
enrichplot::cnetplot(ekegg_up,circular=T,colorEdge = TRUE,showCategory = 100,node_lable="gene")
dev.off()

# deg down
ekegg_down <- setReadable(ekegg_down, 'org.Hs.eg.db', 'ENTREZID')

pdf("GO_results2/KEGG_gene_pathway_down.pdf",w=8,h=5)
enrichplot::cnetplot(ekegg_down,circular=T,colorEdge = TRUE,showCategory = 100,node_lable="gene")
dev.off()

ego_up_all <- setReadable(ego_up_all, 'org.Hs.eg.db', 'ENTREZID')
pdf("GO_results2/GO_MF_up.pdf",w=10,h=7)
enrichplot::cnetplot(ego_up_all,circular=T,colorEdge = TRUE,node_lable="gene")
dev.off()

ego_down_all <- setReadable(ego_down_all, 'org.Hs.eg.db', 'ENTREZID')
pdf("GO_results2/GO_MF_down.pdf",w=10,h=7)
enrichplot::cnetplot(ego_down_all,circular=T,colorEdge = TRUE,node_lable="gene")
dev.off()

save(file= "ibmx.Rda", ibmx_up, ibmx_down)




library(VennDiagram) 
x <- read.csv("DEG results2/HUB genes.csv") %>% as.list()
# Find hub genes
# Venn plot of downer DEG in UC
table <- get.venn.partitions(x)


str_c(table[1,7][[1]],collapse = ", ")











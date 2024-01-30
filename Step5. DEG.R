setwd("C:/Users/may/Desktop/project/diann-4vs4-20230529")
library(tidyverse)
library(DEP)
library(openxlsx)
# BiocManager::install("DEP")
# colnames(da)
da <- read_csv("txt/report.pg_matrix.csv") %>%
  mutate(across(starts_with("Control")|starts_with("OI"), ~log2(.x)))
da$`Genes` %>% duplicated() %>% any()

## if TRUE  :) #######################
# Find the duplicated genes
da %>% group_by(`Genes`) %>% summarise(frequency = n()) %>%
  arrange(desc(frequency)) %>% filter(frequency > 1)

# 确保基因名唯一
# 使用Gene.names中的注释作为主要名称,使用Protein.IDs中的注释作为没有基因名的蛋白质的ID
da <- make_unique(da,"Genes","Protein.Ids",delim = ";")

# 验证
da$name %>% duplicated() %>% any()
# [1] FALSE

da <- da %>%
  select(name, starts_with("Control"), starts_with("OI")) %>%
  select(-Control1) %>%
  select(-OI1) %>%
  rename(Gene="name")


da <- da %>%
  rowwise() %>%
  mutate(Num_valid_Ctrl = sum(!is.na(c_across(starts_with("Control"))))) %>%
  mutate(Num_valid_OI = sum(!is.na(c_across(starts_with("OI"))))) %>% ungroup() %>%
  filter(!(Num_valid_Ctrl+Num_valid_OI==0))



##################Output statistical result#############################
write.csv(deg,file="DEG results/DEG.csv",row.names = F)


##### volcano plot ##########
cut_off_pvalue=0.05
cut_off_logFC = 0.585

data_deg_cp$change = ifelse(data_deg_cp$pvalue < cut_off_pvalue & abs(data_deg_cp$`logfc(Exp/Ctrl)`) >= cut_off_logFC, 
                   ifelse(data_deg_cp$`logfc(Exp/Ctrl)`>= cut_off_logFC,'Up','Down'),'Stable')

dm<-data_deg_cp[data_deg_cp$pvalue<0.05 & abs(data_deg_cp$`logfc(Exp/Ctrl)`) >= 0.585,]



a<-theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=element_text(size=10,face = "bold",hjust = 0.5),
        axis.title = element_text(size=10),
        axis.text = element_text(size=8),
        #legend.position = "right",
        legend.text = element_text(size=8),
        legend.title=element_blank(),
        aspect.ratio = 1,
        legend.box.margin = margin(0,0,0,-0.1,unit = 'cm'),
        legend.justification = c(0,1.1),
        legend.key.width = unit(0.1,"cm")
  ) 
colnames(data_deg_cp)

library(ggrepel)
p<-ggplot(
  data_deg_cp, aes(x = `logfc(Exp/Ctrl)`, y = -log2(pvalue)))+  
  geom_point(aes(color=change),alpha=0.5, size=1) +
  
  #自定义颜色
  scale_color_manual(values = c("#ff4757","#546de5","#d2dae2"), breaks=c("Up","Down","Stable"))+
  
  # Add auxiliary line
  geom_vline(xintercept=c(-0.585,0.585),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(cut_off_pvalue),lty=4,col="black",lwd=0.5) +
  
  
  # Add axis
  labs(x="log2 FoldChange(OI/Control)",
       y="-log2(pval)")+
  
  # Add labels of some position
  geom_text_repel(data=dm,aes(x=`logfc(Exp/Ctrl)`,y= -log2(pvalue),label=Gene,color=change),
                  max.overlaps = 100,    # 最大覆盖率，调大该值则不被覆盖，反之。
                  size = 1,    # 字体大小
                  force=0.05,force_pull=2,
                  box.padding = unit(0.2, "lines"),  # 标记离box距离
                  point.padding = unit(0.03, "lines"), #标记离点距离
                  segment.color = "black",   # 标记线条的颜色
                  min.segment.length=0.2,
                  segment.size = 0.1,
                  show.legend = FALSE)+
  # xlim(-20,20)+
  # ggtitle("Frontal Cortex")+
  a

pdf("DEG results/volcano_plot.pdf",w=5,h=4)
p
dev.off()

colnames(data_deg_cp)
DEG <- data_deg_cp %>%
  filter(change!="Stable") %>%
  select(Gene, pvalue, q.value, change)




## ========================= Start from here ============================ #
up <- (DEG %>% filter(change=="Up")) %>% select(Gene)
down <- (DEG %>% filter(change=="Down")) %>% select(Gene)

# write.csv(rbind(up,up1),"GO enrichment/DEGlist_up.csv",row.names = F)
# write.csv(rbind(down,down1),"GO enrichment/DEGlist_down.csv",row.names = F)


DEG_comp <- list()
DEG_comp[["Up"]] <- up
DEG_comp[["Down"]] <- down

### Differentiate them into categories 
library(openxlsx)

cate <- read.xlsx("category_resource/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx"
) %>% select(Category, Gene.Symbol)
cate_resource <- list()
for (i in unique(cate$Category)) {
  A <- cate %>%
    filter(Category==i)
  cate_resource[[i]] <- A$Gene.Symbol
}

######### 创造一个列表
x <- list()
for (i in names(cate_resource)) {
  for (j in c("Up","Down")) {
    x[[j]][[i]] <- (DEG_comp[[j]] %>%
                      filter(Gene %in% cate_resource[[i]]))$Gene
  }
}
for (j in c("Up","Down")) {
    x[[j]][["Non-matrisome"]] <- (DEG_comp[[j]] %>%
                      filter(!Gene %in% unlist(cate_resource, use.names = FALSE)))$Gene
}


table_dep <- matrix(data=NA,nrow=7,ncol=2)
row.names(table_dep) <- names(x[["Up"]])
colnames(table_dep) <- c("Higher_in_OI_WNT","Higher_in_Control")

for (i in names(x[["Up"]])) {
  table_dep[i,"Higher_in_OI_WNT"] <- str_c(unlist(x[["Up"]][[i]]),collapse = ", ")
  table_dep[i,"Higher_in_Control"] <- str_c(unlist(x[["Down"]][[i]]),collapse = ", ")
  }

write.csv(table_dep,file="DEG results2/DEP category.csv")

################################################################################
# Heatmap based on DEGs ########################################################
da_deg <- da %>% filter(name %in% DEG$Gene)

da_heat <- da_deg %>% 
  select(starts_with("Control") | starts_with("OI"))

rownames(da_heat) <- da_deg$name


da_heat <- data.matrix(da_heat)
library(gplots)
library(RColorBrewer)
library(devtools)

colnames(da_heat) <- sample_names

dist_no_na <- function(da_t) {
  edist <- dist(da_t)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}

pdf("DEG results2/Heatmap of DEP.pdf")
heatmap.2(da_heat,
          trace="none",
          symbreaks = T,
          na.color = "gray",
          scale="row",
          col=bluered,
          density.info = c("none"),
          # col=rev(colorRampPalette(brewer.pal(10,"RdBu"))(20)),
          # cexRow = 0.3,
          cexCol = 1,
          cexRow = 0.3,
          srtCol = 45,
          keysize = 1,
          distfun = dist_no_na,
          
          #offsetCol=-0.5
          
          margins = c(7,16)
          #  key.title=NA
          #key.xlab = NA
)

dev.off()



############# Exclusive Genes in one group

data_deg_ex <- da %>%
  filter(Num_valid_pro_Ctrl >= 2 & Num_valid_pro_OI ==0 | Num_valid_pro_OI >= 2 & Num_valid_pro_Ctrl ==0)
data_deg_ex$change <- ifelse(data_deg_ex$Num_valid_pro_OI >= 2 & data_deg_ex$Num_valid_pro_Ctrl ==0, "Up","Down")

DEG_ex <- data_deg_ex %>%
  select(name, change)
up1 <- (DEG_ex %>% filter(change=="Up")) %>% select(name)
down1 <- (DEG_ex %>% filter(change=="Down")) %>% select(name)

DEG_ex <- list()
DEG_ex[["Up"]] <- up1
DEG_ex[["Down"]] <- down1





######### 创造一个列表
l <- list()
for (i in names(cate_resource)) {
  for (j in c("Up","Down")) {
    l[[j]][[i]] <- (DEG_ex[[j]] %>%
                      filter(name %in% cate_resource[[i]]))$name
  }
}
for (j in c("Up","Down")) {
  l[[j]][["Non-matrisome"]] <- (DEG_ex[[j]] %>%
                                  filter(!name %in% unlist(cate_resource, use.names = FALSE)))$name
}


table_dep_ex <- matrix(data=NA,nrow=7,ncol=2)
row.names(table_dep_ex) <- names(x[["Up"]])
colnames(table_dep_ex) <- c("Higher_in_OI_WNT","Higher_in_Control")

for (i in names(l[["Up"]])) {
  table_dep_ex[i,"Higher_in_OI_WNT"] <- str_c(unlist(l[["Up"]][[i]]),collapse = ", ")
  table_dep_ex[i,"Higher_in_Control"] <- str_c(unlist(l[["Down"]][[i]]),collapse = ", ")
}

write.csv(table_dep_ex,file="DEG results2/DEP(exclusive) category.csv")






# 2 fold change DEP
data_deg_2fold <- data_deg %>%
  filter(Num_valid_Ctrl >= 2 & Num_valid_OI == 1 | Num_valid_Ctrl == 1 & Num_valid_OI >= 2) %>%
  filter(abs(`logfc(Exp/Ctrl)`)>0.585)
data_deg_2fold$change <- ifelse(data_deg_2fold$`logfc(Exp/Ctrl)`>0.585, "Up","Down")

DEG_2fold <- data_deg_2fold %>%
  select(Gene, change)
up <- (DEG_2fold %>% filter(change=="Up")) %>% select(Gene)
down <- (DEG_2fold %>% filter(change=="Down")) %>% select(Gene)

DEG_2fold <- list()
DEG_2fold[["Up"]] <- up
DEG_2fold[["Down"]] <- down

######### 创造一个列表
l <- list()
for (i in names(cate_resource)) {
  for (j in c("Up","Down")) {
    l[[j]][[i]] <- (DEG_ex[[j]] %>%
                      filter(Gene %in% cate_resource[[i]]))$Gene
  }
}
for (j in c("Up","Down")) {
  l[[j]][["Non-matrisome"]] <- (DEG_ex[[j]] %>%
                                  filter(!Gene %in% unlist(cate_resource, use.names = FALSE)))$Gene
}


table_dep_2fold <- matrix(data=NA,nrow=7,ncol=2)
row.names(table_dep_2fold) <- names(x[["Up"]])
colnames(table_dep_2fold) <- c("Higher_in_OI","Higher_in_Control")

for (i in names(l[["Up"]])) {
  table_dep_2fold[i,"Higher_in_OI"] <- str_c(unlist(l[["Up"]][[i]]),collapse = ",")
  table_dep_2fold[i,"Higher_in_Control"] <- str_c(unlist(l[["Down"]][[i]]),collapse = ",")
}

write.csv(table_dep_2fold,file="DEG results/DEP(2 fold) category.csv")

################################################################################3
# Heatmap based on different DEGs ########################################################
# Exclusive DEGs
da_deg_ex <- da %>% filter(name %in% unlist(DEG_ex, use.names = FALSE))

da_heat <- da_deg_ex %>% 
  select(starts_with("Control") | starts_with("OI"))
# da_heat[is.na(da_heat)] = 0
rownames(da_heat) <- da_deg_ex$name

da_heat <- data.matrix(da_heat)
colnames(da_heat) <- sample_names

library(gplots)
library(RColorBrewer)
library(devtools)


dist_no_na <- function(da_t) {
  edist <- dist(da_t)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}

pdf("DEG results2/Heatmap of exclusive DEP.pdf")
heatmap.2(da_heat,
          trace="none",
          symbreaks = T,
          na.color = "gray",
          scale="row",
          col=bluered,
          density.info = c("none"),
          # col=rev(colorRampPalette(brewer.pal(10,"RdBu"))(20)),
          # cexRow = 0.3,
          cexCol = 1,
          cexRow = 0.3,
          srtCol = 45,
          keysize = 1,
          distfun = dist_no_na,
          
          #offsetCol=-0.5
          
          margins = c(7,16)
          #  key.title=NA
          #key.xlab = NA
)

dev.off()


# Heatmap based on different DEGs ########################################################
# 2 fold change DEGs
da_deg_2fold <- da %>% filter(Gene %in% unlist(DEG_2fold, use.names = FALSE))
da_heat <- da_deg_2fold %>% 
  select(starts_with("Control") | starts_with("OI"))
rownames(da_heat) <- da_deg_2fold$Gene
da_heat <- data.matrix(da_heat)

library(gplots)
library(RColorBrewer)
library(devtools)


dist_no_na <- function(da_t) {
  edist <- dist(da_t)
  edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
  return(edist)
}

pdf("DEG results/Heatmap of 2fold DEP.pdf")
heatmap.2(da_heat,
          trace="none",
          symbreaks = T,
          na.color = "gray",
          scale="row",
          col=bluered,
          density.info = c("none"),
          # col=rev(colorRampPalette(brewer.pal(10,"RdBu"))(20)),
          # cexRow = 0.3,
          cexCol = 1,
          cexRow = 0.3,
          srtCol = 45,
          keysize = 1,
          distfun = dist_no_na,
          
          #offsetCol=-0.5
          
          margins = c(7,16)
          #  key.title=NA
          #key.xlab = NA
)

dev.off()





# Heatmap based on all found DEGs ########################################################
da_deg_tot <- da %>% filter(name %in% c(unlist(DEG_comp,use.names = F),unlist(DEG_ex,use.names = F)))
da_heat <- da_deg_tot %>% 
  select(starts_with("Control") | starts_with("OI"))
rownames(da_heat) <- da_deg_tot$name
da_heat <- data.matrix(da_heat)
colnames(da_heat) <- sample_names

pdf("DEG results2/Heatmap of total DEP.pdf")
heatmap.2(da_heat,
          trace="none",
          symbreaks = T,
          na.color = "gray",
          scale="row",
          col=bluered,
          density.info = c("none"),
          # col=rev(colorRampPalette(brewer.pal(10,"RdBu"))(20)),
          # cexRow = 0.3,
          cexCol = 1,
          cexRow = 0.3,
          srtCol = 45,
          keysize = 1,
          distfun = dist_no_na,
          
          #offsetCol=-0.5
          
          margins = c(7,16)
          #  key.title=NA
          #key.xlab = NA
)

dev.off()


DEG_up <- c(unlist(DEG_comp[["Up"]],use.names = F),unlist(DEG_ex[["Up"]],use.names = F))

DEG_down <- c(unlist(DEG_comp[["Down"]],use.names = F),unlist(DEG_ex[["Down"]],use.names = F))

save(file="DEG.Rdata", DEG_up, DEG_down)

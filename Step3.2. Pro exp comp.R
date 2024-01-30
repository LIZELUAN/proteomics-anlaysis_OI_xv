setwd("C:/Users/oixra/Desktop/diann-4vs4-20230529")
library(tidyverse)
library(DEP)
library(data.table)
# BiocManager::install("DEP")
# colnames(da)


# 初步处理数据
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



################
# Protein exp comparison V2

da <- da %>% select(name,starts_with("Control"),starts_with("OI")) %>%
  rename(Gene=name)

sample_name <- c("Control1","Control2","Control3","Control4","OI_IFITM5","OI_WNT1","OI_WNT2","OI_WNT3")

colnames(da) <- c("Gene",sample_name)




da_long <-melt(da,
               id.vars = c('Gene'), 
               measure.va = c("Control1","Control2","Control3","Control4","OI_IFITM5","OI_WNT1","OI_WNT2","OI_WNT3"),
               variable.name='Sample',
               value.name='Exp')

b <- theme_classic() +
  theme(
    legend.position = "None",
    axis.title.x=element_blank()
  )

library(ggpubr)  
library(ggrepel)
library(ggsci)

p <- list()
for (i in c("CTNNB1","DKK1","CAPRIN1","BMP1","ACTA1","GAPDH","SOST","COL1A1","COL1A2","DMP1","RUNX2","SP7","SERPINH1","SERPINF1","CRTAP","FKBP10","BGLAP","SPP1","SPARC","P3H1","ACAN","CTNNB1","SATB2","ATF4","ALPL")) {
  da_i <- da_long %>%
    filter(Gene==i) 
  da_i <- data.frame(da_i,Grp=c("Control","Control","Control","Control","OI","OI","OI"))
  
  p[[i]] <- da_i %>%
    ggplot(aes(x = Grp, y = Exp)) +
    geom_violin(trim=T,aes(fill=Grp,color=Grp),alpha=0.2) +
    stat_compare_means(method = "wilcox.test", paired = FALSE,vjust=10) +
    theme_bw() + b +
    ggtitle(i) +
    geom_boxplot(aes(color=Grp),alpha=0.2,width=0.05, outlier.alpha = 0.01) +
    geom_jitter(shape=16, aes(color=Sample),position=position_jitter(width = 0.2, seed = 1)) +
    geom_text_repel(aes(label=Sample,color=Sample),position=position_jitter(width = 0.2, height = 0, seed = 1),size=2.5,force=0.02) +
    scale_color_aaas()
}


pdf("results/GENES_V2.2.pdf",w=9,h=9)
ggarrange(plotlist=p,ncol=4,nrow = 4)
dev.off()



gene_list <- c("CTNNB1","DKK1","CAPRIN1","BMP1","ACTA1","GAPDH","SOST","COL1A1","COL1A2","DMP1","RUNX2","SP7","SERPINH1","SERPINF1","CRTAP","FKBP10","BGLAP","SPP1","SPARC","P3H1","ACAN","CTNNB1","SATB2","ATF4")


################################ Correlation analysis #################################3
a<-theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=element_text(size=10,face = "bold.italic",hjust = 0.5),
        axis.title = element_text(size=10),
        axis.text = element_text(size=8),
        legend.text = element_text(size=10),
        legend.title=element_blank(),
        legend.position = "top",aspect.ratio = 1.2
  ) 

da_CTNNB1 <- da_long %>%
  filter(Gene=="CTNNB1") 

da_Runx2 <- da_long %>% 
  filter(Gene=="RUNX2") 

da_CRTAP <- da_long %>% 
  filter(Gene=="CRTAP") 

da_BMP1 <- da_long %>% 
  filter(Gene=="BMP1") 

da_SOST <- da_long %>%
  filter(Gene=="SOST") 


da_comp <- data.frame(Sample=da_CTNNB1$Sample,CTNNB1=da_CTNNB1$Exp,RUNX2=da_Runx2$Exp,CRTAP=da_CRTAP$Exp,BMP1=da_BMP1$Exp, SOST=da_SOST$Exp, Grp=c("Control","Control","Control","Control","OI_IFITM5","OI_WNT","OI_WNT","OI_WNT"))

# da_comp[is.na(da_comp)]<- 0


d <- list()


pdf("results/Correlation.pdf",w=4,h=3)
ggplot(data=da_comp, 
       aes(x=RUNX2, y=BMP1))+
  geom_point(aes(color=Grp,shape=Grp),size=3,alpha=0.6)+
  # xlab("Proportion of neurons")+ylab("Overall editing level")+
  geom_smooth(method="lm",formula=y~x)+
  # scale_color_discrete(labels=c("ALS","Ctrl","OND"))+
  # ggtitle("Correlation between CTNNB1 and GAPDH")+
  geom_text_repel(aes(label=Sample),size=2) +
  stat_cor(method = "pearson",size=3,vjust=  0.18,show.legend = FALSE)+
  a

dev.off()

d[[2]] <- ggplot(data=da_comp, 
       aes(x=CTNNB1, y=CAPRIN1))+
  geom_point(aes(color=Grp,shape=Grp),size=3,alpha=0.6)+
  # xlab("Proportion of neurons")+ylab("Overall editing level")+
  geom_smooth(method="lm",formula=y~x)+
  # scale_color_discrete(labels=c("ALS","Ctrl","OND"))+
  # ggtitle("Correlation between CTNNB1 and GAPDH")+
  geom_text_repel(aes(label=Sample),size=2) +
  stat_cor(method = "pearson",size=3,vjust=  0.18,show.legend = FALSE)+
  a


d[[3]] <- ggplot(data=da_comp, 
       aes(x=GAPDH, y=CAPRIN1))+
  geom_point(aes(color=Grp,shape=Grp),size=3,alpha=0.6)+
  # xlab("Proportion of neurons")+ylab("Overall editing level")+
  geom_smooth(method="lm",formula=y~x)+
  # scale_color_discrete(labels=c("ALS","Ctrl","OND"))+
  # ggtitle("Correlation between CTNNB1 and GAPDH")+
  geom_text_repel(aes(label=Sample),size=2) +
  stat_cor(method = "pearson",size=3,vjust=  0.18,show.legend = FALSE)+
  a

pdf("scatterplot.pdf",h=3)
ggarrange(plotlist=d,ncol=3,nrow = 1)
dev.off()

########################################################
#########################################################
# Overall exp level of different categories


library(openxlsx)

cate <- read.xlsx("category_resource/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx"
) %>% select(Category, Gene.Symbol)
cate_resource <- list()
for (i in unique(cate$Category)) {
  A <- cate %>%
    filter(Category==i)
  cate_resource[[i]] <- A$Gene.Symbol
}


sample_name <- c("Control1","Control2","Control3","Control4","OI_IFITM5","OI_WNT1","OI_WNT2","OI_WNT3")
Group <- c("Control","Control","Control","Control","OI","OI","OI") 
df <- da

######### 创造一个列表
x <- list()
for (i in names(cate_resource)) {
    xi <- df %>% filter(Gene %in% cate_resource[[i]]) %>%
      mutate(Category=i)
    
   xi1 <-melt(xi,
                  id.vars = c('Category'), 
                             measure.va = c("Control1","Control2","Control3","Control4"),
                             variable.name='Sample',
                             value.name='Log2_LFQ')
   xi1 <- data.frame(xi1,Grp="Control")
   xi2 <-melt(xi,
              id.vars = c('Category'), 
              measure.va = c("OI_IFITM5","OI_WNT1","OI_WNT2","OI_WNT3"),
              variable.name='Sample',
              value.name='Log2_LFQ')
   xi2 <- data.frame(xi2,Grp="OI_WNT")
   
   x[[i]] <- rbind(xi1, xi2)
   
   }
    
nm <- df %>% filter(!Gene %in% unlist(cate_resource, use.names = FALSE)) %>% 
      mutate(Category="Non-matrisome")

nm1 <-melt(nm,
           id.vars = c('Category'), 
           measure.va = c("Control1","Control2","Control3","Control4"),
           variable.name='Sample',
           value.name='Log2_LFQ')
nm1 <- data.frame(nm1,Grp="Control")
nm2 <-melt(nm,
           id.vars = c('Category'), 
           measure.va = c("OI_IFITM5","OI_WNT1","OI_WNT2","OI_WNT3"),
           variable.name='Sample',
           value.name='Log2_LFQ')
nm2 <- data.frame(nm2,Grp="OI_WNT")

x[["Non-matrisome"]] <- rbind(nm1, nm2)

# Core matrisome
cm <- df %>% filter(!Gene %in% c(cate_resource[["ECM Glycoproteins"]],cate_resource[["ECM Glycoproteins"]],cate_resource[["ECM Glycoproteins"]])) %>% 
  mutate(Category="Core-matrisome")

cm1 <-melt(cm,
           id.vars = c('Category'), 
           measure.va = c("Control1","Control2","Control3","Control4"),
           variable.name='Sample',
           value.name='Log2_LFQ')
cm1 <- data.frame(cm1,Grp="Control")
cm2 <-melt(cm,
           id.vars = c('Category'), 
           measure.va = c("OI_IFITM5","OI_WNT1","OI_WNT2","OI_WNT3"),
           variable.name='Sample',
           value.name='Log2_LFQ')
cm2 <- data.frame(cm2,Grp="OI_WNT")

x[["Core-matrisome"]] <- rbind(cm1, cm2)



# Non-core matrisome
ncm <- df %>% filter(!Gene %in% c(cate_resource[["ECM-affiliated Proteins"]],cate_resource[["ECM Regulators"]],cate_resource[["Secreted Factors"]])) %>% 
  mutate(Category="Core-matrisome")

ncm1 <-melt(ncm,
           id.vars = c('Category'), 
           measure.va = c("Control1","Control2","Control3","Control4"),
           variable.name='Sample',
           value.name='Log2_LFQ')
ncm1 <- data.frame(ncm1,Grp="Control")
ncm2 <-melt(ncm,
           id.vars = c('Category'), 
           measure.va = c("OI_IFITM5","OI_WNT1","OI_WNT2","OI_WNT3"),
           variable.name='Sample',
           value.name='Log2_LFQ')
ncm2 <- data.frame(ncm2,Grp="OI_WNT")

x[["Non-Core-matrisome"]] <- rbind(ncm1, ncm2)





b <- theme_classic2() +
  theme(legend.position = "None",
        axis.title.x=element_blank(),
        plot.title=element_text(size=12,face = "bold.italic",hjust = 0.5)
  )

library(ggpubr)  
library(ggrepel)
 p <- list()
 for (i in names(x)) {
 p[[i]] <- x[[i]] %>%
     filter(!(Sample=="OI_IFITM5")) %>%
     ggplot(aes(x = Sample, y =Log2_LFQ)) +
     geom_violin(trim=T,aes(fill=Grp,color=Grp),alpha=0.2) +
     # stat_compare_means(method = "wilcox.test", paired = FALSE) +
     theme_bw() + b +
     ggtitle(i) +
     geom_boxplot(aes(color=Grp),alpha=0.2,width=0.3) +
     # geom_jitter(shape=16, aes(color=Sample),position=position_jitter(width = 0.2, height=0, seed = 1)) +
     # geom_text_repel(aes(label=Sample,color=Sample),position=position_jitter(0.2, height=0, seed = 1),size=2.5,force=0.05)+
   scale_color_aaas()
 }

pdf("log2LFQ for each category.pdf",w=7,h=7) 
ggarrange(plotlist=p,ncol=3,nrow=3)
dev.off()






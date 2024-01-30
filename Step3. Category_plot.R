setwd("C:/Users/oixra/Desktop/OI_WNT")

library(tidyverse)
library(DEP)

# BiocManager::install("DEP")

#####################################################
## Preprocess of data with Gene name ################
#####################################################

df_ori <- read_csv("txt/report.pg_matrix.csv")
# %>% drop_na(`Gene names`)
df_ori$`Genes` %>% duplicated() %>% any()

## if TRUE 
# Find the duplicated genes
df_ori %>% group_by(Genes) %>% summarise(frequency = n()) %>%
  arrange(desc(frequency)) %>% filter(frequency > 1)

# 确保基因名唯一
# 使用Gene.names中的注释作为主要名称,使用Protein.IDs中的注释作为没有基因名的蛋白质的ID
df_ori <- make_unique(df_ori,"Genes","Protein.Ids",delim = ";") # 2941 25

# 验证
df_ori$name %>% duplicated() %>% any()

# 获得新列 Gene
df_ori <- df_ori %>% select(-Genes) %>% rename(Gene=name)


#去掉一些不想要的列,只保留重要的???
df <- df_ori %>%
  select(Gene|starts_with("Control")|starts_with("OI"))
sample_name <- c("Control1","Control2","Control3","Control4","OI_IFITM5","OI_WNT1","OI_WNT2","OI_WNT3")

colnames(df) <- c("Gene",sample_name)





#######################################################
##### Load category resource to classify proteins #############
#######################################################


# Create category resource

library(openxlsx)

cate <- read.xlsx("category_resource/Hs_Matrisome_Masterlist_Naba et al_2012.xlsx"
) %>% select(Category, Gene.Symbol)
cate_resource <- list()
for (i in unique(cate$Category)) {
  A <- cate %>%
    filter(Category==i)
  cate_resource[[i]] <- A$Gene.Symbol
}



x <- matrix(data=NA,nrow=6,ncol=8)
row.names(x) <- names(cate_resource)
colnames(x) <- sample_name
for (i in sample_name) { 
  for (j in names(cate_resource)) {
    x[j,i] <- df %>% drop_na(i) %>% 
      filter(Gene %in% cate_resource[[j]]) %>%
      nrow()
  }
}
x <- as.data.frame(x)
xsum <- x %>%
  summarise(across(everything(),~sum(.x)))

df_sum <- df %>%  
  summarize(across(!Gene, ~ sum(!is.na(.x))))

row.names(df_sum)<-"total"
x <- rbind(x,df_sum)


non_matrisome=df_sum-xsum
row.names(non_matrisome)<-"non-matrisome"
x <- rbind(x,non_matrisome)

core_matrisome <- data.frame(x) %>%
  filter(rownames(x) %in% c("ECM Glycoproteins", "Collagens", "Proteoglycans")) %>%  
  summarize(across(everything(), ~sum(.x))) 

row.names(core_matrisome)<-"core-matrisome"

x <- rbind(x,core_matrisome)


non_core_matrisome <- data.frame(x) %>%
  filter(rownames(x) %in% c("ECM-affiliated Proteins", "ECM Regulators", "Secreted Factors")) %>%  
  summarize(across(everything(), ~sum(.x)))

row.names(non_core_matrisome)<-"non_core-matrisome"

x <- rbind(x,non_core_matrisome)


x_test <- x %>% 
  mutate(Control1 = Control1/1265*100) %>% 
  mutate(Control2 = Control2/2618*100) %>% 
  mutate(Control3 = Control3/3806*100) %>% 
  mutate(Control4 = Control4/3796*100) %>% 
  mutate(OI_IFITM5 = OI_IFITM5/1913*100) %>% 
  mutate(OI_WNT1 = OI_WNT1/3732*100) %>% 
  mutate(OI_WNT2 = OI_WNT2/3557*100) %>% 
  mutate(OI_WNT3 = OI_WNT3/2383*100) 





library(data.table)
d <- setDT(x, keep.rownames = "Category")

d_long <-melt(d,
              id.vars = c('Category'), 
              measure.va = sample_name,
              variable.name='Group',
              value.name='Count')

plist <- c("non-matrisome","Collagens","Proteoglycans","ECM Glycoproteins","ECM-affiliated Proteins","ECM Regulators","Secreted Factors")
d_long$Category <- factor(d_long$Category,level=(plist))

a<-theme_classic()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.title=element_text(size=10,face = "bold.italic",hjust = 0.7),
        # axis.title.y = element_text(size=12),
        # axis.text = element_text(size=12),
        # axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        # legend.text = element_text(size=10),
        # legend.title=element_blank(),
        legend.justification = "top",
        legend.key.size = unit(0.4 , 'cm'), #change legend key size
        legend.key.height = unit(0.4, 'cm'), #change legend key height
        legend.key.width = unit(0.4, 'cm'), #change legend key width
        legend.title = element_text(size=9), #change legend title font size
        legend.text = element_text(size=9)
  ) 

pdf("results/number of proteins detected in categories.pdf",w=5,h=5)
# only real changes are geom bar and percent y axis
d_long %>% ggplot(aes(x=Group,y=Count,fill=Category)) + 
  geom_bar(stat="identity")+
  ylab("Number of proteins detected") +
  xlab("Samples")+
  scale_fill_manual("Category", values = c("non-matrisome" = "gray", "Collagens" = "lightgreen", "Proteoglycans" = "pink","ECM Glycoproteins"="yellow","ECM-affiliated Proteins"="blue","ECM Regulators"="lightblue","Secreted Factors"="purple")) +
  # scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, 50)) +
  ggtitle("Number of detected proteins in categories") +
  theme_bw() + a
# scale_fill_ucscgb()+scale_color_ucscgb()+

dev.off()


d_long_ECM <- d_long %>%
  filter(Category!="non-matrisome")


pdf("results/number of ECM proteins detected.pdf",w=5,h=5)
# only real changes are geom bar and percent y axis
d_long_ECM %>% ggplot(aes(x=Group,y=Count,fill=Category)) + 
  geom_bar(stat="identity")+
  ylab("Number of proteins detected") +
  xlab("Samples")+
  scale_fill_manual("Category", values = c("Collagens" = "lightgreen", "Proteoglycans" = "pink","ECM Glycoproteins"="yellow","ECM-affiliated Proteins"="blue","ECM Regulators"="lightblue","Secreted Factors"="purple")) +
  # scale_y_continuous(limits = c(0, 300), breaks = seq(0, 300, 50)) +
  ggtitle("Number of detected ECM proteins in categories") +
  # geom_text(aes(label=Count),position=position_stack(0.5)) +
  theme_bw() + a
# scale_fill_ucscgb()+scale_color_ucscgb()+

dev.off()

write.csv(d,file="results/protein category detected.csv",row.names = F)




d_long <- data.frame(d_long,Grp=rep(c("Control","OI"),each=40)) %>%
  rename(Sample=Group) %>%
  filter(!(Category=="total"))

b <- theme_classic2() +
theme(legend.position = "None",
         axis.title.x=element_blank(),
      plot.title=element_text(size=12,face = "bold.italic",hjust = 0.5)
   )
  
  
library(ggpubr)  
library(ggrepel)
library(ggsci)

pdf("results/violinplot for detected protein number in different categories.pdf")
p <- list()
for (i in unique(d_long$Category)) {
p[[i]] <- d_long %>% filter(Category==i) %>% filter(Sample!="Control1") %>%
ggplot(aes(x = Grp, y = Count)) +
  geom_violin(trim=T,aes(fill=Grp,color=Grp),alpha=0.2) +
  stat_compare_means(method = "wilcox.test", paired = FALSE) +
  theme_bw() + b +
  ylab("Percentage") +
  ggtitle(i) +
  geom_boxplot(aes(color=Grp),alpha=0.2,width=0.1) +
geom_jitter(shape=16, aes(color=Sample), position=position_jitter(width = 0.2, height = 0, seed = 1)) +
geom_text_repel(aes(label=Sample,color=Sample),position=position_jitter(0.2, height = 0, seed = 1),size=2.5,force=0.05) +
  scale_color_aaas()
}
pdf("results/number of proteins detected in categories_4vs4.pdf",w=10,h=6)
ggarrange(plotlist=p,ncol=4,nrow=2)
dev.off()

setwd("C:/Users/oixra/Desktop/project/diann-4vs4-20230529")

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
sample_name <- c("Control1","Control2","Control3","Control4","OI_IFTM5","OI_WNT1","OI_WNT2","OI_WNT3")

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

######### 创造一个列表
x <- list()
for (i in names(cate_resource)) {
  for (j in sample_name) {
    x[[i]][[j]] <- (df %>% drop_na(j) %>% 
                 filter(Gene %in% cate_resource[[i]]))$Gene
  }
}

for (j in sample_name) {
x[["Non-matrisome"]][[j]] <- (df %>% drop_na(j) %>% 
                                filter(!Gene %in% unlist(cate_resource, use.names = FALSE)))$Gene
}

Y <- list()
  for (j in sample_name) {
    Y[[j]] <- (df %>% drop_na(j))$Gene
  }

# ########################### Venn plot 旧版本 ####################################3
# # Import required package
# library("ggvenn")
# p_control <- list()
# 
# for (i in names(cate_resource)) {
#   #Make the plot
#   k <- ggvenn(x[[i]],c("Control1","Control2","Control3","Control4"),show_percentage = F,fill_alpha = 0,
#               set_name_size=3)
#   p_control[[i]] <- annotate_figure(k, bottom = text_grob(i,color = "black", face = "bold", size = 10))
#   }
# 
# ################ TOTAL #################33
# k <- ggvenn(Y,c("Control1","Control2","Control3","Control4"),show_percentage = F,fill_alpha = 0,
#             set_name_size=3)
# p_control[["Total"]] <- annotate_figure(k, bottom = text_grob("All proteins",color = "black", face = "bold", size = 10))
# 
# library(ggpubr)
# pdf("Category_Results/venn_plot_Control.pdf")
# ggarrange(plotlist = p_control,ncol = 2, nrow = 2)
# 
# dev.off()
# 
# ################################ Venn plot版本2 #########################################

# Import required package
library(VennDiagram)
p_ctrl<-list()
for (i in names(x)) {
  A <- x[[i]][["Control1"]]
  B <- x[[i]][["Control2"]]
  C <- x[[i]][["Control3"]]
  D <- x[[i]][["Control4"]]
  #Make the plot
  p_ctrl[[i]] <- venn.diagram(x = list(A,B,C,D),filename = NULL,
                         category.names = c("Ctrl1" , "Ctrl2","Ctrl3","Ctrl4"),
                         output=TRUE,main.fontface = "bold",
                         compression = "lzw", main=i,
                         lwd = 2)}
################ TOTAL #################
A <- Y[["Control1"]]
B <- Y[["Control2"]]
C <- Y[["Control3"]]
D <- Y[["Control4"]]
p_ctrl[["Total"]] <- venn.diagram(x = list(A,B,C,D),filename = NULL,
                                                      category.names = c("Ctrl1" , "Ctrl2","Ctrl3","Ctrl4"),
                                                      output=TRUE,main.fontface = "bold",
                                                      compression = "lzw", main="All proteins",
                                                      lwd = 2)
library(ggpubr)
pdf("Category_Results/venn_plot_Control.pdf")
ggarrange(plotlist = p_ctrl,ncol = 3, nrow = 3)

dev.off()


# ################################ Venn plot版本2 ——Exp Group #########################################

# Import required package
library(VennDiagram)
p_OI<-list()
for (i in names(x)) {
  A <- x[[i]][["OI_IFTM5"]]
  B <- x[[i]][["OI_WNT1"]]
  C <- x[[i]][["OI_WNT2"]]
  D <- x[[i]][["OI_WNT3"]]
  #Make the plot
  p_OI[[i]] <-  venn.diagram(x = list(A,B,C,D),filename = NULL,
                             category.names = c("OI_IFTM5", "OI_WNT1","OI_WNT2","OI_WNT3"),
                             cat.cex = 0.6,
                             main.fontfamily = "sans",
                             cat.fontfamily = "sans",
                             output=TRUE,main.fontface = "bold",
                             cat.fontface = "bold",
                             compression = "lzw", main=i,
                             main.pos = c(0.5,0.69),
                             lwd = 2,
                             margin=0.8) 
  } 
################ TOTAL #################
A <- Y[["OI_IFTM5"]]
B <- Y[["OI_WNT1"]]
C <- Y[["OI_WNT2"]]
D <- Y[["OI_WNT3"]]
p_OI[["Total"]] <- venn.diagram(x = list(A,B,C,D),filename = NULL,
                                  category.names = c("OI_IFTM5", "OI_WNT1","OI_WNT2","OI_WNT3"),
                                cat.cex = 0.6,
                                main.fontface = "bold",
                                  main.fontfamily = "sans",
                                  cat.fontfamily = "sans",
                                  output=TRUE,
                                cat.fontface = "bold",
                                  compression = "lzw", main="All proteins",
                                  lwd = 2,
                                main.pos = c(0.5,0.69),
                                margin=0.8)
library(ggpubr)
pdf("Category_Results/venn_plot_OI.pdf",w=7,h=7)
ggarrange(plotlist = p_OI,ncol = 3, nrow = 3)
dev.off()


######## 再创造一个列表
l <- list()
for (i in names(cate_resource)) {
    l[[i]][["Ctrl"]] <- (df %>% rowwise() %>% 
                      filter(!sum(is.na(c_across(starts_with("Control"))))==4) %>% ungroup() %>%
                      filter(Gene %in% cate_resource[[i]]))$Gene
    l[[i]][["OI_WNT"]] <- (df %>% rowwise() %>% 
                           filter(!sum(is.na(c_across(starts_with("OI_WNT"))))==3) %>% ungroup() %>%
                           filter(Gene %in% cate_resource[[i]]))$Gene
    
    l[[i]][["OI_IFTM5"]] <- (df %>% 
                         filter(!is.na(OI_IFTM5)) %>%
                         filter(Gene %in% cate_resource[[i]]))$Gene
  }

# Core matrisome









################## Non-matrisome
l[["Non-matrisome"]][["Ctrl"]] <- (df %>% rowwise() %>% 
                                    filter(!sum(is.na(c_across(starts_with("Control"))))==4) %>% ungroup() %>%
                                     filter(!Gene %in% unlist(cate_resource, use.names = FALSE)))$Gene

l[["Non-matrisome"]][["OI_WNT"]] <- (df %>% select(-OI_IFTM5) %>% rowwise() %>% 
                                    filter(!sum(is.na(c_across(starts_with("OI_WNT"))))==3) %>% ungroup() %>%
                                     filter(!Gene %in% unlist(cate_resource, use.names = FALSE)))$Gene

l[["Non-matrisome"]][["OI_IFTM5"]] <- (df %>%
                                      filter(!is.na(OI_IFTM5)) %>% 
                                       filter(!Gene %in% unlist(cate_resource, use.names = FALSE)))$Gene

###  For all proteins
l[["All proteins"]][["Ctrl"]] <- (df %>% rowwise() %>% 
                       filter(!sum(is.na(c_across(starts_with("Control"))))==4)
                       )$Gene
l[["All proteins"]][["OI_WNT"]] <- (df %>% rowwise() %>% 
                     filter(!sum(is.na(c_across(starts_with("OI_WNT"))))==3)
                    )$Gene
l[["All proteins"]][["OI_IFTM5"]] <- (df %>%
                                  filter(!is.na(OI_IFTM5)))$Gene


p<-list()
t <- list()
for (i in names(l)) {
  A <- l[[i]][["Ctrl"]]
  B <- l[[i]][["OI_WNT"]]
  C <- l[[i]][["OI_IFITM5"]]
  #Make the plot
  p[[i]] <- venn.diagram(x = list(A,B),filename = NULL,
                         category.names = c("Ctrl", "OI_WNT"),
                         output=TRUE,main.fontface = "bold",
                         compression = "lzw", main=i,col=c("#440154ff", '#21908dff'),
                         lwd = 1.5,scaled = FALSE, main.pos = c(0.5,0.605),margin=1.85,
                         main.fontfamily = "sans",
                         cat.fontfamily = "sans",
                         cat.fontface = "bold")
  
  
  
  t[[i]] <- get.venn.partitions(x = list(A,B), force.unique = TRUE, keep.elements = TRUE,
                           hierarchical = FALSE
                           )
}






pdf("Category_Results/venn_plot_OIvsCtrl2.pdf")
ggarrange(plotlist = p,ncol = 2, nrow = 4)

dev.off()
names(t)
A <- t[["ECM Glycoproteins"]]
B <- t[["Collagens"]]
C <- t[["Proteoglycans"]]
D <- t[["ECM-affiliated Proteins"]]
    E <- t[["ECM Regulators"]]
    f <- t[["Secreted Factors"]]


f[2,4]


f[3,4]
#

















############ 之前计数
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

non_matrisome=df_sum-xsum
row.names(non_matrisome)<-"non-matrisome"
x <- rbind(x,non_matrisome)

library(data.table)
d <- setDT(x, keep.rownames = "Category")

d_long <-melt(d,
              id.vars = c('Category'), 
              measure.va = c("Control1","Control2","Control3","Control4","OI1","OI2","OI3","OI4"),
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
        legend.key.size = unit(0.4
                               , 'cm'), #change legend key size
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














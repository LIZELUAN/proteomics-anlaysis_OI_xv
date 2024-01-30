setwd("C:/Users/oixra/Desktop/diann-4vs4-20230529")
library(tidyverse)
library(DEP)
library(openxlsx)
# BiocManager::install("DEP")
# colnames(da)
da <- read_csv("txt/report.pg_matrix.csv") %>%
  mutate(across(starts_with("Control")|starts_with("OI"), ~log2(.x))) %>%
  select(c(-Control1,-OI1))
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

sample_names = c("Control2","Control3","Control4","OI_WNT1","OI_WNT2","OI_WNT3")

# =========== NA statistics ====================
lfq1 <- da %>%
  select(starts_with("Control") | starts_with("OI"))
val1 <- sum(!is.na(lfq1))
na1 <- sum(is.na(lfq1))


# ========== Data Winsorisation ==========
da <- da %>%
  rowwise() %>%
  mutate(Num_valid_pro_Ctrl = sum(!is.na(c_across(starts_with("Control"))))) %>%
  mutate(Num_valid_pro_OI = sum(!is.na(c_across(starts_with("OI"))))) %>%
  mutate(sum_of_valid_pro = sum(!is.na(c_across(starts_with("Control") | starts_with("OI"))))) %>%
  arrange(-sum_of_valid_pro) %>% ungroup() %>%
  filter(!(sum_of_valid_pro==0))



da <- da %>% 
  rowwise() %>% 
  mutate(mean_exp=mean(c_across(starts_with("Control") | starts_with("OI")),na.rm=T)) %>% ungroup()

cut_point_top <- quantile(da$mean_exp, 0.95)
cut_point_bottom <- quantile(da$mean_exp, 0.05)



# =====================过滤缺失值=====================
da_filtered <- da %>%
  filter(Num_valid_pro_Ctrl>=2, Num_valid_pro_OI>=2, mean_exp < cut_point_top, mean_exp > cut_point_bottom)
# Select protein subset that capture the detected values#

# =========== NA statistics ====================
lfq2 <- da_filtered %>%
  select(starts_with("Control") | starts_with("OI"))
val2 <- sum(!is.na(lfq2))
na2 <- sum(is.na(lfq2))
valid_percent <- val2/val1*100
na_percent <- na2/na1*100


# =========== Da imputation by MICE ====================
library(mice)
lfq <- da_filtered %>%
  select(starts_with("Control") | starts_with("OI"))

# colnames(lfq) <- sample_names

imputed <- mice(lfq, m=5, maxit=50, meth = "pmm", seed=500)  ####### 生成五个数据集 
lfq <- complete(imputed,1)
# 更新数据
da_filtered[,colnames(lfq)]<-lfq

# =========== Check the normality of data ====================
densityplot(imputed)









# =========Data scaling & PCA ============
data_normalized <- da_filtered %>%
  filter(name %in% DEG$Gene) %>%
scale(center=T, scale=T) %>% select(starts_with("Control") | starts_with("OI"))
# da[,colnames(data_normalized)]<-data_normalized
# 

colnames(data_normalized) <- sample_names

# ============ PCA analysis ================
# library("FactoMineR")
# res.pca <- PCA(data_normalized, graph = TRUE)
pca <- prcomp(t(data_normalized), center=F, scale=F)
# plot(pca$x, main="after PCA")
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

library(ggplot2)
library(ggrepel)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("Control","Control","Control","OI","OI","OI")) 

pdf("step4_results/PCA.pdf",w=5,h=5)
ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(color = Group),size=3) +
  geom_text_repel(aes(label=Sample, color=Group), show.legend = F,force=1,force_pull = 1) +
  scale_color_manual(values = c("orange","purple")) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_classic() + theme(legend.position="None") +
  ggtitle("PCA graph")  
  #stat_ellipse(aes(fill = Group), level = 0.95, geom = 'polygon', alpha=0.1, show.legend = FALSE) +
  #scale_fill_manual(values = c('orange', 'purple'))

dev.off()

# ===== PCA & SVM ========================

# SVM
library(e1071)
library(grid)
library(lattice)
library(gridExtra)
library(stringr)

# Read in the data
dat <- data.frame(pca$x[,c(1,2)],Group=c(-1,-1,-1,1,1,1))
# Create a logistic rgression model on the reduced data
svmfit = svm(Group~.,data=dat,kernel="polynomial",coef0=1)

# plot(svmfit,dat)
# Visualize the decision boundaries for training set
X1 = seq(range(pca$x[,1])[1]*2, range(pca$x[,1])[2], len=20)
X2 = seq(range(pca$x[,2])[1]*1.5, range(pca$x[,2])[2]*1.05, len=20)
xgrid = expand.grid(PC1=X1, PC2=X2)
ygrid = predict(svmfit, xgrid)

# ygrid <- predict(svmfit, xgrid, decision.values=T)
pdf("PCA_SVM_3vs3_deg.pdf",w=4.5,h=4.5)
plot(dat[,-3], cex=2,col= alpha(ifelse(dat$Group == -1, "lightblue","red"),0.5),pch = 19,
     main = 'SVM after PCA',
     xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) ,
     ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")))
contour(X1, X2, matrix(as.numeric(ygrid),20,20), add=TRUE,level=0)
contour(X1, X2, matrix(as.numeric(ygrid), 20,20), add=TRUE,level=0.5,lwd=2,lty=2,col="red")
contour(X1, X2, matrix(as.numeric(ygrid), 20,20), add=TRUE,level=-0.5,lwd=2,lty=2,col="blue")
text(dat[,-3],labels=row.names(dat),cex=0.3)
dev.off()

# =================== DEG 分析 ===========================================
data <- da_filtered %>% 
  select(starts_with("Control") | starts_with("OI")) %>%
  as.data.frame()

rownames(data) <- da_filtered$name



# limma differential expression analysis
library(limma)

group_factor <- factor(c("Control","Control","Control","OI","OI","OI"), levels=c("Control","OI"),ordered = F)

# Experimental group: OI   Control group: Control
design <-model.matrix(~group_factor)
fit<-lmFit(data,design)
fit<-eBayes(fit)
options(digits = 4)
deg<-topTable(fit,adjust='BH',number = Inf)

# Get the up, down, stable genes
cut_off_pvalue <- 0.05
cut_off_logFC <- 0.585

deg$change <- ifelse(deg$P.Value < cut_off_pvalue & abs(deg$logFC) >= cut_off_logFC, 
                     ifelse(deg$logFC>= cut_off_logFC,'Up','Down'),'Stable')

# Identify the up and down genes as DEGs




##################Output statistical result#############################
write.csv(deg,file="step4_results2/DEG2.csv",row.names = T)

# =====volcano plot ====================================================================================
cut_off_pvalue=0.05
cut_off_logFC = 0.585
deg$Gene <- row.names(deg) 
dm <- deg[deg$change!="Stable",]


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
colnames(deg)

p<-ggplot(
  deg, aes(x = logFC, y = -log2(P.Value)))+  
  geom_point(aes(color=change),alpha=0.5, size=1) +
  
  #自定义颜色
  scale_color_manual(values = c("#ff4757","#546de5","#d2dae2"), breaks=c("Up","Down","Stable"))+
  
  # Add auxiliary line
  geom_vline(xintercept=c(-0.585,0.585),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(cut_off_pvalue),lty=4,col="black",lwd=0.5) +
  
  
  # Add axis
  labs(x="log2 FoldChange(OI_WNT/Control)",
       y="-log2(pval)")+
  
  # Add labels of some position
  geom_text_repel(data=dm,aes(x=logFC,y= -log2(P.Value),label=Gene,color=change),
                  max.overlaps = 100,    # 最大覆盖率，调大该值则不被覆盖，反之。
                  size = 2,    # 字体大小
                  box.padding = unit(0.2, "lines"),  # 标记离box距离
                  point.padding = unit(0.03, "lines"), #标记离点距离
                  segment.color = "black",   # 标记线条的颜色
                  show.legend = FALSE)+
  # xlim(-20,20)+
  # ggtitle("Frontal Cortex")+
  a

pdf("step4_results2/volcano_plot.pdf",w=5,h=5)
p
dev.off()

DEG <- deg %>%
  filter(change!="Stable") %>%
  select(Gene, P.Value, adj.P.Val, change)


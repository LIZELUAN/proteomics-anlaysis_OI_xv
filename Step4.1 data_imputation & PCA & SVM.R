setwd("C:/Users/oixra/Desktop/diann-4vs4-20230529")
library(tidyverse)
library(DEP)
library(openxlsx)

da <- read_csv("txt/report.pg_matrix.csv") %>%
  mutate(across(starts_with("Control")|starts_with("OI"), ~log2(.x)))
  # select(c(-Control1,-OI1))
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

da <- da %>% select(name,starts_with("Control"),starts_with("OI")) %>%
  rename(Gene=name)

# sample_name <- c("Control1","Control2","Control3","Control4","OI_IFITM5","OI_WNT1","OI_WNT2","OI_WNT3")

# colnames(da) <- c("Gene",sample_name)



# =========== NA statistics ====================
lfq1 <- da %>%
  select(-Gene)
val1 <- sum(!is.na(lfq1))
na1 <- sum(is.na(lfq1))

# =====================过滤缺失值=====================
da <- da %>%
  rowwise() %>%
  mutate(Num_valid_Ctrl = sum(!is.na(c_across(starts_with("Control"))))) %>%
  mutate(Num_valid_OI = sum(!is.na(c_across(starts_with("OI"))))) %>% ungroup() %>%
  filter(!(Num_valid_Ctrl+Num_valid_OI==0))
# Select protein subset that capture the detected values#

# ========== Data Winsorisation ==========
da <- da %>% 
  rowwise() %>% 
  mutate(mean_exp=mean(c_across(starts_with("Control") | starts_with("OI")),na.rm=T)) %>% ungroup()

cut_point_top <- quantile(da$mean_exp, 0.95)
cut_point_bottom <- quantile(da$mean_exp, 0.05)

da_filtered <- da %>% 
  filter(mean_exp < cut_point_top & mean_exp > cut_point_bottom) %>%
  filter(Num_valid_Ctrl>=4 & Num_valid_OI>=4) # # Select protein subset that capture the detected values#



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

imputed <- mice(lfq, m=5, maxit=50, meth = "pmm", seed=500)  ####### 生成五个数据集 
lfq <- complete(imputed,1)
# 更新数据
da_filtered[,colnames(lfq)]<-lfq




# =========Data scaling & PCA ============
data_normalized <- da_filtered %>%
  select(starts_with("Control") | starts_with("OI"))
# da[,colnames(data_normalized)]<-data_normalized
# 

# ============ PCA analysis ================
pca <- prcomp(t(data_normalized))
# plot(pca$x, main="after PCA")
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

library(ggplot2)
library(ggrepel)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("Control","Control","Control","Control","OI","OI","OI","OI")) 

pdf("PCA.pdf",w=3,h=3)
ggplot(data=pca.data, aes(x=X, y=Y)) +
  geom_point(aes(color = Group),size=3) +
  geom_text_repel(aes(label=Sample, color=Group), show.legend = F) +
  scale_color_manual(values = c("orange","purple")) +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
  theme_classic() + theme(legend.position="None") +
  ggtitle("PCA graph") +  
  stat_ellipse(level = 0.95, geom = 'polygon', alpha=0.1, show.legend = FALSE)
  #scale_fill_manual(values = c('orange', 'purple'))

dev.off()


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
pdf("PCA_SVM_3vs3.pdf",w=4.5,h=4.5)
plot(dat[,-3], cex=2,col= alpha(ifelse(dat$Group == -1, "lightblue","red"),0.5),pch = 19,
     main = 'SVM after PCA',
     xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) ,
       ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")))
contour(X1, X2, matrix(as.numeric(ygrid),20,20), add=TRUE,level=0)
contour(X1, X2, matrix(as.numeric(ygrid), 20,20), add=TRUE,level=0.5,lwd=2,lty=2,col="red")
contour(X1, X2, matrix(as.numeric(ygrid), 20,20), add=TRUE,level=-0.5,lwd=2,lty=2,col="blue")
text(dat[,-3],labels=row.names(dat),cex=0.3)
dev.off()


########## SVM without OI1, Control1 #######################################################3

# =========Data scaling & PCA ============
data_normalized <- da_filtered %>%
  select(starts_with("Control") | starts_with("OI"))
# da[,colnames(data_normalized)]<-data_normalized
# 

# ============ PCA analysis ================
pca <- prcomp(t(data_normalized %>% select(-OI1) %>% select(-Control1)))
# plot(pca$x, main="after PCA")
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

library(ggplot2)
library(ggrepel)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("Control","Control","Control","OI","OI","OI")) 





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
pdf("PCA_SVM_3vs3.pdf",w=5,h=5)
plot(dat[,-3], cex=2,col= ifelse(dat$Group == -1, "lightblue","red"),pch = 19,
     main = 'SVM after Principal Component Analysis',
     xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) ,
     ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")))
contour(X1, X2, matrix(as.numeric(ygrid),20,20), add=TRUE,level=0)
contour(X1, X2, matrix(as.numeric(ygrid), 20,20), add=TRUE,level=0.5,lwd=2,lty=2,col="red")
contour(X1, X2, matrix(as.numeric(ygrid), 20,20), add=TRUE,level=-0.5,lwd=2,lty=2,col="blue")
text(dat[,-3],labels=row.names(dat))
dev.off()



# ############ Select the Control Group ##################################################################3
# setwd("C:/Users/may/Desktop/project/diann-4vs4-20230529")
# library(tidyverse)
# library(DEP)
# library(openxlsx)
# da <- read_csv("txt/report.pg_matrix.csv") %>%
#   mutate(across(starts_with("Control")|starts_with("OI"), ~log2(.x)))
# da$`Genes` %>% duplicated() %>% any()
# 
# ## if TRUE  :) #######################
# # Find the duplicated genes
# da %>% group_by(`Genes`) %>% summarise(frequency = n()) %>%
#   arrange(desc(frequency)) %>% filter(frequency > 1)
# 
# # 确保基因名唯一
# # 使用Gene.names中的注释作为主要名称,使用Protein.IDs中的注释作为没有基因名的蛋白质的ID
# da <- make_unique(da,"Genes","Protein.Ids",delim = ";")
# 
# # 验证
# da$name %>% duplicated() %>% any()
# # [1] FALSE
# 
# 
# #选择对照组
# da_con <- da %>% 
#   select("name",starts_with("Control")) %>% drop_na()
# 
# 
# # =====================过滤缺失值=====================
# da_con <- da_con %>%
#   rowwise() %>%
#   mutate(Num_valid_pro_Ctrl = sum(!is.na(c_across(starts_with("Control"))))) 
# 
# 
# # ========== Data Winsorisation ==========
# da_con <- da_con %>% 
#   rowwise() %>% 
#   mutate(mean_exp=mean(c_across(starts_with("Control")),na.rm=T)) %>% ungroup()
# 
# cut_point_top <- quantile(da_con$mean_exp, 0.95)
# cut_point_bottom <- quantile(da_con$mean_exp, 0.05)
# 
# da_con_filtered <- da_con %>% 
#   filter(mean_exp < cut_point_top & mean_exp > cut_point_bottom) %>%
#   filter(Num_valid_pro_Ctrl>=2)# # Select protein subset that capture the detected values#
# 
# # =========== da imputation by MICE ====================
# library(mice)
# lfq <- da_con_filtered %>%
#   select(starts_with("Control"))
# 
# imputed <- mice(lfq, m=5, maxit=50, meth = "pmm", seed=500)  ####### 生成五个数据集 
# lfq <- complete(imputed,1)
# # 更新数据
# da_con_filtered[,colnames(lfq)]<-lfq
# 
# 
# # =========data scaling & PCA ============
# da_con_normalized <- da_con_filtered %>%
#   select(starts_with("Control"))
# # da_con[,colnames(da_conta_normalized)]<-da_conta_normalized
# # 
# 
# # ============ PCA analysis ================
# pca <- prcomp(t(da_con_normalized))
# # plot(pca$x, main="after PCA")
# pca.var <- pca$sdev^2
# pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
# 
# 
# library(ggplot2)
# library(ggrepel)
# pca.data <- data.frame(Sample=rownames(pca$x),
#                        X=pca$x[,1],
#                        Y=pca$x[,2],
#                        Group=c("Control","Control","Control","Control")) 
# 
# 
# 
# pdf("PCA-control.pdf",w=3,h=3)
# ggplot(data=pca.data, aes(x=X, y=Y)) +
#   geom_point(aes(color = Group),size=3) +
#   geom_text_repel(aes(label=Sample, color=Group), show.legend = F) +
#   scale_color_manual(values = c("orange","purple")) +
#   xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
#   ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
#   theme_classic() + theme(legend.position="None") +
#   ggtitle("PCA graph")  
# #stat_ellipse(aes(fill = Group), level = 0.95, geom = 'polygon', alpha=0.1, show.legend = FALSE) +
# #scale_fill_manual(values = c('orange', 'purple'))
# 
# dev.off()
# 
# #############################################################
# ##############################################################
# # PC1 PC3
# 
# pca.data <- data.frame(Sample=rownames(pca$x),
#                        X=pca$x[,1],
#                        Y=pca$x[,3],
#                        Group=c("Control","Control","Control","Control")) 
# 
# 
# 
# pdf("PCA-control_v2.pdf",w=3,h=3)
# ggplot(data=pca.data, aes(x=X, y=Y)) +
#   geom_point(aes(color = Group),size=3) +
#   geom_text_repel(aes(label=Sample, color=Group), show.legend = F) +
#   scale_color_manual(values = c("orange","purple")) +
#   xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) +
#   ylab(paste("PC3 - ", pca.var.per[3], "%", sep = "")) +
#   theme_classic() + theme(legend.position="None") +
#   ggtitle("PCA graph")  
# #stat_ellipse(aes(fill = Group), level = 0.95, geom = 'polygon', alpha=0.1, show.legend = FALSE) +
# #scale_fill_manual(values = c('orange', 'purple'))
# 
# dev.off()
# 
# 
# #############################################################
# ##############################################################
# # PC2 PC3
# 
# pca.data <- data.frame(Sample=rownames(pca$x),
#                        X=pca$x[,2],
#                        Y=pca$x[,3],
#                        Group=c("Control","Control","Control","Control")) 
# 
# 
# 
# pdf("PCA-control_v3.pdf",w=3,h=3)
# ggplot(data=pca.data, aes(x=X, y=Y)) +
#   geom_point(aes(color = Group),size=3) +
#   geom_text_repel(aes(label=Sample, color=Group), show.legend = F) +
#   scale_color_manual(values = c("orange","purple")) +
#   xlab(paste("PC2 - ", pca.var.per[2], "%", sep = "")) +
#   ylab(paste("PC3 - ", pca.var.per[3], "%", sep = "")) +
#   theme_classic() + theme(legend.position="None") +
#   ggtitle("PCA graph")  
# #stat_ellipse(aes(fill = Group), level = 0.95, geom = 'polygon', alpha=0.1, show.legend = FALSE) +
# #scale_fill_manual(values = c('orange', 'purple'))
# 
# dev.off()


# ================PCA analysis with DEP=========================
########## SVM without OI1, Control1 #######################################################3
data_normalized <- da_filtered %>%
  select(name, starts_with("Control"), starts_with("OI")) %>% 

# ============ PCA analysis ================
pca <- prcomp(t(data_normalized %>% select(-OI1) %>% select(-Control1)))
# plot(pca$x, main="after PCA")
pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

library(ggplot2)
library(ggrepel)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       Group=c("Control","Control","Control","OI","OI","OI")) 





# Read in the data
dat <- data.frame(pca$x[,c(1,2)][c("Control2","Control3","Control4","OI2","OI3","OI4"),],Group=c(-1,-1,-1,1,1,1))
# Create a logistic rgression model on the reduced data
svmfit = svm(Group~.,data=dat,kernel="polynomial",coef0=1)

# plot(svmfit,dat)
# Visualize the decision boundaries for training set
X1 = seq(range(pca$x[,1])[1]*2, range(pca$x[,1])[2], len=20)
X2 = seq(range(pca$x[,2])[1]*1.5, range(pca$x[,2])[2]*1.05, len=20)
xgrid = expand.grid(PC1=X1, PC2=X2)
ygrid = predict(svmfit, xgrid)

# ygrid <- predict(svmfit, xgrid, decision.values=T)
pdf("PCA_SVM_3vs3.pdf",w=5,h=5)
plot(dat[,-3], cex=2,col= ifelse(dat$Group == -1, "lightblue","red"),pch = 19,
     main = 'SVM after Principal Component Analysis',
     xlab(paste("PC1 - ", pca.var.per[1], "%", sep = "")) ,
     ylab(paste("PC2 - ", pca.var.per[2], "%", sep = "")))
contour(X1, X2, matrix(as.numeric(ygrid),20,20), add=TRUE,level=0)
contour(X1, X2, matrix(as.numeric(ygrid), 20,20), add=TRUE,level=0.5,lwd=2,lty=2,col="red")
contour(X1, X2, matrix(as.numeric(ygrid), 20,20), add=TRUE,level=-0.5,lwd=2,lty=2,col="blue")
text(dat[,-3],labels=row.names(dat))
dev.off()  







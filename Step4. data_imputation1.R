setwd("C:/Users/oixra/Desktop/diann-4vs4-20230529")
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

sample_names = c("Control1","Control2","Control3","Control4","OI_IFITM5","OI_WNT1","OI_WNT2","OI_WNT3")

# =========== NA statistics ====================
lfq1 <- da %>%
  select(starts_with("Control") | starts_with("OI"))
val1 <- sum(!is.na(lfq1))
na1 <- sum(is.na(lfq1))


# ========== Data Winsorisation ==========
da <- da %>% 
  rowwise() %>% 
  mutate(mean_exp=mean(c_across(starts_with("Control") | starts_with("OI")),na.rm=T)) %>% ungroup()

cut_point_top <- quantile(da$mean_exp, 0.95)
cut_point_bottom <- quantile(da$mean_exp, 0.05)



# =====================过滤缺失值=====================
da <- da %>%
  rowwise() %>%
  mutate(Num_valid_pro_Ctrl = sum(!is.na(c_across(starts_with("Control"))))) %>%
  mutate(Num_valid_pro_OI = sum(!is.na(c_across(starts_with("OI"))))) %>%
  mutate(sum_of_valid_pro = sum(!is.na(c_across(starts_with("Control") | starts_with("OI"))))) %>%
  arrange(-sum_of_valid_pro) %>% ungroup()

da_filtered <- da %>%
  filter(sum_of_valid_pro>4,mean_exp < cut_point_top, mean_exp > cut_point_bottom)
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
plot(imputed)


# =========Data scaling & PCA ============
data_normalized <- da_filtered %>%
  select(starts_with("Control") | starts_with("OI")) %>%
scale(center=T, scale=T)
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
                       Group=c("Control","Control","Control","Control","OI","OI","OI","OI")) 

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


data <- da_filtered %>% 
  select(name | starts_with("Control") | starts_with("OI"))

# DEG
data_deg <- data %>% rowwise() %>%
  mutate(Ctrl_avg = mean(c_across(starts_with("Control")))) %>%
  mutate(Exp_avg = mean(c_across(starts_with("OI")))) %>%
  mutate(`logfc(Exp/Ctrl)`=Exp_avg-Ctrl_avg) 


data_deg$pvalue <- sapply(1:nrow(data_deg), function(i) t.test(data_deg[i,2:5],
                                                               data_deg[i,7:9],var.equal = T,paired=F)$p.value)

data_deg$q.value <- p.adjust(data_deg$pvalue,method = "BH")


##################Output statistical result#############################
write.csv(data_deg,file="step4_results/DEG.csv",row.names = F)



# Strange for p value

# group_list <- factor(group[,1],ordered = F)
# design <- model.matrix(~0+group_list)
# colnames(design) <- levels(group_list)
# rownames(design) <- colnames(expr.log)


# # 线性模型拟合
# fit <- lmFit(expr.log,design)
# 
# # 构建比较矩阵（分析的时候最好输出一下构建的比较矩阵，反了的话有点麻烦）
# cont.matrix <- makeContrasts(contrasts = paste0(unique(group_list),collapse = "-"),levels = design)
# 
# 
# fit2 <- contrasts.fit(fit,cont.matrix) # 构建数据线性模型,计算估计的相关系数和标准差
# fit2 <- eBayes(fit2) # 贝叶斯检验
# tmpOut <- topTable(fit2,coef = 1,n = Inf,adjust.method = "BH",sort.by = "logFC") # 结果
# 
# # sort.by设置结果的排序
# 
# limma.na <- na.omit(tmpOut) # 去除带有NA值的结果
# 
# # 筛选差异表达蛋白质（看自己想用哪个P值）
# dif <- limma.na[limma.na$P.Value <= 0.05 & abs(limma.na$logFC) > log2(2),]或
# dif <- limma.na[limma.na$adj.P.Val <= 0.05 & abs(limma.na$logFC) > log2(2),]
# 
# 
# 
# 
# 
# 
# 








##### volcano plot ##########
cut_off_pvalue=0.05
cut_off_logFC = 0.585

data_deg$change = ifelse(data_deg$pvalue < cut_off_pvalue & abs(data_deg$`logfc(Exp/Ctrl)`) >= cut_off_logFC, 
                   ifelse(data_deg$`logfc(Exp/Ctrl)`>= cut_off_logFC,'Up','Down'),'Stable')

dm<-data_deg[data_deg$pvalue<0.05 & abs(data_deg$`logfc(Exp/Ctrl)`) >= 0.585,]



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
colnames(data_deg)

p<-ggplot(
  data_deg, aes(x = `logfc(Exp/Ctrl)`, y = -log2(pvalue)))+  
  geom_point(aes(color=change),alpha=0.5, size=1) +
  
  #自定义颜色
  scale_color_manual(values = c("#ff4757","#546de5","#d2dae2"), breaks=c("Up","Down","Stable"))+
  
  # Add auxiliary line
  geom_vline(xintercept=c(-0.585,0.585),lty=4,col="black",lwd=0.5) +
  geom_hline(yintercept = -log2(cut_off_pvalue),lty=4,col="black",lwd=0.5) +
  
  
  # Add axis
  labs(x="log2 FoldChange(Exp/Ctrl)",
       y="-log2(pval)")+
  
  # Add labels of some position
  geom_text_repel(data=dm,aes(x=`logfc(Exp/Ctrl)`,y= -log2(pvalue),label=name,color=change),
                  max.overlaps = 100,    # 最大覆盖率，调大该值则不被覆盖，反之。
                  size = 2,    # 字体大小
                  box.padding = unit(0.2, "lines"),  # 标记离box距离
                  point.padding = unit(0.03, "lines"), #标记离点距离
                  segment.color = "black",   # 标记线条的颜色
                  show.legend = FALSE)+
  # xlim(-20,20)+
  # ggtitle("Frontal Cortex")+
  a

pdf("step4_results/volcano_plot.pdf",w=5,h=5)
p
dev.off()

colnames(data_deg)
DEG <- data_deg %>%
  filter(change!="Stable") %>%
  select(name, pvalue, q.value, change)



# =========== PCA版本二 =========================
# #将基因表达值矩阵作个转置，使行为样本，列为基因
# data <- t(data_normalized)
# 
# #我们使用 FactoMineR 包中的方法，实现 PCA 分析和聚类添加
# library(FactoMineR)

# #样本中基因表达值的 PCA 分析
# gene.pca <- PCA(gene, ncp = 2, scale.unit = TRUE, graph = FALSE)
# plot(gene.pca) 




# ============ heat map ================
library(pheatmap)
pheatmap(data_normalized,
         cluster_cols = T,
         angle_col=45)









# 
# library(MSnSet.utils)
# plot_pca(data_normalized)
# # install.packages("devtools")
# # devtools::install_github("PNNL-Comp-Mass-Spec/MSnSet.utils")
# library(MSnSet.utils)


########### The "MSnSet" ##################


# BiocManager::install("MSnbase")


'# missForest
(xmis = da,maxiter = 10,verbose = TRUE)$ximp
# index <- which(methods == "RF")
# results_da[,index] <- imputed_da[,index]

# Write to a new file
library(openxlsx)
sheets = list("log2" = df_log2,"quant_filled" = df_final)
write.xlsx(sheets,"./quant.xlsx")

# write to a new file
# write.xlsx(df_final, "summary.proteins.xlsx",sheetName = "quant",row.names = F)

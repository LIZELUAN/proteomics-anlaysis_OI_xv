# Quality assessment
setwd("C:/Users/oixra/Desktop/diann-4vs4-20230529")
library(openxlsx)
library(ggplot2)
library(tidyverse)


df_ori <- read_csv("txt/report.pg_matrix.csv")

# generate a percentage statistical table for protein
df_sum <- df_ori %>%  
  summarize(across(Genes | starts_with("Control") | starts_with("OI"), ~ sum(!is.na(.x)))) %>% 
  rename(Total=Genes)

library(stringr)

df_ori2 <-read_csv("txt/report.pr_matrix.csv") %>%
  mutate(length=str_length(Precursor.Id))

# generate a percentage statistical table for peptide.
df_sum2 <- df_ori2 %>% 
  summarize(across(Genes | starts_with("Control") | starts_with("OI"), ~ sum(!is.na(.x)))) %>%
  rename(Total=Genes)

# Merge protein and peptide into one table
new_col <- data.frame(Type=c("protein","precursor"))
sta <- cbind(new_col,rbind(df_sum,df_sum2))

write.csv(sta,file="results/table_of_protein_and_peptide.csv",row.names=F)





# transform the table
library(reshape2)
sta_long1<-melt(sta,
                id.vars = c('Type'), 
                measure.va = c("Control1","Control2","Control3","Control4","OI1","OI2","OI3","OI4","Total"),
                variable.name='Group',
                value.name='Count')
sta_long1
write.csv(sta_long1,file="results/table_of_protein_and_peptide2.csv",row.names=F)




# Change the order
glist <- c("Control1","Control2","Control3","Control4","OI1","OI2","OI3","OI4","Total")
sta_long1$Group <- factor(sta_long1$Group,level=rev(glist))

###############################################################
##  The number of protein and peptide types detected  #########
#################### Draw output  #############################
library("ggsci")
pdf("results/identification_types.pdf",h=4,w=6.9)
sta_long1 %>%
  
  ggplot(aes(x=Group, y=Count,fill=Type)) + 
  geom_bar(stat="identity",position = "dodge") +
  geom_text(aes(label=Count), position=position_dodge(width = 0.9),size = 2.6,hjust = -0.05) + 
  # scale_fill_manual(values=c("lightblue", "orange")) +
  scale_fill_simpsons() +
  coord_flip() + theme_classic() + theme(legend.position="top")
dev.off()
###


#####################################################
########### Peptide number per protein  ##############
#####################################################
pep_num_pp <- df_ori %>%
  select(starts_with("Razor + unique peptides ")) %>%
  as.matrix() %>%
  as.vector() %>%
  as.data.frame()
colnames(pep_num_pp) <- "x"
pdf("results/Distribution of peptide count of protein group.pdf",h=5,w=5)
ggplot(pep_num_pp, aes(x = x)) +                           # Modify filling of bars
  geom_histogram(binwidth = 1,col = "black", fill = "red") +
  labs(x = "Peptides count per protein group", y = "Protein groups count") +
  scale_y_continuous(limits = c(0, 2200), breaks = seq(0, 2200, 200)) +
  scale_x_continuous(limits = c(0, 70), breaks = seq(1, 70, 3)) +
  theme_bw()
dev.off()

# hist(pep_num_pp,breaks=70,col = "lightblue",
#   xlab = "peptides count per protein group", ylab = "protein groups count",
#   main = "Distribution of peptide count of protein group")




#######################################################
############## Peptide length distribution ############
######################################################
pdf("results/Peptide length distribution.pdf",h=5,w=8)
ggplot(df_ori2, aes(x = length)) +                           # Modify filling of bars
  geom_histogram(binwidth = 1,col = "black", fill = "red") +
  labs(x = "Peptide length", y = "Count") +
  # scale_y_continuous(limits = c(0, 2200), breaks = seq(0, 2200, 200)) +
  scale_x_continuous(breaks = seq(7,92, 2)) +
  theme_bw()
dev.off()




#############################################
#### NA protein number of each sample   ####
############################################


# generate a percentage statistical table for protein
dfna_sum <- df_ori %>%  
  summarize(across(starts_with("Control") | starts_with("OI"), ~ sum(is.na(.x))))
tot_pro_num <- df_ori %>% nrow()


# generate a percentage statistical table for peptide
dfna_sum2 <- df_ori2 %>% 
  summarize(across(starts_with("Control") | starts_with("OI"), ~ sum(is.na(.x))))

# Merge protein and peptide into one table
new_col <- data.frame(Type=c("protein","precursor"))
staNA <- cbind(new_col,rbind(dfna_sum, dfna_sum2))


# transform the table
library(reshape2)
staNA_long1<-melt(staNA,
                id.vars = c('Type'), 
                measure.va = c("Control1","Control2","Control3","Control4","OI1","OI2","OI3","OI4"),
                variable.name='Group',
                value.name='Count')
staNA_long1

# Change the order
glist <- c("Control1","Control2","Control3","Control4","OI1","OI2","OI3","OI4")
staNA_long1$Group <- factor(staNA_long1$Group,level=glist)


### Draw output########
pdf("results/NA protein number of each sample.pdf",h=4,w=8.1)
staNA_long1 %>%
  filter(Type == "protein") %>%
  ggplot(aes(x=Group,y=Count)) + 
  scale_fill_gradient(low = "lightblue",high = "#0071AA") +
  geom_bar(aes(fill = Count), stat="identity",  position = "dodge", width=0.7) +
  geom_text(aes(label = paste0(round(Count/tot_pro_num*100, digits=2), "%")),
            size=2.6,
            hjust = -0.03) + 
  coord_flip() + 
  ylab("Count of NA proteins") +
  theme_classic() +
  theme(legend.position="None") +
  ggtitle("NA protein number of each sample")+
  theme(plot.title=element_text(face="bold",hjust=0.5))
  

dev.off()



  





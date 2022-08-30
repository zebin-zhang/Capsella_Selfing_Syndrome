#########################################################################
#          PCA -- Zebin Zhang   (ZEBIN_ZHANG@OUTLOOK.COM)
#########################################################################

##  Set work directly and load files.
setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project")
# This file contains reads counts infromation of each gene in all samples
rawTMM <- read.table("InputData/Cbp_diploids_scaled_FRL_masked_RNA_ASE_rename.txt", header=TRUE, row.names=1)

head(rawTMM[,1:36])
diploidsData <- rawTMM[,1:36]
#write.table(diploidsData, file = "InputData/Cbp_diploids_raw_reads_count.txt", quote = F, sep = "\t")
head(rawTMM)
dim(rawTMM)
# This file contails informations of population and tissue in each accession
Accession <- read.table("InputData/DiploidsPhenotypicFile.txt", header = T)
head(Accession)

### Define species
CR_F <- c(Accession$species=="CR" & Accession$tissue=="flower")
CG_F <- c(Accession$species=="CG" & Accession$tissue=="flower")
CO_F <- c(Accession$species=="CO" & Accession$tissue=="flower")

CR_L <- c(Accession$species=="CR" & Accession$tissue=="leaf")
CG_L <- c(Accession$species=="CG" & Accession$tissue=="leaf")
CO_L <- c(Accession$species=="CO" & Accession$tissue=="leaf")

CR_R <- c(Accession$species=="CR" & Accession$tissue=="root")
CG_R <- c(Accession$species=="CG" & Accession$tissue=="root")
CO_R <- c(Accession$species=="CO" & Accession$tissue=="root")

dim(rawTMM)
TMM <- na.omit(rawTMM)
dim(TMM)

## Flower of Cr, Cg, Co
#CPMnorm.Pop$CR_F <- rowMeans(CPMnorm[,CR_F], na.rm = T)
#CPMnorm.Pop$CG_F <- rowMeans(CPMnorm[,CG_F], na.rm = T)
#CPMnorm.Pop$CO_F <- rowMeans(CPMnorm[,CO_F], na.rm = T)
## Leaf of Cr, Cg, Co
#CPMnorm.Pop$CR_L <- rowMeans(CPMnorm[,CR_L], na.rm = T)
#CPMnorm.Pop$CG_L <- rowMeans(CPMnorm[,CG_L], na.rm = T)
#CPMnorm.Pop$CO_L <- rowMeans(CPMnorm[,CO_L], na.rm = T)
## Root of Cr, Cg, Co
#CPMnorm.Pop$CR_R <- rowMeans(CPMnorm[,CR_R], na.rm = T)
#CPMnorm.Pop$CG_R <- rowMeans(CPMnorm[,CG_R], na.rm = T)
#CPMnorm.Pop$CO_R <- rowMeans(CPMnorm[,CO_R], na.rm = T)
#
#row.names(CPMnorm.Pop) <- as.character(CPMnorm.Pop$Gene)
#CPMnorm.Pop$Gene <- NULL
#
#head(CPMnorm.Pop)
#dim(CPMnorm.Pop)
#dim(na.omit(CPMnorm.Pop))
###############################
# PCA in all tissues: Flower, Leaf, Root
# Perform PCA
pca.all <- prcomp(t(TMM))
# NOTE: By default, prcomp() expects the samples to be rows and the genes to be columns
# So the data have to be transpose the matrix using the t() function
# If don't transpose the matrix, you will ultimately get a graph that shows how the genes are related to each other.

# prcomp() returns thress things:
# 1) x -- x contains the principal components (PCs) for drawing a graph.
## plot PC1 and PC2
plot(pca.all$x[,1],pca.all$x[,2])
# 2) sdev -- "standard deviation", to calculate how much variation in the original data each PC accounts for.
## how much variation in the original data PC1 accounts for.
pca.all.var <- pca.all$sdev^2
head(pca.all.var)
## Since the percentage of variation that each PC accounts for is way more interesting than the acutal value,
## we calculate the percentage..
pca.all.var.per <- round(pca.all.var/sum(pca.all.var)*100,1)
## Plotting the percentage is easy with barplot()
barplot(pca.all.var.per, main = "Contribuation of PCs", xlab = "Pricipal Component", ylab = "Percent Variation")
# 3) rotation -- the loading scores rotation
# Format the data the way ggplot2 likes it
dim(pca.all$x)[1]
pca.all.data <- data.frame(ID = 1:dim(pca.all$x)[1])
for (i in 1:dim(pca.all$x)[1]){
  pcs <- paste0("PC", i)
  pca.all.data[[pcs]] <- pca.all$x[,i]
}

pca.all.data$ID <- row.names(pca.all$x)

pca.all.data$ID <- row.names(pca.all.data)
pca.all.data$Tissue <- c(rep("Flower",12), rep("Leaf",12), rep("Root",12))
pca.all.data$Species <- rep(c(rep("CR",4), rep("CG",4), rep("CO", 4)), 3)
pca.all.data$Species <- factor(pca.all.data$Species, levels = c("CG","CR","CO"))
pca.all.data$MatingSystem <- rep(c(rep("Selfer",4), rep("Outcrosser",4), rep("Selfer", 4)), 3)


ZscoreAll <- pca.all.data
head(ZscoreAll)
for(i in 2:13){
  ZscoreAll[,i] <- scale(ZscoreAll[,i])
}
head(ZscoreAll)

# Tissue
p1 <- ggplot(ZscoreAll, aes(PC1, PC2, col = Tissue, fill = Species))  +
  geom_point(size = 6, stroke = 2, aes(shape = Species)) +
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("#fa990e","#155800","#089400"))  +
  scale_fill_manual(values = c("#f90052","darkviolet","#006beb")) +
  # stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.all.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.all.var.per[2], "%", sep = "")) +
  # stat_ellipse(geom="polygon", aes(fill = tissue), # add frame
  #              alpha = 0.05, 
  #              show.legend = FALSE, 
  #              level = 0.90) + 
  scale_x_continuous(breaks=seq(-3, 2, 1), limits=c(-3, 2)) +
  scale_y_continuous(breaks=seq(-2, 2, 1), limits=c(-2, 2)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold")) +
  theme(legend.position = "none") +
  ggtitle("by Tissue") +
  theme(plot.title = element_text(size=16,face="bold"))

p1_i <- ggplot(ZscoreFlower, aes(PC1, PC2, col = Species, fill = Species))  +
  geom_point(size = 5, stroke = 1, aes(shape = Species)) +
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("#fa990e","#fa990e","#fa990e"))  +
  scale_fill_manual(values = c("#fa990e","#fa990e","#fa990e")) +
  xlim(-2,2) + ylim(-2,2)+
  # stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  #xlab("PC1") +
  #ylab("PC2") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

p_tissue <- p1 + annotation_custom(ggplotGrob(p1_i), xmin = -3.41, xmax = -1, 
                  ymin = -2.35, ymax = 0)




# Species
p2 <- ggplot(ZscoreAll, aes(PC1, PC2, col = Species, fill = Species))  +
  geom_point(size = 6, stroke = 2, aes(shape = Species)) +
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("#f90052","darkviolet","#006beb"))  +
  scale_fill_manual(values = c("#f90052","darkviolet","#006beb")) +
  # stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.all.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.all.var.per[2], "%", sep = "")) +
  # stat_ellipse(geom="polygon", aes(fill = tissue), # add frame
  #              alpha = 0.05, 
  #              show.legend = FALSE, 
  #              level = 0.90) + 
  scale_x_continuous(breaks=seq(-3, 2, 1), limits=c(-3, 2)) +
  scale_y_continuous(breaks=seq(-2, 2, 1), limits=c(-2, 2)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold")) +
  theme(legend.position = "none") +
  ggtitle("by Species") +
  theme(plot.title = element_text(size=16,face="bold"))

p2_i <- ggplot(ZscoreFlower, aes(PC1, PC2, col = Species, fill = Species))  +
  geom_point(size = 5, stroke = 1, aes(shape = Species)) +
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("#f90052","darkviolet","#006beb"))  +
  scale_fill_manual(values = c("#f90052","darkviolet","#006beb")) +
  xlim(-2,2) + ylim(-2,2)+
  # stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  #xlab("PC1") +
  #ylab("PC2") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

p_species <- p2 + annotation_custom(ggplotGrob(p2_i), xmin = -3.41, xmax = -1, 
                                   ymin = -2.35, ymax = 0)

# Mating System
p3 <- ggplot(ZscoreAll, aes(PC1, PC2, col = MatingSystem, fill = Species))  +
  geom_point(size = 6, stroke = 2, aes(shape = Species)) +
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("#007a79","#976f4f"))  +
  #scale_fill_manual(values = c("#f90052","darkviolet","#006beb")) +
  # stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.all.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.all.var.per[2], "%", sep = "")) +
  scale_x_continuous(breaks=seq(-3, 2, 1), limits=c(-3, 2)) +
  scale_y_continuous(breaks=seq(-2, 2, 1), limits=c(-2, 2)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold")) +
  theme(legend.position = "none") +
  ggtitle("by Mating System") +
  theme(plot.title = element_text(size=16,face="bold"))

p3_i <- ggplot(ZscoreFlower, aes(PC1, PC2, col = Species, fill = Species))  +
  geom_point(size = 5, stroke = 1, aes(shape = Species)) +
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("#007a79","#976f4f","#976f4f"))  +
  scale_fill_manual(values = c("#007a79","#976f4f","#976f4f")) +
  xlim(-2,2) + ylim(-2,2)+
  # stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  #xlab("PC1") +
  #ylab("PC2") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

p_MST <- p3 + annotation_custom(ggplotGrob(p3_i), xmin = -3.41, xmax = -1, 
                                    ymin = -2.35, ymax = 0)


library(grid)
grid.newpage()
library(ggpubr)
ggarrange(p_tissue, p_species, p_MST, ncol=3, nrow=1) # set the common legend

###############################
## PCA in Flowers
### Define flower
getwd()
TMM.Flower <- read.table("InputData/TMM_Diploids_F.txt", header = T, sep = "\t")
head(TMM.Flower)
dim(TMM.Flower)
dim(na.omit(TMM.Flower))
### perform PCA in Flower
pca.flower <- prcomp(t(na.omit(TMM.Flower)))
# 1) X
## plot PC1 and PC2
plot(pca.flower$x[,1],pca.flower$x[,2])
# 2) sdev 
pca.flower.var <- pca.flower$sdev^2
head(pca.flower.var)
## calculate the percentage..
pca.flower.var.per <- round(pca.flower.var/sum(pca.flower.var)*100,1)
## Plotting the percentage is easy with barplot()
barplot(pca.flower.var.per, main = "Contribuation of PCs in Flowers", xlab = "Pricipal Component", ylab = "Percent Variation")
# 3) rotation -- the loading scores rotation

### Draw graph in Flowers
library(ggplot2)
# Format the data the way ggplot2 likes it
dim(pca.flower$x)[1]
pca.flower.data <- data.frame(ID = 1:dim(pca.flower$x)[1])
for (i in 1:dim(pca.flower$x)[1]){
  pcs <- paste0("PC", i)
  pca.flower.data[[pcs]] <- pca.flower$x[,i]
}

pca.flower.data$ID <- row.names(pca.flower$x)

pca.flower.data$ID <- row.names(pca.flower.data)
pca.flower.data$tissue <- "Flower"
pca.flower.data$Species <- c(rep("CG",4), rep("CR",4), rep("CO", 4))
pca.flower.data$Species <- factor(pca.flower.data$Species, levels = c("CG","CR","CO"))
#pca.flower.data$species <- c(c("CR","CG","CO"))

ZscoreFlower <- pca.flower.data
head(ZscoreFlower)
for(i in 2:13){
  ZscoreFlower[,i] <- scale(ZscoreFlower[,i])
}
head(ZscoreFlower)


# graph parttern 1
ggplot(data = pca.flower.data, aes(x=PC1, y=PC2, label=ID)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.flower.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.flower.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Reads count in Flowers, Leaves, and Roots")
# graph parttern 2
# PC1 vs PC2
ggplot(ZscoreFlower, aes(PC1, PC2, col = Species, fill = Species))  +
  geom_point(size = 6, stroke = 2, aes(shape = Species)) +
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("#f90052","darkviolet","#006beb"))  +
  scale_fill_manual(values = c("#f90052","darkviolet","#006beb")) +
  # stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.flower.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.flower.var.per[2], "%", sep = "")) +
  stat_ellipse(geom="polygon", aes(fill = Species), # add frame
               alpha = 0.05, 
               show.legend = FALSE, 
               level = 0.90) + 
  scale_x_continuous(breaks=seq(-2, 2, 1), limits=c(-2.1, 2.1)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold")) +
  theme(legend.position = "none") +
  ggtitle("Flower") +
  theme(plot.title = element_text(size=16,face="bold"))


###############################
## PCA in Leaves
### Define leaf
TMM.Leaf <- read.table("InputData/TMM_Diploids_L.txt", header = T, sep = "\t")
head(TMM.Leaf)
dim(TMM.Leaf)
dim(na.omit(TMM.Leaf))
### perform PCA in Leaf
pca.leaf <- prcomp(t(na.omit(TMM.Leaf)))
summary(pca.leaf)
# 1) X
## plot PC1 and PC2
plot(pca.leaf$x[,1],pca.leaf$x[,2])
# 2) sdev 
pca.leaf.var <- pca.leaf$sdev^2
head(pca.leaf.var)
## calculate the percentage..
pca.leaf.var.per <- round(pca.leaf.var/sum(pca.leaf.var)*100,1)
## Plotting the percentage is easy with barplot()
barplot(pca.leaf.var.per, main = "Contribuation of PCs in Leaves", xlab = "Pricipal Component", ylab = "Percent Variation")
# 3) rotation -- the loading scores rotation

### Draw graph in Leaves
library(ggplot2)
# Format the data the way ggplot2 likes it
dim(pca.leaf$x)[1]
pca.leaf.data <- data.frame(ID = 1:dim(pca.leaf$x)[1])
for (i in 1:dim(pca.leaf$x)[1]){
  pcs <- paste0("PC", i)
  pca.leaf.data[[pcs]] <- pca.leaf$x[,i]
}

pca.leaf.data$ID <- row.names(pca.leaf$x)

pca.leaf.data$tissue <- "Leaf"
pca.leaf.data$Species <- c(rep("CG",4), rep("CR",4), rep("CO", 4))
pca.leaf.data$Species <- factor(pca.leaf.data$Species, levels = c("CG","CR","CO"))

ZscoreLeaf <- pca.leaf.data
head(ZscoreLeaf)
for(i in 2:13){
  ZscoreLeaf[,i] <- scale(ZscoreLeaf[,i])
}
head(ZscoreLeaf)

# graph parttern 1
ggplot(data = pca.leaf.data, aes(x=PC1, y=PC2, label=ID)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.leaf.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.leaf.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Reads count in Leaves, Leaves, and Leaves")
# graph parttern 2
# PC1 vs PC2
p2 <- ggplot(ZscoreLeaf, aes(PC1, PC2, col = Species, fill = Species))  +
  geom_point(size = 6, stroke = 2, aes(shape = Species)) +
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("#f90052","darkviolet","#006beb"))  +
  scale_fill_manual(values = c("#f90052","darkviolet","#006beb")) +
  #stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.leaf.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.leaf.var.per[2], "%", sep = "")) +
  stat_ellipse(geom="polygon", aes(fill = Species), # add frame
               alpha = 0.05, 
               show.legend = FALSE, 
               level = 0.85) + 
  xlim(-4,4) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold")) +
  theme(legend.position = "none") +
  ggtitle("Leaf")+
  theme(plot.title = element_text(size=16,face="bold"))


###############################
## PCA in Roots
### Define root

TMM.Root <- read.table("InputData/TMM_Diploids_R.txt", header = T, sep = "\t")
head(TMM.Root)
dim(TMM.Root)
dim(na.omit(TMM.Root))
### perform PCA in Root
pca.root <- prcomp(t(na.omit(TMM.Root)))
# 1) X
## plot PC1 and PC2
plot(pca.root$x[,1],pca.root$x[,2])
# 2) sdev 
pca.root.var <- pca.root$sdev^2
head(pca.root.var)
## calculate the percentage..
pca.root.var.per <- round(pca.root.var/sum(pca.root.var)*100,1)
## Plotting the percentage is easy with barplot()
barplot(pca.root.var.per, main = "Contribuation of PCs in Roots", xlab = "Pricipal Component", ylab = "Percent Variation")
# 3) rotation -- the loading scores rotation

### Draw graph in Roots
library(ggplot2)
# Format the data the way ggplot2 likes it
dim(pca.root$x)[1]
pca.root.data <- data.frame(ID = 1:dim(pca.root$x)[1])
for (i in 1:dim(pca.root$x)[1]){
  pcs <- paste0("PC", i)
  pca.root.data[[pcs]] <- pca.root$x[,i]
}

pca.root.data$ID <- row.names(pca.root$x)
pca.root.data$tissue <- "Root"
pca.root.data$Species <- c(rep("CG",4), rep("CR",4), rep("CO", 4))
pca.root.data$Species <- factor(pca.root.data$Species, levels = c("CG","CR","CO"))

ZscoreRoot <- pca.root.data
head(ZscoreRoot)
for(i in 2:13){
  ZscoreRoot[,i] <- scale(ZscoreRoot[,i])
}
head(ZscoreRoot)
# graph parttern 1
ggplot(data = ZscoreRoot, aes(x=PC1, y=PC2, label=ID)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.root.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.root.var.per[2], "%", sep = "")) +
  theme_bw() +
  ggtitle("Reads count in Roots")
# graph parttern 2
p3 <- ggplot(ZscoreRoot, aes(PC1, PC2, col = Species, fill = Species))  +
  geom_point(size = 6, stroke = 2, aes(shape = Species)) +
  scale_shape_manual(values=c(8,10,4))+
  scale_color_manual(values = c("#f90052","darkviolet","#006beb"))  +
  scale_fill_manual(values = c("#f90052","darkviolet","#006beb")) +
  #stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
  xlab(paste("PC1 - ", pca.root.var.per[1], "%", sep = "")) +
  ylab(paste("PC2 - ", pca.root.var.per[2], "%", sep = "")) +
  stat_ellipse(geom="polygon", aes(fill = Species), # add frame
               alpha = 0.05, 
               show.legend = FALSE, 
               level = 0.85) + 
  xlim(-2,2) +
  #labs(fill = "Species") +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        panel.border = element_rect(fill= "transparent")) +
  theme(axis.text=element_text(size=16,face="bold"),
        axis.title=element_text(size=16,face="bold"),
        legend.title = element_text(size=16, face = "bold"),
        legend.text = element_text(size=16)) +
  theme(legend.position = "none") +
  ggtitle("Root") +
  theme(plot.title = element_text(size=16,face="bold"))
  


# # only plot legend
# library(ggplot2) 
# library(grid)
# library(gridExtra) 
# my_plot <- ggplot(pca.root.data, aes(PC1, PC2, col = population, fill = population))  +
#   geom_point(size = 3, aes(shape = population)) +
#   scale_shape_manual(values=c(8,10,4))+
#   scale_color_manual(values = c("indianred1","darkviolet","blue"))  +
#   scale_fill_manual(values = c("indianred1","darkviolet","blue")) +
#   #stat_ellipse(aes(x= PC1, y=PC2, group = species), geom = "polygon", alpha = 0.1) +
#   xlab(paste("PC1 - ", pca.root.var.per[1], "%", sep = "")) +
#   ylab(paste("PC2 - ", pca.root.var.per[2], "%", sep = "")) +
#   #stat_ellipse(geom="polygon", aes(fill = population), # add frame
#   #             alpha = 0.05, 
#   #             show.legend = FALSE, 
#   #             level = 0.90) + 
#   theme_minimal() +
#   theme(panel.grid = element_blank(), 
#         panel.border = element_rect(fill= "transparent")) +
#   ggtitle("Gene expression in Roots") 
# # Using the cowplot package
# legend <- cowplot::get_legend(my_plot)
# p4 <- grid.draw(legend)

library(grid)
grid.newpage()
library(ggpubr)
ggarrange(p3, p2, p1, ncol=3, nrow=1, common.legend = TRUE, legend="right") # set the common legend




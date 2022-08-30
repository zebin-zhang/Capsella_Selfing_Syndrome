#########################################################

#  D index pairwise comparison amongst Cg, Cr, and Co.  #

#########################################################

## Set work directory and load files
setwd("~/Desktop/My_Computer/Uppsala_University/Capsella_Project/")
#### This file contains reads counts infromation of each gene in all samples
raw.RC <- read.table("InputData/Cbp_diploids_scaled_FRL_masked_RNA_ASE.csv", header=TRUE, row.names="gene")
raw.RC <- read.table("InputData/Cbp_diploids_scaled_FRL_masked_RNA_ASE_rename.txt", header=T, row.names="gene")
head(raw.RC)
dim(raw.RC)
# CPMnorm <- read.table("Capsella_ReadsCount_FRL_masked_RNA_ASE_FILTERED_CPM.csv", header = T, sep = "\t")
# head(CPMnorm)
# dim(CPMnorm)
##### This file contains information of population and tissue in each sample
Accession <- read.csv("InputData/28genomes_annot_R.csv")
head(Accession)

### Define populations
CR_F <- c(Accession$populations=="CR" & Accession$tissue=="flower")
CG_F <- c(Accession$populations=="CG" & Accession$tissue=="flower")
CO_F <- c(Accession$populations=="CO" & Accession$tissue=="flower")

CR_L <- c(Accession$populations=="CR" & Accession$tissue=="leaf")
CG_L <- c(Accession$populations=="CG" & Accession$tissue=="leaf")
CO_L <- c(Accession$populations=="CO" & Accession$tissue=="leaf")

CR_R <- c(Accession$populations=="CR" & Accession$tissue=="root")
CG_R <- c(Accession$populations=="CG" & Accession$tissue=="root")
CO_R <- c(Accession$populations=="CO" & Accession$tissue=="root")

geneMis <- 1 
geneKeep <- cbind(  ((rowSums(raw.RC[,CG_F], na.rm = T) >= 1) >= geneMis) &
                    ((rowSums(raw.RC[,CR_F], na.rm = T) >= 1) >= geneMis) &
                    ((rowSums(raw.RC[,CO_F], na.rm = T) >= 1) >= geneMis) &
                    
                    ((rowSums(raw.RC[,CG_L], na.rm = T) >= 1) >= geneMis) &
                    ((rowSums(raw.RC[,CR_L], na.rm = T) >= 1) >= geneMis) &
                    ((rowSums(raw.RC[,CO_L], na.rm = T) >= 1) >= geneMis) &
                    
                    ((rowSums(raw.RC[,CG_R], na.rm = T) >= 1) >= geneMis) &
                    ((rowSums(raw.RC[,CR_R], na.rm = T) >= 1) >= geneMis) &
                    ((rowSums(raw.RC[,CO_R], na.rm = T) >= 1) >= geneMis) 
                     )
raw.RC.Keep <- raw.RC[geneKeep,]
message(paste("Number of genes before", dim(raw.RC)[1], "\nNumber of genes after", dim(raw.RC.Keep)[1]))

### Calculate the library sizes 
# Library size for Cbp = (Co + Cg)*0.5:
libColsums <- colSums(raw.RC.Keep, na.rm = T) # add one because later we will have to log values.
data.lib.size <- c()

for (i in c(1:c(length(libColsums)))){
    data.lib.size[i] <- libColsums[i]
}

barplot(data.lib.size/1000000, names.arg=Accession$accession, las =2, border = F,
        col= as.character(Accession$col), ylab = "Library size (million reads)", cex.names=0.5)
# pdf("Lib-size_norm_FRL_All2.pdf", height=2, width=5)
# par(mar=c(5, 4, 1, 0), cex=0.5)
# barplot(data.lib.size/1000000, names.arg=Accession$accession, las =2, border = F,
#      col= tissuesCol, ylab = "Library size (million reads)", cex.names=0.7)
# dev.off()


### Count per million (CPM) normilazation
CPMnorm <- raw.RC.Keep
logCPMnorm <- raw.RC.Keep
for (i in c(1:c(length(data.lib.size)))){
  message(data.lib.size[i])
  CPMnorm[,i] <- (raw.RC.Keep[,i] / data.lib.size[i]) * 1000000 # why needs to +1 in here?
  logCPMnorm[,i] <- log((raw.RC.Keep[,i] /data.lib.size[i])*1000000) # why needs to +1 in here?
}
head(CPMnorm)
dim(CPMnorm)
head(logCPMnorm)
dim(logCPMnorm)

##############################
# Define populations
CR_F <- c(Accession$populations=="CR" & Accession$tissue=="flower")
CG_F <- c(Accession$populations=="CG" & Accession$tissue=="flower")
CO_F <- c(Accession$populations=="CO" & Accession$tissue=="flower")

CR_L <- c(Accession$populations=="CR" & Accession$tissue=="leaf")
CG_L <- c(Accession$populations=="CG" & Accession$tissue=="leaf")
CO_L <- c(Accession$populations=="CO" & Accession$tissue=="leaf")

CR_R <- c(Accession$populations=="CR" & Accession$tissue=="root")
CG_R <- c(Accession$populations=="CG" & Accession$tissue=="root")
CO_R <- c(Accession$populations=="CO" & Accession$tissue=="root")

##############################

CPMnorm.Pop <- data.frame(Gene=row.names(CPMnorm))
head(CPMnorm.Pop)
# Use the MEAN value of each sample to represents as population value
# Flower of Cr, Cg, Co
CPMnorm.Pop$CR_F <- rowMeans(CPMnorm[,CR_F], na.rm = T)
CPMnorm.Pop$CG_F <- rowMeans(CPMnorm[,CG_F], na.rm = T)
CPMnorm.Pop$CO_F <- rowMeans(CPMnorm[,CO_F], na.rm = T)
# Leaf of Cr, Cg, Co
CPMnorm.Pop$CR_L <- rowMeans(CPMnorm[,CR_L], na.rm = T)
CPMnorm.Pop$CG_L <- rowMeans(CPMnorm[,CG_L], na.rm = T)
CPMnorm.Pop$CO_L <- rowMeans(CPMnorm[,CO_L], na.rm = T)
# Root of Cr, Cg, Co
CPMnorm.Pop$CR_R <- rowMeans(CPMnorm[,CR_R], na.rm = T)
CPMnorm.Pop$CG_R <- rowMeans(CPMnorm[,CG_R], na.rm = T)
CPMnorm.Pop$CO_R <- rowMeans(CPMnorm[,CO_R], na.rm = T)

row.names(CPMnorm.Pop) <- as.character(CPMnorm.Pop$Gene)
CPMnorm.Pop$Gene <- NULL
head(CPMnorm.Pop)
dim(CPMnorm.Pop)
dim(na.omit(CPMnorm.Pop))
str(CPMnorm.Pop)

## Comparasion 
CPMnorm.Pop.Comp <- CPMnorm.Pop[,c("CG_F","CR_F","CO_F","CG_L","CR_L","CO_L","CG_R","CR_R","CO_R")]
head(CPMnorm.Pop.Comp)
# First calcuate the absolute value cross species.
### Flower
CPMnorm.Pop.Comp$abs.F_CgCo <- abs(CPMnorm.Pop.Comp[,"CG_F"] - CPMnorm.Pop.Comp[,"CO_F"])
CPMnorm.Pop.Comp$abs.F_CgCr <- abs(CPMnorm.Pop.Comp[,"CG_F"] - CPMnorm.Pop.Comp[,"CR_F"])
CPMnorm.Pop.Comp$abs.F_CrCo <- abs(CPMnorm.Pop.Comp[,"CR_F"] - CPMnorm.Pop.Comp[,"CO_F"])
### Leaf
CPMnorm.Pop.Comp$abs.L_CgCo <- abs(CPMnorm.Pop.Comp[,"CG_L"] - CPMnorm.Pop.Comp[,"CO_L"])
CPMnorm.Pop.Comp$abs.L_CgCr <- abs(CPMnorm.Pop.Comp[,"CG_L"] - CPMnorm.Pop.Comp[,"CR_L"])
CPMnorm.Pop.Comp$abs.L_CrCo <- abs(CPMnorm.Pop.Comp[,"CR_L"] - CPMnorm.Pop.Comp[,"CO_L"])
### Root
CPMnorm.Pop.Comp$abs.R_CgCo <- abs(CPMnorm.Pop.Comp[,"CG_R"] - CPMnorm.Pop.Comp[,"CO_R"])
CPMnorm.Pop.Comp$abs.R_CgCr <- abs(CPMnorm.Pop.Comp[,"CG_R"] - CPMnorm.Pop.Comp[,"CR_R"])
CPMnorm.Pop.Comp$abs.R_CrCo <- abs(CPMnorm.Pop.Comp[,"CR_R"] - CPMnorm.Pop.Comp[,"CO_R"])

## Pick the samllar value of CgCr & CgCo, as the differences from outcorsser to selfer
## Use the value of CrCo as the difference within selfer
# install.packages("matrixStats")
library(matrixStats)
# Flower
CPMnorm.Pop.Comp$F.oc2sf <- rowMins(as.matrix(CPMnorm.Pop.Comp[,c("abs.F_CgCo","abs.F_CgCr")]), na.rm = TRUE)
CPMnorm.Pop.Comp$F.sf2sf <- CPMnorm.Pop.Comp$abs.F_CrCo
# Leaf
CPMnorm.Pop.Comp$L.oc2sf <- rowMins(as.matrix(CPMnorm.Pop.Comp[,c("abs.L_CgCo","abs.L_CgCr")]), na.rm = TRUE)
CPMnorm.Pop.Comp$L.sf2sf <- CPMnorm.Pop.Comp$abs.L_CrCo
# Root
CPMnorm.Pop.Comp$R.oc2sf <- rowMins(as.matrix(CPMnorm.Pop.Comp[,c("abs.R_CgCo","abs.R_CgCr")]), na.rm = TRUE)
CPMnorm.Pop.Comp$R.sf2sf <- CPMnorm.Pop.Comp$abs.R_CrCo
# Since we want to find some genes simliar expressed within selfer but different from selfer to outcrosser.
# So we require the value of sf2sf the smaller the better and value of oc2sf the bigger the better
# Here I creat a index of differences between two mating system:
# D.mate =  oc2sf - sf2sf
CPMnorm.Pop.Comp$F_D.mate <- CPMnorm.Pop.Comp$F.oc2sf - CPMnorm.Pop.Comp$F.sf2sf
CPMnorm.Pop.Comp$L_D.mate <- CPMnorm.Pop.Comp$L.oc2sf - CPMnorm.Pop.Comp$L.sf2sf
CPMnorm.Pop.Comp$R_D.mate <- CPMnorm.Pop.Comp$R.oc2sf - CPMnorm.Pop.Comp$R.sf2sf
# Considering the reads count mapped to each gene are different, it's hard to speak the same D.mate will have the same effect on
# mating type transformation from outcrosser to selfer. For instance:
# gene1  Cg = 10  Cr = 5  Co = 6, so the oc2sf=10-6=4, sf2sf=6-5=1, the D.mate = 4-1 =3
# gene2  Cg = 1000  Cr = 995  Co = 996, so the oc2sf=1000-996=4, sf2sf=996-995=1, the D.mate = 4-1 =3
# For gene1, it does have differences between outcrosser and selfer, but not gene2.
# So, to minimize the bias form reads count, we shink the D value to from -1 to 1 by devide to the Max read count of that gene across species.
# D = D.mate / CPM.species,  where CPM.sepcies = Max(CPM.Cg, CPM.Cr, CPM.Co)
# When D = 0, means totally nothing difference between outcrosser vs. selfer and selfer vs. selfer;
# when 0 < D <= 1, means the difference of outcrosser vs. selfer is bigger than the difference of selfer vs. selfer;
# when -1 <= D < 0, means the difference within selfer is bigger than differences between selfer and outcrosser. 

# Maximum CPM amongst thress species per tissue.
CPMnorm.Pop.Comp$F_CPM.species <- rowMaxs(as.matrix(CPMnorm.Pop.Comp[,c("CG_F","CR_F","CO_F")]), na.rm = TRUE)
CPMnorm.Pop.Comp$L_CPM.species <- rowMaxs(as.matrix(CPMnorm.Pop.Comp[,c("CG_L","CR_L","CO_L")]), na.rm = TRUE)
CPMnorm.Pop.Comp$R_CPM.species <- rowMaxs(as.matrix(CPMnorm.Pop.Comp[,c("CG_R","CR_R","CO_R")]), na.rm = TRUE)

# Perform the calculation
### Flower
CPMnorm.Pop.Comp$F.D <- CPMnorm.Pop.Comp$F_D.mate / CPMnorm.Pop.Comp$F_CPM.species
CPMnorm.Pop.Comp$L.D <- CPMnorm.Pop.Comp$L_D.mate / CPMnorm.Pop.Comp$L_CPM.species
CPMnorm.Pop.Comp$R.D <- CPMnorm.Pop.Comp$R_D.mate / CPMnorm.Pop.Comp$R_CPM.species

# HIST
hist(CPMnorm.Pop.Comp$F.D, breaks = 100)
hist(CPMnorm.Pop.Comp$L.D, breaks = 100)
hist(CPMnorm.Pop.Comp$R.D, breaks = 100)

head(CPMnorm.Pop.Comp)
# Statistic how many genes with D value bigger than 0
dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$F.D > 0,])[1]
dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$L.D > 0,])[1]
dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$R.D > 0,])[1]
message("In Flower, The number of genes with higher differences between outcrosser and selfer than differnces within selfer are: ", dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$F.D > 0,])[1])
message("In Leaf, The number of genes with higher differences between outcrosser and selfer than differnces within selfer are: ", dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$L.D > 0,])[1])
message("In Root, The number of genes with higher differences between outcrosser and selfer than differnces within selfer are: ", dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$R.D > 0,])[1])

# Pick the top 5% gene
library(ggplot2)
library(magrittr)
library(dplyr)
#library(ps)
#install.packages("hrbrthemes")
#library(hrbrthemes)
CPMnorm.Pop.Comp.Flower.top0.05 <- CPMnorm.Pop.Comp %>%
  subset(F.D > 0) %>%
  subset(F.D > quantile(F.D, prob = 1 - 5/100))
CPMnorm.Pop.Comp.Leaf.top0.05 <- CPMnorm.Pop.Comp %>%
  subset(L.D > 0) %>%
  subset(L.D > quantile(L.D, prob = 1 - 5/100))
CPMnorm.Pop.Comp.Root.top0.05 <- CPMnorm.Pop.Comp %>%
  subset(R.D > 0) %>%
  subset(R.D > quantile(R.D, prob = 1 - 5/100))
# Pick the cutoff value
min(CPMnorm.Pop.Comp.Flower.top0.05$F.D)  
min(CPMnorm.Pop.Comp.Leaf.top0.05$L.D)  
min(CPMnorm.Pop.Comp.Root.top0.05$R.D)  
message("In flower, The cutoff value of top 5% is: ", min(CPMnorm.Pop.Comp.Flower.top0.05$F.D))
message("In leaf, The cutoff value of top 5% is: ", min(CPMnorm.Pop.Comp.Leaf.top0.05$L.D))
message("In root, The cutoff value of top 5% is: ", min(CPMnorm.Pop.Comp.Root.top0.05$R.D))

# Group genes by D index
CPMnorm.Pop.Comp$F.D.group <- "D.outcrosser < D.selfer"
CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$F.D > 0,]$F.D.group <- "D.outcrosser > D.selfer"
CPMnorm.Pop.Comp[(CPMnorm.Pop.Comp$F.D > min(CPMnorm.Pop.Comp.Flower.top0.05$F.D)),]$F.D.group <- "Top 5%"
CPMnorm.Pop.Comp$F.D.group <- factor(CPMnorm.Pop.Comp$F.D.group, levels = c("Top 5%", "D.outcrosser > D.selfer", "D.outcrosser < D.selfer"))
levels(CPMnorm.Pop.Comp$F.D.group)

CPMnorm.Pop.Comp$L.D.group <- "D.outcrosser < D.selfer"
CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$L.D > 0,]$L.D.group <- "D.outcrosser > D.selfer"
CPMnorm.Pop.Comp[(CPMnorm.Pop.Comp$L.D > min(CPMnorm.Pop.Comp.Leaf.top0.05$L.D)),]$L.D.group <- "Top 5%"
CPMnorm.Pop.Comp$L.D.group <- factor(CPMnorm.Pop.Comp$L.D.group, levels = c("Top 5%", "D.outcrosser > D.selfer", "D.outcrosser < D.selfer"))
levels(CPMnorm.Pop.Comp$L.D.group)

CPMnorm.Pop.Comp$R.D.group <- "D.outcrosser < D.selfer"
CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$R.D > 0,]$R.D.group <- "D.outcrosser > D.selfer"
CPMnorm.Pop.Comp[(CPMnorm.Pop.Comp$R.D > min(CPMnorm.Pop.Comp.Root.top0.05$R.D)),]$R.D.group <- "Top 5%"
CPMnorm.Pop.Comp$R.D.group <- factor(CPMnorm.Pop.Comp$R.D.group, levels = c("Top 5%", "D.outcrosser > D.selfer", "D.outcrosser < D.selfer"))
levels(CPMnorm.Pop.Comp$R.D.group)

head(CPMnorm.Pop.Comp)
mean(CPMnorm.Pop.Comp$F.D)
median(CPMnorm.Pop.Comp$F.D)
mean(CPMnorm.Pop.Comp$L.D)
median(CPMnorm.Pop.Comp$L.D)
mean(CPMnorm.Pop.Comp$R.D)
median(CPMnorm.Pop.Comp$R.D)


p1 <- ggplot(CPMnorm.Pop.Comp, aes(x=F.D, color = "black", fill = F.D.group)) +
  geom_histogram( color="black", alpha=0.6, binwidth = 0.04) +
  scale_fill_manual(values=c("red", "lightpink","gray90")) +
  labs(x = "Changes of differences from (Outcrosser vs. selfer) to (within selfer)",
       y = "Gene counts", title = "Flower") +
  labs(fill="") + # remove the title for legend
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), # put title to centre
        legend.justification=c(1,1), legend.position=c(0.99,0.99)) # put the legend to top right (1,1), left bottom(0,0)
p2 <- ggplot(CPMnorm.Pop.Comp, aes(x=L.D, color = "black", fill = L.D.group)) +
  geom_histogram( color="black", alpha=0.6, binwidth = 0.04) +
  scale_fill_manual(values=c("green4", "#D7FFD7","gray90")) +
  labs(x = "Changes of differences from (Outcrosser vs. selfer) to (within selfer)",
       y = "Gene counts", title = "Leaf") +
  labs(fill="") + # remove the title for legend
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), # put title to centre
        legend.justification=c(1,1), legend.position=c(0.99,0.99)) # put the legend to top right (1,1), left bottom(0,0)
p3 <- ggplot(CPMnorm.Pop.Comp, aes(x=R.D, color = "black", fill = R.D.group)) +
  geom_histogram( color="black", alpha=0.6, binwidth = 0.04) +
  scale_fill_manual(values=c("purple4", "#EACDFE","gray90")) +
  labs(x = "Changes of differences from (Outcrosser vs. selfer) to (within selfer)",
       y = "Gene counts", title = "Root") +
  labs(fill="") + # remove the title for legend
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), # put title to centre
        legend.justification=c(1,1), legend.position=c(0.99,0.99)) # put the legend to top right (1,1), left bottom(0,0)
library("grid")
grid.newpage()
pushViewport(viewport(layout = grid.layout(1,3)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(1,1))
print(p2, vp = vplayout(1,2))
print(p3, vp = vplayout(1,3))


# Plot Ridgeline chart
# library
library(ggridges)
library(ggplot2)
library(tidyverse)
library(viridis)

rgF <- CPMnorm.Pop.Comp[,c("F.D","F.D.group")]
rgF$Tissue <- "Flower"
names(rgF) <- c("D_index","D_group","Tissue")
rgF$Gene.name <- rownames(rgF)
head(rgF)

rgL <- CPMnorm.Pop.Comp[,c("L.D","L.D.group")]
rgL$Tissue <- "Leaf"
names(rgL) <- c("D_index","D_group","Tissue")
rgL$Gene.name <- rownames(rgL)
head(rgL)

rgR <- CPMnorm.Pop.Comp[,c("R.D","R.D.group")]
rgR$Tissue <- "Root"
names(rgR) <- c("D_index","D_group","Tissue")
rgR$Gene.name <- rownames(rgR)
head(rgR)
dim(rgR)

Ridgeline <- rbind(rgF,rgL, rgR)
head(Ridgeline)
dim(Ridgeline)
str(Ridgeline)
Ridgeline$Tissue <- factor(Ridgeline$Tissue, levels = c("Root","Leaf","Flower"))

ggplot(Ridgeline, aes(x=D_index, y=Tissue,  fill=Tissue)) +
  geom_density_ridges(alpha=0.8,
                      quantile_lines=TRUE,
                      quantile_fun=function(x,...)median(x)) + # add median value
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  theme_ridges() +
  theme(
    legend.position="none"
  ) +
  xlab("") +
  ylab("D index distribution across tissues")

#####################################################
# Compare D index result with EdgeR result

DE.Flower.MTT_Uniq <- read.table("OutputData/DE_genes_of_MTT_unique_in_Flower.txt", sep = "\t" , header=TRUE, stringsAsFactors = F)

F_MTT_Keep <- rgF$Gene.name %in% DE.Flower.MTT_Uniq$Gene.name
F_MTT <- rgF[F_MTT_Keep,]
F_MTT$unique <- "MTT_unique"
rgF$unique <- "Overall"
rg_Flower <- rbind(rgF, F_MTT)

L_MTT_Keep <- rgL$Gene.name %in% DE.Leaf.MTT_Uniq$Gene.name
L_MTT <- rgL[L_MTT_Keep,]
L_MTT$unique <- "MTT_unique"
rgL$unique <- "Overall"
rg_Leaf <- rbind(rgL, L_MTT)

R_MTT_Keep <- rgR$Gene.name %in% DE.Root.MTT_Uniq$Gene.name
R_MTT <- rgR[R_MTT_Keep,]
R_MTT$unique <- "MTT_unique"
rgR$unique <- "Overall"
rg_Root <- rbind(rgR, R_MTT)
dim(R_MTT)

Ridgeline_FLR <- rbind(rg_Flower, rg_Leaf, rg_Root)
head(Ridgeline_FLR)
str(Ridgeline_FLR)
Ridgeline_FLR$Tissue <- factor(Ridgeline_FLR$Tissue, levels = c("Root","Leaf","Flower"))
Ridgeline_FLR$unique <- factor(Ridgeline_FLR$unique, levels = c("Overall","MTT_unique"))

library(dplyr)
library(forcats)


ggplot(Ridgeline_FLR, aes(x=D_index, y=Tissue,  fill=Tissue)) +
  geom_density_ridges(alpha=0.9,
                      bandwidth = c(0.05, 0.03),
                      quantile_lines=TRUE,
                      scale = 1.2,
                      quantile_fun=function(x,...)median(x)) + # add median value
  scale_fill_manual(values = c( "#E7B800", "#00AFBB","#FC4E07")) +
 # scale_fill_viridis(discrete=T) +
 # scale_color_viridis(discrete=F) +
  scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +
  facet_wrap(~unique) +
  theme_ridges() +
  theme(
    legend.position="none"
  ) +
  xlab("") +
  ylab("D index distribution across tissues")


ggplot(Ridgeline_FLR, aes(y = Tissue)) +
  geom_density_ridges(
    aes(x = D_index, fill = paste(Tissue, unique)), 
    scale = 1.2,
    alpha = .8, color = "white", from = -1, to = 1
  ) +
  labs(
    x = "Vote (%)",
    y = "Election Year",
    title = "Indy vs Unionist vote in Catalan elections",
    subtitle = "Analysis unit: municipalities (n = 949)",
    caption = "Marc Belzunces (@marcbeldata) | Source: Idescat"
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_cyclical(
    breaks = c("1980 Indy", "1980 Unionist"),
    labels = c(`1980 Indy` = "Indy", `1980 Unionist` = "Unionist"),
    values = c("#ff0000", "#0000ff", "#ff8080", "#8080ff"),
    name = "Option", guide = "legend"
  ) +
  coord_cartesian(clip = "off") +
  theme_ridges(grid = FALSE)










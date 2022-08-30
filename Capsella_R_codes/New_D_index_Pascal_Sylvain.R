###############################################################
#        Expression index proposed by Pascal & Sylvain        #
###############################################################

library(stringr)
library(matrixStats)
library(ggplot2)
library(cowplot)

## Set work directory and load files
setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project/")
for (tissue in c("F", "L", "R")) {
  Reads <- read.table(paste("InputData/TMM_Diploids_", tissue, ".txt", sep = ""), header = T, sep = "\t")
  head(Reads)
  CGColumns <- str_detect(names(Reads), "CG") # detect which columns contain the "CG" 
  CRColumns <- str_detect(names(Reads), "CR") # detect which columns contain the "CG" 
  COColumns <- str_detect(names(Reads), "CO") # detect which columns contain the "CG" 
  
  Reads$MeanCG <- rowMeans(Reads[, CGColumns], na.rm = T)
  Reads$MeanCR <- rowMeans(Reads[, CRColumns], na.rm = T)
  Reads$MeanCO <- rowMeans(Reads[, COColumns], na.rm = T)
  
  Reads$D.RO <- abs(Reads$MeanCR - Reads$MeanCO) # Distance from CO to CR
  Reads$D.RG <- abs(Reads$MeanCR - Reads$MeanCG) # Distance from CO to CG
  Reads$D.OR <- abs(Reads$MeanCO - Reads$MeanCR) # Distance from CO to CR
  Reads$D.OG <- abs(Reads$MeanCO - Reads$MeanCG) # Distance from CO to CG
  Reads$D.R_Max <- rowMaxs(as.matrix(Reads[,c("D.RO","D.RG")]), na.rm = TRUE)
  Reads$D.O_Max <- rowMaxs(as.matrix(Reads[,c("D.OR","D.OG")]), na.rm = TRUE)
  
  Reads$D.R <- (Reads$D.RO - Reads$D.RG) / Reads$D.R_Max
  Reads$D.O <- (Reads$D.OR - Reads$D.OG) / Reads$D.O_Max
  
  assign(paste0("Reads_", tissue), Reads)
}


head(Reads_F)
output_F <- Reads_F[,c("D.R", "D.O")]
names(output_F) <- c("D_CR", "D_CO")
output_F$Tissue <- "Flower"
write.table(output_F, file = "OutputData/Flower_D_indcies_all_genes.txt",  quote = F, sep = "\t")

output_L <- Reads_L[,c("D.R", "D.O")]
names(output_L) <- c("D_CR", "D_CO")
output_L$Tissue <- "Leaf"
write.table(output_L, file = "OutputData/Leaf_D_indcies_all_genes.txt",  quote = F, sep = "\t")

output_R <- Reads_R[,c("D.R", "D.O")]
names(output_R) <- c("D_CR", "D_CO")
output_R$Tissue <- "Root"
write.table(output_R, file = "OutputData/Root_D_indcies_all_genes.txt",  quote = F, sep = "\t")

dim(Reads_F)

D.R_F <- Reads_F$D.R
D.R_LR <- c(Reads_L$D.R, Reads_R$D.R)

ks.test(D.R_F, D.R_LR)
# set function for notin
'%!in%' <- function(x,y)!('%in%'(x,y))
################################### 1
tissue<-c("F","L","R") # set tissue
# Loop start
for(i in tissue){ # loop 1 start
  print(i)
  # set a sub loop for viable name
  if (i == "F"){ # loop 2 start
    Reads <- Reads_F
  }else{
    if (i == "L"){
      Reads <- Reads_L
    }else{
      Reads <- Reads_R
    }
  } # loop 2 end
  
  for (index in c("D.R", "D.O")) {
    if(index == "D.R"){
      for (com in c("All", "CGCR", "CRCO", "CGCO")) { # loop 3 start
        if (com == "All"){ # if_else start
          Density <- with(density(Reads$D.R, from = -1, to = 1), data.frame(x, y))
          N <- dim(Reads)[1]
          Pos <- round(median(Reads[Reads$D.R > 0,]$D.R), 2)
          Neg <- round(median(Reads[Reads$D.R < 0,]$D.R), 2)
          perPos <- dim(Reads[Reads$D.R > 0,])[1] /  dim(Reads)[1]
          perNeg <- dim(Reads[Reads$D.R < 0,])[1] /  dim(Reads)[1]
          de <- "all genes"
          MyMedian <- round(median(Reads$D.R), 2)
          # write.table(Reads, file = paste("OutputData/DE_", i, "_", com, "_",index,"_index.txt", sep = ""), quote = F, sep = "\t")
        }else{
          DE <- read.csv(paste("OutputData/DE_", i, "_", com, "_FDR.05.txt", sep = ""), header = T, row.names = 1, sep = "\t")
          keep <- row.names(Reads) %in% row.names(DE)
          subReads <- Reads[keep,] 
          N <- dim(subReads)[1]
          assign(paste0("N_",com,"_",i), N)
          Density <- with(density(subReads$D.R, from = -1, to = 1), data.frame(x, y))
          Pos <- round(median(subReads[subReads$D.R > 0,]$D.R), 2)
          Neg <- round(median(subReads[subReads$D.R < 0,]$D.R), 2)
          perPos <- dim(subReads[subReads$D.R > 0,])[1] /  dim(subReads)[1]
          perNeg <- dim(subReads[subReads$D.R < 0,])[1] /  dim(subReads)[1]
          de <- "DE genes"
          MyMedian <- round(median(subReads$D.R), 2)
          #write.table(subReads, file = paste("OutputData/DE_", i, "_", com, "_FDR.05_",index,"_index.txt", sep = ""), quote = F, sep = "\t")
        } # if_else end
        
        assign(paste0(index, "_",i, "_", com, "_Median"), MyMedian)
        assign(paste0(index, "_",i, "_", com, "_Pos"), Pos)
        assign(paste0(index, "_",i, "_", com, "_Neg"), Neg)
        assign(paste0(index, "_",i, "_", com, "_perPos"), perPos)
        assign(paste0(index, "_",i, "_", com, "_perNeg"), perNeg)
        assign(paste0(index, "_",i, "_", com, "_Density"), Density)
      } # loop 3 end
    }else{
      for (com in c("All", "CGCR", "CRCO", "CGCO")) { # loop 3 start
        if (com == "All"){ # if_else start
          Density <- with(density(Reads$D.O, from = -1, to = 1), data.frame(x, y))
          N <- dim(Reads)[1]
          Pos <- round(median(Reads[Reads$D.O > 0,]$D.O), 2)
          Neg <- round(median(Reads[Reads$D.O < 0,]$D.O), 2)
          perPos <- dim(Reads[Reads$D.O > 0,])[1] /  dim(Reads)[1]
          perNeg <- dim(Reads[Reads$D.O < 0,])[1] /  dim(Reads)[1]
          de <- "all genes"
          MyMedian <- round(median(Reads$D.O), 2)
          # write.table(Reads, file = paste("OutputData/DE_", i, "_", com, "_",index,"_index.txt", sep = ""), quote = F, sep = "\t")
        }else{
          DE <- read.csv(paste("OutputData/DE_", i, "_", com, "_FDR.05.txt", sep = ""), header = T, row.names = 1, sep = "\t")
          keep <- row.names(Reads) %in% row.names(DE)
          subReads <- Reads[keep,] 
          N <- dim(subReads)[1]
          assign(paste0("N_",com,"_",i), N)
          Density <- with(density(subReads$D.O, from = -1, to = 1), data.frame(x, y))
          Pos <- round(median(subReads[subReads$D.O > 0,]$D.O), 2)
          Neg <- round(median(subReads[subReads$D.O < 0,]$D.O), 2)
          perPos <- dim(subReads[subReads$D.O > 0,])[1] /  dim(subReads)[1]
          perNeg <- dim(subReads[subReads$D.O < 0,])[1] /  dim(subReads)[1]
          de <- "DE genes"
          MyMedian <- round(median(subReads$D.O), 2)
          #write.table(subReads, file = paste("OutputData/DE_", i, "_", com, "_FDR.05_",index,"_index.txt", sep = ""), quote = F, sep = "\t")
        } # if_else end
        
        assign(paste0(index, "_",i, "_", com, "_Median"), MyMedian)
        assign(paste0(index, "_",i, "_", com, "_Pos"), Pos)
        assign(paste0(index, "_",i, "_", com, "_Neg"), Neg)
        assign(paste0(index, "_",i, "_", com, "_perPos"), perPos)
        assign(paste0(index, "_",i, "_", com, "_perNeg"), perNeg)
        assign(paste0(index, "_",i, "_", com, "_Density"), Density)
      } # loop 3 end
    }
  }
  
} # loop 1 end

D.R_Density <- cbind(D.R_F_All_Density, D.R_L_All_Density$y, D.R_R_All_Density$y)
names(D.R_Density) <- c("D.R","Flowers","Leaves","Roots")
head(D.R_Density)

D.R_CGCRDensity <- cbind(D.R_F_CGCR_Density, D.R_L_CGCR_Density$y, D.R_R_CGCR_Density$y)
names(D.R_CGCRDensity) <- c("D.R","Flowers","Leaves","Roots")
head(D.R_CGCRDensity)

D.R_CRCODensity <- cbind(D.R_F_CRCO_Density, D.R_L_CRCO_Density$y, D.R_R_CRCO_Density$y)
names(D.R_CRCODensity) <- c("D.R","Flowers","Leaves","Roots")
head(D.R_CRCODensity)

D.R_CGCODensity <- cbind(D.R_F_CGCO_Density, D.R_L_CGCO_Density$y, D.R_R_CGCO_Density$y)
names(D.R_CGCODensity) <- c("D.R","Flowers","Leaves","Roots")
head(D.R_CGCODensity)
###
D.O_Density <- cbind(D.O_F_All_Density, D.O_L_All_Density$y, D.O_R_All_Density$y)
names(D.O_Density) <- c("D.O","Flowers","Leaves","Roots")
head(D.O_Density)

D.O_CGCRDensity <- cbind(D.O_F_CGCR_Density, D.O_L_CGCR_Density$y, D.O_R_CGCR_Density$y)
names(D.O_CGCRDensity) <- c("D.O","Flowers","Leaves","Roots")
head(D.O_CGCRDensity)

D.O_CRCODensity <- cbind(D.O_F_CRCO_Density, D.O_L_CRCO_Density$y, D.O_R_CRCO_Density$y)
names(D.O_CRCODensity) <- c("D.O","Flowers","Leaves","Roots")
head(D.O_CRCODensity)

D.O_CGCODensity <- cbind(D.O_F_CGCO_Density, D.O_L_CGCO_Density$y, D.O_R_CGCO_Density$y)
names(D.O_CGCODensity) <- c("D.O","Flowers","Leaves","Roots")
head(D.O_CGCODensity)

ylim <- 1.8
#### Pattern 2

D.overlap.p1 <- -0.264

d1 <- ggplot(data = D.R_Density) +
  geom_area(mapping = aes(x = ifelse(D.R < 0, D.R, NA), y = ylim -0.1), fill = "#f4f4f2", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D.R > 0, D.R, NA), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.R > 0, D.R, NA), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.R > 0, D.R, NA), y = Flowers), fill = "white", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D.R < 0, D.R, NA), y = Flowers), fill = "#fecc47", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D.R < 0, D.R, NA), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.R < 0, D.R, NA), y = Leaves), fill = "#f4f4f2") + 
  geom_line(mapping = aes(x = D.R, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D.R, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D.R, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = D.R_F_All_Median, y = 0, xend = D.R_F_All_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.R_L_All_Median, y = 0, xend = D.R_L_All_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.R_R_All_Median, y = 0, xend = D.R_R_All_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=17307)"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=17307)"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=17307)"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "All genes",
       x = expression(D[italic(CR)]),
       y = "Density"
  ) 

D.overlap <- -0.299
d2 <- ggplot(data = D.R_CGCRDensity) +
  geom_area(mapping = aes(x = ifelse(D.R < 0, D.R, NA), y = ylim -0.1), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.R > D.overlap, D.R, NA), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.R > D.overlap, D.R, NA), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.R > 0, D.R, NA), y = Flowers), fill = "white") + 
  geom_area(mapping = aes(x = ifelse(D.R < D.overlap, D.R, NA), y = Flowers), fill = "#fecc47") + 
  geom_area(mapping = aes(x = ifelse(D.R < D.overlap, D.R, NA), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.R < D.overlap, D.R, NA), y = Leaves), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.R >= D.overlap & D.R <= 0, D.R, NA), y = Flowers), fill = "#f4f4f2") + 
  geom_line(mapping = aes(x = D.R, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D.R, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D.R, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = D.R_F_CGCR_Median, y = 0, xend = D.R_F_CGCR_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.R_L_CGCR_Median, y = 0, xend = D.R_L_CGCR_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.R_R_CGCR_Median, y = 0, xend = D.R_R_CGCR_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=", N_CGCR_F, ")"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=", N_CGCR_L, ")"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=", N_CGCR_R, ")"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "CGCR DE genes",
       x = expression(D[italic(CR)]),
       y = "Density"
  ) 

d3 <- ggplot(data = D.R_CRCODensity) +
  geom_area(mapping = aes(x = ifelse(D.R < 0, D.R, NA), y = ylim -0.1), fill = "#f4f4f2", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D.R > 0, D.R, NA), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.R > 0, D.R, NA), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.R > 0, D.R, NA), y = Flowers), fill = "white", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D.R < 0, D.R, NA), y = Flowers), fill = "#fecc47", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D.R < 0, D.R, NA), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.R < 0, D.R, NA), y = Leaves), fill = "#f4f4f2") + 
  geom_line(mapping = aes(x = D.R, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D.R, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D.R, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = D.R_F_CRCO_Median, y = 0, xend = D.R_F_CRCO_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.R_L_CRCO_Median, y = 0, xend = D.R_L_CRCO_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.R_R_CRCO_Median, y = 0, xend = D.R_R_CRCO_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=", N_CRCO_F, ")"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=", N_CRCO_L, ")"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=", N_CRCO_R, ")"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "CRCO DE genes",
       x = expression(D[italic(CR)]),
       y = "Density"
  ) 

CGCO.overlap <- 0.135
d4 <- ggplot(data = D.R_CGCODensity) +
  geom_area(mapping = aes(x = ifelse(D.R < 0, D.R, NA), y = ylim -0.1), fill = "#f4f4f2", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D.R > 0, D.R, NA), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.R > 0, D.R, NA), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.R > 0, D.R, NA), y = Flowers), fill = "white", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D.R < CGCO.overlap, D.R, NA), y = Flowers), fill = "#fecc47", size = 2) + 
  geom_area(mapping = aes(x = ifelse(D.R <= 0, D.R, NA), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.R <= 0, D.R, NA), y = Leaves), fill = "#f4f4f2") + 
  #geom_area(mapping = aes(x = ifelse(D >= 0 & D <= CGCO.overlap, D, NA), y = Flowers), fill = "#fecc47") + 
  geom_area(mapping = aes(x = ifelse(D.R >= -0.002 & D.R <= CGCO.overlap, D.R, NA), y = Leaves), fill = "white") + 
  geom_line(mapping = aes(x = D.R, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D.R, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D.R, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = D.R_F_CGCO_Median, y = 0, xend = D.R_F_CGCO_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.R_L_CGCO_Median, y = 0, xend = D.R_L_CGCO_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.R_R_CGCO_Median, y = 0, xend = D.R_R_CGCO_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=", N_CGCO_F, ")"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=", N_CGCO_L, ")"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=", N_CGCO_R, ")"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "CGCO DE genes",
       x = expression(D[italic(CR)]),
       y = "Density"
  ) 

newd1 <- ggplot(data = D.O_Density) +
  geom_area(mapping = aes(x = ifelse(D.O < 0, D.O, NA), y = ylim -0.1), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p1, D.O, NA), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p1, D.O, NA), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p1, D.O, NA), y = Flowers), fill = "white") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p1, D.O, NA), y = Flowers), fill = "#fecc47") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p1, D.O, NA), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p1, D.O, NA), y = Leaves), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O >= D.overlap.p1 & D.O <= 0, D.O, NA), y = Flowers), fill = "#f4f4f2") + 
  geom_line(mapping = aes(x = D.O, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D.O, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D.O, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = D.O_F_All_Median, y = 0, xend = D.O_F_All_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.O_L_All_Median, y = 0, xend = D.O_L_All_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.O_R_All_Median, y = 0, xend = D.O_R_All_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=", 17307, ")"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=", 17307, ")"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=", 17307, ")"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "All genes",
       x = expression(D[italic(CO)]),
       y = "Density"
  ) 

D.overlap.p2 <- -0.143
newd2 <- ggplot(data = D.O_CGCRDensity) +
  geom_area(mapping = aes(x = ifelse(D.O < 0, D.O, NA), y = ylim -0.1), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p2, D.O, NA), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p2, D.O, NA), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p2, D.O, NA), y = Flowers), fill = "white") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p2, D.O, NA), y = Flowers), fill = "#fecc47") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p2, D.O, NA), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p2, D.O, NA), y = Leaves), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O >= D.overlap.p2 & D.O <= 0, D.O, NA), y = Flowers), fill = "#f4f4f2") + 
  geom_line(mapping = aes(x = D.O, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D.O, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D.O, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = D.O_F_CGCR_Median, y = 0, xend = D.O_F_CGCR_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.O_L_CGCR_Median, y = 0, xend = D.O_L_CGCR_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.O_R_CGCR_Median, y = 0, xend = D.O_R_CGCR_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=", N_CGCR_F, ")"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=", N_CGCR_L, ")"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=", N_CGCR_R, ")"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "CGCR DE genes",
       x = expression(D[italic(CO)]),
       y = "Density"
  ) 


D.overlap.p3 <- -0.050
newd3 <- ggplot(data = D.O_CRCODensity) +
  geom_area(mapping = aes(x = ifelse(D.O < 0, D.O, NA), y = ylim -0.1), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p3, D.O, NA), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p3, D.O, NA), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p3, D.O, NA), y = Flowers), fill = "white") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p3, D.O, NA), y = Flowers), fill = "#fecc47") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p3, D.O, NA), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p3, D.O, NA), y = Leaves), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O >= D.overlap.p3 & D.O <= 0, D.O, NA), y = Flowers), fill = "#f4f4f2") + 
  geom_line(mapping = aes(x = D.O, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D.O, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D.O, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = D.O_F_CRCO_Median, y = 0, xend = D.O_F_CRCO_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.O_L_CRCO_Median, y = 0, xend = D.O_L_CRCO_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.O_R_CRCO_Median, y = 0, xend = D.O_R_CRCO_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=", N_CRCO_F, ")"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=", N_CRCO_L, ")"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=", N_CRCO_R, ")"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "CRCO DE genes",
       x = expression(D[italic(CO)]),
       y = "Density"
  ) 


D.overlap.p4 <- -0.256
newd4 <- ggplot(data = D.O_CGCODensity) +
  geom_area(mapping = aes(x = ifelse(D.O < 0, D.O, NA), y = ylim -0.1), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p4, D.O, NA), y = Roots), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p4, D.O, NA), y = Leaves), fill = "#a4dbe4") + 
  geom_area(mapping = aes(x = ifelse(D.O > D.overlap.p4, D.O, NA), y = Flowers), fill = "white") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p4, D.O, NA), y = Flowers), fill = "#fecc47") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p4, D.O, NA), y = Roots), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O < D.overlap.p4, D.O, NA), y = Leaves), fill = "#f4f4f2") + 
  geom_area(mapping = aes(x = ifelse(D.O >= D.overlap.p4 & D.O <= 0, D.O, NA), y = Flowers), fill = "#f4f4f2") + 
  geom_line(mapping = aes(x = D.O, y = Leaves), color = "#155800", size=1 )+
  geom_line(mapping = aes(x = D.O, y = Roots), color = "#089400", size=1 ) +
  geom_line(mapping = aes(x = D.O, y = Flowers), color = "#fa990e",size=1 )+
  geom_segment(aes(x = D.O_F_CGCO_Median, y = 0, xend = D.O_F_CGCO_Median, yend = ylim - 0.1), color = "#fa990e", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.O_L_CGCO_Median, y = 0, xend = D.O_L_CGCO_Median, yend = ylim - 0.1), color = "#155800", linetype = "dashed") + # Flowers
  geom_segment(aes(x = D.O_R_CGCO_Median, y = 0, xend = D.O_R_CGCO_Median, yend = ylim - 0.1), color = "#089400", linetype = "dashed") + # Flowers
  annotate(geom="text", x=-1, y=ylim, label= paste0("F: (n=", N_CGCO_F, ")"), color = "#fa990e", hjust = 0) +
  annotate(geom="text", x=0,  y=ylim, label= paste0("L: (n=", N_CGCO_L, ")"), color = "#155800") +
  annotate(geom="text", x=1,  y=ylim, label= paste0("R: (n=", N_CGCO_R, ")"), color = "#089400", hjust = 1) +
  xlim(-1,1) + 
  ylim(0,ylim) +
  theme_classic() +
  labs(title = "CGCO DE genes",
       x = expression(D[italic(CO)]),
       y = "Density"
  ) 

library(ggpubr)
library("grid")
library(ggthemes) # Load
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,1)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(newd1, vp = vplayout(1,1))
print(newd2, vp = vplayout(2,1))
print(newd3, vp = vplayout(3,1))
print(newd4, vp = vplayout(4,1))


p_b1 <- ggplot(OverallOverlap) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "#fde6b6")
p_b2 <- ggplot(OverallOverlap) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "#d3e0ce")
p_b3 <- ggplot(OverallOverlap) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "#cbd7b5")

ggarrange(d1, newd1, d4,  newd2, 
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2,2)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(d1, vp = vplayout(1,2))
print(d4, vp = vplayout(2,2))
print(newd1, vp = vplayout(1,1))
print(newd2, vp = vplayout(2,1))
#print(ggdraw(p_b1), vp = vplayout(1,3:4))
#print(ggdraw(p_b2), vp = vplayout(2,3:4))
#print(ggdraw(p_b3), vp = vplayout(3,3:4))
#print(p1, vp = vplayout(1,3))
#print(p2, vp = vplayout(1,4))
#print(p3, vp = vplayout(2,3))
#print(p4, vp = vplayout(2,4))
#print(p5, vp = vplayout(3,3))
#print(p6, vp = vplayout(3,4))
print(p7, vp = vplayout(3,1:2))
#print(p8, vp = vplayout(3,2))
#print(p10, vp = vplayout(3,2))
#print(p10, vp = vplayout(4,2))


# Contingency test
N.mst <- c(786, 85, 215)
N.cgco <- c(412, 162, 468)
myTable <- cbind(N.mst, N.cgco)
row.names(myTable) <- c("F","L","R")
table(myTable)


library(MASS)
Cars93
Contingency.f1 <- c(rep("case1",307), rep("overlap", 786), rep("case2", 1419),rep("case1",307), rep("overlap", 412), rep("case2", 886))
Contingency.f2 <- c(rep("mst",2512), rep("cgco", 1605))
table(Contingency.f1,Contingency.f2)
fisher.test(Contingency.f1,Contingency.f2)

Contingency.l1 <- c(rep("case1",62), rep("overlap", 86), rep("case2", 562),rep("case1",62), rep("overlap", 162), rep("case2", 951))
Contingency.l2 <- c(rep("mst",710), rep("cgco", 1175))
table(Contingency.l1,Contingency.l2)
fisher.test(Contingency.l1,Contingency.l2)

Contingency.r1 <- c(rep("case1",153), rep("overlap", 215), rep("case2", 720),rep("case1",153), rep("overlap", 468), rep("case2", 2202))
Contingency.r2 <- c(rep("mst",1088), rep("cgco", 2823))
table(Contingency.r1,Contingency.r2)
fisher.test(Contingency.r1,Contingency.r2)

Contingency.t1 <- c(rep("case1",621), rep("overlap", 215), rep("case2", 2965),rep("case1",368), rep("overlap", 468), rep("case2", 4447))
Contingency.t2 <- c(rep("mst",3801), rep("cgco", 5283))
table(Contingency.t1,Contingency.t2)
fisher.test(Contingency.t1,Contingency.t2)

Contingency.f1 <- c(rep("case1",2000), rep("overlap", 786), rep("case1",2000), rep("overlap", 412))
Contingency.f2 <- c(rep("mst",2512), rep("cgco", 1605))
table(Contingency.f1,Contingency.f2)
fisher.test(Contingency.f1,Contingency.f2)


Flower.mst <- c(786, (719+3217))
res <- chisq.test(Flower.mst, p = c( 412/4189, (1093+2684)/4189))
res

Leaf.mst <- c(86, (224 + 1729))
res <- chisq.test(Leaf.mst, p = c(162/2428, 2266/2428))
res

Root.mst <- c(215, 3990)
res <- chisq.test(Root.mst, p = c(468/5472, 5004/5472))
res


Contingency.f1 <- c(rep("case1",4211), rep("overlap", 786), rep("case1", 4052), rep("overlap", 412))
Contingency.f2 <- c(rep("mst",(4211 + 786) ), rep("cgco", (4052+412) ))
table(Contingency.f1,Contingency.f2)
fisher.test(Contingency.f1,Contingency.f2)

Contingency.l1 <- c(rep("case1",1992), rep("overlap", 86), rep("case1",2305), rep("overlap", 162))
Contingency.l2 <- c(rep("mst",(1992+86)), rep("cgco", (2305+162)))
table(Contingency.l1,Contingency.l2)
fisher.test(Contingency.l1,Contingency.l2)

Contingency.r1 <- c(rep("case1",153), rep("overlap", 215), rep("case2", 720),rep("case1",153), rep("overlap", 468), rep("case2", 2202))
Contingency.r2 <- c(rep("mst",1088), rep("cgco", 2823))
table(Contingency.r1,Contingency.r2)
fisher.test(Contingency.r1,Contingency.r2)






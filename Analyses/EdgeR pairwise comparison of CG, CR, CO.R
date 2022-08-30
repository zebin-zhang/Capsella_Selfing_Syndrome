
########    EdgeR Pairwise comparison of CG, CR, CO    ########

setwd("~/Desktop/My_Computer/Uppsala_University/Capsella_Project/")

# Install requested libraries from Bioconductor
# try http:// if https:// URLs are not supported
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# BiocManager::install("edgeR")
# BiocManager::install("limma")
# # Load required libraries 
library(DESeq2)
library(edgeR)
library(limma)
library(Biobase)

scaleRow <- function(x) {
  rm <- rowMeans(x, na.rm = T)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd, na.rm = T)
  x <- sweep(x, 1, sx, "/")
  return(round(x, 6))
}

#### This file contains reads counts infromation of each gene in all samples
# Phased data
raw.RC <- read.table("InputData/Cbp_diploids_scaled_FRL_masked_RNA_ASE_rename.txt", header=T, row.names="gene")
head(raw.RC)
dim(raw.RC)
# CPMnorm <- read.table("Capsella_ReadsCount_FRL_masked_RNA_ASE_FILTERED_CPM.csv", header = T, sep = "\t")
# head(CPMnorm)
# dim(CPMnorm)
##### This file contains information of population and tissue in each sample
Accession <- read.csv("InputData/28genomes_annot_R.csv", sep = ";")
head(Accession)

### Define populations
######################################################
CR_F <- c(Accession$populations=="CR" & Accession$tissue=="flower")
CG_F <- c(Accession$populations=="CG" & Accession$tissue=="flower")
CO_F <- c(Accession$populations=="CO" & Accession$tissue=="flower")

CR_L <- c(Accession$populations=="CR" & Accession$tissue=="leaf")
CG_L <- c(Accession$populations=="CG" & Accession$tissue=="leaf")
CO_L <- c(Accession$populations=="CO" & Accession$tissue=="leaf")

CR_R <- c(Accession$populations=="CR" & Accession$tissue=="root")
CG_R <- c(Accession$populations=="CG" & Accession$tissue=="root")
CO_R <- c(Accession$populations=="CO" & Accession$tissue=="root")

Cbp_Co <- c("ASI_Co", "EUR_Co", "ME_Co", "CASI_Co")
Cbp_Cg <- c("ASI_Cg", "EUR_Cg", "ME_Cg", "CASI_Cg")

Cbp_Co_F <- c(Accession$populations %in% Cbp_Co & Accession$tissue=="flower")
Cbp_Co_L <- c(Accession$populations %in% Cbp_Co & Accession$tissue=="leaf")
Cbp_Co_R <- c(Accession$populations %in% Cbp_Co & Accession$tissue=="root")
Cbp_Cg_F <- c(Accession$populations %in% Cbp_Cg & Accession$tissue=="flower")
Cbp_Cg_L <- c(Accession$populations %in% Cbp_Cg & Accession$tissue=="leaf")
Cbp_Cg_R <- c(Accession$populations %in% Cbp_Cg & Accession$tissue=="root")

ASI_Co_F <- c(Accession$populations=="ASI_Co" & Accession$tissue=="flower")
EUR_Co_F <- c(Accession$populations=="EUR_Co" & Accession$tissue=="flower")
ME_Co_F <- c(Accession$populations=="ME_Co" & Accession$tissue=="flower")
CASI_Co_F <- c(Accession$populations=="CASI_Co" & Accession$tissue=="flower")

ASI_Cg_F <- c(Accession$populations=="ASI_Cg" & Accession$tissue=="flower")
EUR_Cg_F <- c(Accession$populations=="EUR_Cg" & Accession$tissue=="flower")
ME_Cg_F <- c(Accession$populations=="ME_Cg" & Accession$tissue=="flower")
CASI_Cg_F <- c(Accession$populations=="CASI_Cg" & Accession$tissue=="flower")

ASI_Co_L <- c(Accession$populations=="ASI_Co" & Accession$tissue=="leaf")
EUR_Co_L <- c(Accession$populations=="EUR_Co" & Accession$tissue=="leaf")
ME_Co_L <- c(Accession$populations=="ME_Co" & Accession$tissue=="leaf")
CASI_Co_L <- c(Accession$populations=="CASI_Co" & Accession$tissue=="leaf")

ASI_Cg_L <- c(Accession$populations=="ASI_Cg" & Accession$tissue=="leaf")
EUR_Cg_L <- c(Accession$populations=="EUR_Cg" & Accession$tissue=="leaf")
ME_Cg_L <- c(Accession$populations=="ME_Cg" & Accession$tissue=="leaf")
CASI_Cg_L <- c(Accession$populations=="CASI_Cg" & Accession$tissue=="leaf")

ASI_Co_R <- c(Accession$populations=="ASI_Co" & Accession$tissue=="root")
EUR_Co_R <- c(Accession$populations=="EUR_Co" & Accession$tissue=="root")
ME_Co_R <- c(Accession$populations=="ME_Co" & Accession$tissue=="root")
CASI_Co_R <- c(Accession$populations=="CASI_Co" & Accession$tissue=="root")

ASI_Cg_R <- c(Accession$populations=="ASI_Cg" & Accession$tissue=="root")
EUR_Cg_R <- c(Accession$populations=="EUR_Cg" & Accession$tissue=="root")
ME_Cg_R <- c(Accession$populations=="ME_Cg" & Accession$tissue=="root")
CASI_Cg_R <- c(Accession$populations=="CASI_Cg" & Accession$tissue=="root")
######################################################
raw.RC <- raw.RC[,-(37:132)] # Remove data from Cbp
head(raw.RC)
dim(raw.RC)
# Filter reads
geneMis <- 1 
geneKeep <- cbind((rowSums(raw.RC[,CG_F] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CR_F] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CO_F] >= 1, na.rm = T) >= geneMis) &
                    
                    (rowSums(raw.RC[,CG_L] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CR_L] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CO_L] >= 1, na.rm = T) >= geneMis) &
                    
                    (rowSums(raw.RC[,CG_R] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CR_R] >= 1, na.rm = T) >= geneMis) &
                    (rowSums(raw.RC[,CO_R] >= 1, na.rm = T) >= geneMis) 
)

raw.RC.filter1 <- raw.RC[geneKeep,]
raw.RC.remove <- raw.RC[!geneKeep,]
dim(raw.RC.remove)
hea
message(paste("Number of genes before", dim(raw.RC)[1], "\nNumber of genes after", dim(raw.RC.filter1)[1]))
head(RC)

## Filter HSE
# Filter by NA number
# For each gene in each tissue in parental species, at least 5 samples contain none-zero value.
Fkeep <- rowSums(raw.RC.filter1[,c(1:12)] >= 1, na.rm = T) >= 5
Lkeep <- rowSums(raw.RC.filter1[,c(13:24)] >= 1, na.rm = T) >= 5
Rkeep <- rowSums(raw.RC.filter1[,c(25:36)] >= 1, na.rm = T) >= 5

# merge three logical vectors into one column by cbind(## & ##)
# cbind(## , ##) will combine several columns into one dataframe. 
FLRkeep <- cbind(Fkeep & Lkeep & Rkeep) 

RC <- raw.RC.filter1[FLRkeep,] 
dim(RC)
message("Number of genes BEFORE filter2: ", dim(raw.RC.filter1)[1], 
        "\nNumber of genes AFTER filter2: ",  dim(RC)[1])

tissue<-c("F","L","R") # set tissue

# Loop start
for (i in tissue){
  print(i)
  # set a sub loop for viable name
  if (i == "F"){
    subData <- "Flower"
  }else{
    if (i == "L"){
      subData <- "Leaf"
    }else{
      subData <- "Root"
    }
  }
  
  reads <- RC[, which(sub(".*_","",colnames(RC)) == i)]
  reads <- na.omit(reads)
  reads <- reads[,c(1:4, 9:12, 5:8)] # Order columns by CR CO CG
  head(reads)
  myGroup <- sub("[0-9]_.*", "", colnames(reads))
  myGroup <- factor(myGroup, levels = c("CR","CO","CG"))
  message("my group in ", i, " is: ", myGroup)
  
  myDE <- DGEList(counts = reads, group = myGroup)
  # Normalizing the data
  myDE <- calcNormFactors(myDE, method = "TMM")
  # Estimating the Dispersion
  # This is the first major step in the analysis of DGE data
  myDE <- estimateCommonDisp(myDE)
  # Extract TMM values.
  TMM.Keep <- myDE$pseudo.counts
  # Output data
  # write.table(TMM.Keep, paste("InputData/TMM_Diploids_", i, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
  
  # Create the contrast matrix
  myDE.mat <- model.matrix(~ 0 + myDE$samples$group)
  colnames(myDE.mat) <- levels(myDE$samples$group)
  message("levels of my group in ", i, " is: ", levels(myDE$samples$group)) # check the group order.
  
  # Estimate dispersion parameter for GLM
  myDE2 <- estimateGLMCommonDisp(myDE, myDE.mat)
  myDE2 <- estimateGLMTrendedDisp(myDE2, myDE.mat)
  # You can change method to "auto", "bin.spline", "power", "spline", "bin.loess".
  # The default is "auto" which chooses "bin.spline" when > 200 tags and "power" otherwise.
  # I chose the method="power" after looking at the methods that are offered: bin.spline (default if number of 
  # tags is > 200), power (default otherwise), bin.loess, and spline. We have 17,396 tags, so the default is 
  # bin.spline. When I used the bin.spline method, it was no better than estimating a common dispersion so I 
  # instead used power.
  myDE2 <- estimateGLMTagwiseDisp(myDE2, myDE.mat)
  
  ## DEG analysis 
  # Model fitting
  fit <- glmFit(myDE2, myDE.mat)
  
  # Loop 2
  for (com in c("CRCG", "COCG", "CRCO")) {
    
    # compare (group 1 <CR> - group 3 <CG>) to 0: 
    # this is equivalent to comparing group 1 to group 3
    if (com == "CRCG"){
      lrt <- glmLRT(fit, contrast=c(1,0,-1)) # group 1 to group 3, CG to CR
      vs <- "CR/CG"
    }else{
      if (com == "COCG"){
        lrt <- glmLRT(fit, contrast=c(0,1,-1)) # group 2 to group 3, CG to CO
        vs <- "CO/CG"
      }else{
        lrt <- glmLRT(fit, contrast=c(1,-1,0)) # group 1 to group 2, CR to CO
        vs <- "CR/CO"
      }
    }

    # Access results tables
    Table <- lrt$table

    # Toptags Result with FDR
    TopTags <- topTags(lrt, n = nrow(lrt$table))$table

    # Annotate data according to up-regulated and down-regulated
    #### ggscatter ####
    library(ggpubr)
    #install.packages("ggExtra")
    library(ggExtra)
    
    TopTags$logFDR <- -log10(TopTags$FDR)
    
    # Loop for  TopTags$Tissue
    if (i == "F"){
      TopTags$Tissue <- "Flower"
    }else{
      if (i == "L"){
        TopTags$Tissue <- "Leaf"
      }else{
        TopTags$Tissue <- "Root"
      }
    }
    
    TopTags$DE <- "Not DE"
    TopTags[(TopTags$logFC > 0) & (TopTags$FDR <= 0.05),]$DE <- "Up" # ！not set for FC
    TopTags[(TopTags$logFC < 0) & (TopTags$FDR <= 0.05),]$DE <- "Down"
    TopTags$DE <- factor(TopTags$DE, levels = c("Up", "Down", "Not DE"))
    TopTags$Comparison <- vs
    # write data out
    # All DE result
    # write.table(TopTags, paste("OutputData/DE_", i, "_", com, ".txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
    # Extract the DE genes with FDR < 0.05
    # FC > 2
    DE <- TopTags[(abs(TopTags$logFC) > 1),]
    # write.table(DE, paste("OutputData/DE_", i, "_", com, "_FC2.txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
    # FDR ≤ 0.05
    DE <- TopTags[(TopTags$FDR <= 0.05),]
    # write.table(DE, paste("OutputData/DE_", i, "_", com, "_FDR.05.txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
    assign(paste0("DE_",i,"_",com), DE)
    # FC > 2 & FDR ≤ 0.05
    DE <- TopTags[(abs(TopTags$logFC) > 1) & (TopTags$FDR <= 0.05),]
    # write.table(DE, paste("OutputData/DE_", i, "_", com, "_FC2FDR.05.txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
    #assign(paste0("DE_",i,"_",com), DE)
    # FC > 2 & FDR <= 0.01
    DE <- TopTags[(abs(TopTags$logFC) > 1) & (TopTags$FDR <= 0.01),]
    # write.table(DE, paste("OutputData/DE_", i, "_", com, "_FC2FDR.01.txt", sep = ""), quote = F, sep = "\t", row.names = T, col.names = T)
    

    ### Message out 
    count_R <- dim(TopTags[TopTags$DE == "Up",])[1]
    count_L <- dim(TopTags[TopTags$DE == "Down",])[1]
    median_R <- round(median(TopTags[TopTags$DE == "Up",]$logFC), 2)
    median_L <- round(median(TopTags[TopTags$DE == "Down",]$logFC), 2)
    percentage_R <- paste0(round((dim(TopTags[TopTags$DE == "Up",])[1] / dim(TopTags)[1])*100, 2), "%")
    percentage_L <- paste0(round((dim(TopTags[TopTags$DE == "Down",])[1] / dim(TopTags)[1])*100, 2), "%")
    
    message("Numbers of Up regulated genes of ", com, " in ", i, " is: ", dim(TopTags[TopTags$DE == "Up",])[1])
    message("Numbers of Down regulated genes of ", com, " in ", i, " is: ", dim(TopTags[TopTags$DE == "Down",])[1])
    message("Percentage of Up regulated genes of ", com, " in ", i, " is: ", (dim(TopTags[TopTags$DE == "Up",])[1] / dim(TopTags)[1])*100)
    message("Percentage of Down regulated genes of ", com, " in ", i, " is: ", (dim(TopTags[TopTags$DE == "Down",])[1] / dim(TopTags)[1])*100)
    
    ### Second pattern
    library(ggplot2)
    library(cowplot) 
   
    # Main plot
    # Loop for Main plot according Comparisons
    if (com == "CRCO"){
      Main <- ggscatter(TopTags, x = "logFC", y = "logFDR",
                        color = "DE",       # Color by groups "DE"
                        palette = c("#1B9E77", "#1B9E77", "gray"), # reset the colour dark green     
                        size = 0.2, alpha = 0.5,      # change the plot size and transparency
                        title = paste(vs, "in", subData), # Add title
                        ggtheme = theme_bw() # change background to black pattern
      ) +
        xlab(expression(log[2]~FC(CR/CO))) + 
        ylab(expression(-log[10]~(FDR))) + 
        rremove("legend") + # remove figure legend
        scale_x_continuous(limits = c(-15, 15),  # set the x axis limit
                           breaks = get_breaks(by = 5, from = -15)) + # set the breaks in x axis
        scale_y_continuous(limits = c(0, 50),  # set the x axis limit
                           breaks = get_breaks(by = 10, from = 0)) + # set the breaks in x axis
        annotate(geom="text", x=-15, y=40, label=count_L, color = "#1B9E77", fontface =2, hjust=0) + # add annotation for median down
        annotate(geom="text", x=15, y=40, label=count_R, color = "#1B9E77", fontface =2, hjust=1)  # median up
      # Marginal densities along x axis
      # Need to set coord_flip = TRUE, if you plan to use coord_flip()
      Xdens <- axis_canvas(Main, axis = "x") + # marginal will add to main plot along x axis.
        geom_density(data = TopTags, aes(x = logFC, fill = DE),
                     alpha = 0.7, size = 0.2)+
        scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
        ggpubr::fill_palette(c("#1B9E77", "#1B9E77", "gray")) + 
        annotate(geom="text", x=-15, y=0.5, label=percentage_L, color = "#1B9E77", fontface =2, hjust=0) + # add annotation for median down
        annotate(geom="text", x=15, y=0.5, label=percentage_R, color = "#1B9E77", fontface =2, hjust=1) + # add annotation for median down
        rremove("legend")
      
      Xboxs <- ggboxplot(TopTags, x = "DE", y = "logFC",
                                color = "DE", fill = "DE", palette = c("#1B9E77", "#1B9E77", "gray"),
                                ylim = c(-15,15),
                                alpha = 0.5,
                         rotate= TRUE,  # rotate the boxplot by horizontal
                                ggtheme = theme_bw()) +
        scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
        stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") + # add mean value
        annotate(geom="text", x="Down", y=15, label=paste(median_R, ",", percentage_R), color = "#1B9E77", fontface =2, hjust=1) + # add annotation for median down
        annotate(geom="text", x="Up", y=-15, label=paste(median_L, ",", percentage_L), color = "#1B9E77", fontface =2, hjust=0) + # add annotation for median down
        rremove("legend")
      
      #p1 <- insert_xaxis_grob(Main, Xdens, grid::unit(.5, "null"), position = "top")
      p2 <- insert_xaxis_grob(Main, Xboxs, grid::unit(.5, "null"), position = "top")
      ggdraw(p2)
      assign(paste0(i, "_", com, "_p2"), p2)
    }else{
      if (com == "CRCG"){
        Main <- ggscatter(TopTags, x = "logFC", y = "logFDR",
                          color = "DE",       # Color by groups "DE"
                          palette = c("#D95F02", "#D95F02", "gray"), # reset the colour        
                          size = 0.2, alpha = 0.5,      # change the plot size and transparency
                          title = paste(vs, "in", subData), # Add title
                          ggtheme = theme_bw() # change background to black pattern
        ) + 
          xlab(expression(log[2]~FC(CR/CG))) + 
          ylab(expression(-log[10]~(FDR))) + 
          rremove("legend") + # remove figure legend
          scale_x_continuous(limits = c(-15, 15),  # set the x axis limit
                             breaks = get_breaks(by = 5, from = -15)) + # set the breaks in x axis
          scale_y_continuous(limits = c(0, 50),  # set the x axis limit
                             breaks = get_breaks(by = 10, from = 0)) + # set the breaks in x axis
          annotate(geom="text", x=-15, y=40, label=count_L, color = "#D95F02", fontface =2, hjust=0) + # add annotation for median down
          annotate(geom="text", x=15, y=40, label=count_R, color = "#D95F02", fontface =2, hjust=1)  # median up
        # Marginal densities along x axis
        # Need to set coord_flip = TRUE, if you plan to use coord_flip()
        Xdens <- axis_canvas(Main, axis = "x") + # marginal will add to main plot along x axis.
          geom_density(data = TopTags, aes(x = logFC, fill = DE),
                       alpha = 0.7, size = 0.2)+
          scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
          ggpubr::fill_palette(c("#D95F02", "#D95F02", "gray")) + 
          annotate(geom="text", x=-15, y=0.5, label=percentage_L, color = "#D95F02", fontface =2, hjust=0) + # add annotation for median down
          annotate(geom="text", x=15, y=0.5, label=percentage_R, color = "#D95F02", fontface =2, hjust=1) + # add annotation for median down
          rremove("legend")
        
        Xboxs <- ggboxplot(TopTags, x = "DE", y = "logFC",
                           color = "DE", fill = "DE", palette = c("#D95F02", "#D95F02", "gray"),
                           ylim = c(-15,15),
                           alpha = 0.5,
                           rotate= TRUE,  # rotate the boxplot by horizontal
                           ggtheme = theme_bw()) +
          scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
          stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") + # add mean value
          annotate(geom="text", x="Down", y=15, label=paste(median_R, ",", percentage_R), color = "#D95F02", fontface =2, hjust=1) + # add annotation for median down
          annotate(geom="text", x="Up", y=-15, label=paste(median_L, ",", percentage_L), color = "#D95F02", fontface =2, hjust=0) + # add annotation for median down
          rremove("legend")
        
        #p1 <- insert_xaxis_grob(Main, Xdens, grid::unit(.5, "null"), position = "top")
        p2 <- insert_xaxis_grob(Main, Xboxs, grid::unit(.5, "null"), position = "top")
        assign(paste0(i, "_", com, "_p2"), p2)
        ggdraw(p2)
      }else{
        Main <- ggscatter(TopTags, x = "logFC", y = "logFDR",
                          color = "DE",       # Color by groups "DE"
                          palette = c("#7570B3", "#7570B3", "gray"), # reset the colour        
                          size = 0.2, alpha = 0.5,      # change the plot size and transparency
                          title = paste(vs, "in", subData), # Add title
                          ggtheme = theme_bw() # change background to black pattern
        ) + 
          xlab(expression(log[2]~FC(CO/CG))) + 
          ylab(expression(-log[10]~(FDR))) + 
          rremove("legend") + # remove figure legend
          scale_x_continuous(limits = c(-15, 15),  # set the x axis limit
                             breaks = get_breaks(by = 5, from = -15)) + # set the breaks in x axis
          scale_y_continuous(limits = c(0, 50),  # set the x axis limit
                             breaks = get_breaks(by = 10, from = 0)) + # set the breaks in x axis
          annotate(geom="text", x=-15, y=40, label=count_L, color = "#7570B3", fontface =2, hjust=0) + # add annotation for median down
          annotate(geom="text", x=15, y=40, label=count_R, color = "#7570B3", fontface =2, hjust=1)  # median up
        # Marginal densities along x axis
        # Need to set coord_flip = TRUE, if you plan to use coord_flip()
        Xdens <- axis_canvas(Main, axis = "x") + # marginal will add to main plot along x axis.
          geom_density(data = TopTags, aes(x = logFC, fill = DE),
                       alpha = 0.7, size = 0.2)+
          scale_x_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
          ggpubr::fill_palette(c("#7570B3", "#7570B3", "gray")) + 
          annotate(geom="text", x=-15, y=0.5, label=percentage_L, color = "#7570B3", fontface =2, hjust=0) + # add annotation for median down
          annotate(geom="text", x=15, y=0.5, label=percentage_R, color = "#7570B3", fontface =2, hjust=1) + # add annotation for median down
          rremove("legend")
        
        Xboxs <- ggboxplot(TopTags, x = "DE", y = "logFC",
                           color = "DE", fill = "DE", palette = c("#7570B3", "#7570B3", "gray"),
                           ylim = c(-15,15),
                           alpha = 0.5,
                           rotate= TRUE,  # rotate the boxplot by horizontal
                           ggtheme = theme_bw()) +
          scale_y_continuous( breaks = get_breaks(by = 5, from = -15), limits = c(-15, 15)) +
          stat_summary(fun = mean, geom="point", shape=18, size = 2, color="black", fill="black") + # add mean value
          annotate(geom="text", x="Down", y=15, label=paste(median_R, ",", percentage_R), color = "#7570B3", fontface =2, hjust=1) + # add annotation for median down
          annotate(geom="text", x="Up", y=-15, label=paste(median_L, ",", percentage_L), color = "#7570B3", fontface =2, hjust=0) + # add annotation for median down
          rremove("legend")
        
        #p1 <- insert_xaxis_grob(Main, Xdens, grid::unit(.5, "null"), position = "top")
        p2 <- insert_xaxis_grob(Main, Xboxs, grid::unit(.5, "null"), position = "top")
        ggdraw(p2)
        assign(paste0(i, "_", com, "_p2"), p2)
      }
    }
    
  }
  
}

# combine all plot together
ggarrange(F_CRCO_p2, F_CRCG_p2, F_COCG_p2,
          L_CRCO_p2, L_CRCG_p2, L_COCG_p2,
          R_CRCO_p2, R_CRCG_p2, R_COCG_p2,
          labels = c("A","B","C","D","E","F","G","H","I"),
          nrow = 3, ncol = 3)

# Venn plot
# Install the latest development version:
# if (!require(devtools)) install.packages("devtools")
# devtools::install_github("yanlinlin82/ggvenn")
library("ggvenn")

venn_F <- list(
  CRCG = row.names(DE_F_CRCG),
  COCG = row.names(DE_F_COCG),
  CRCO = row.names(DE_F_CRCO)
)

venn_L <- list(
  CRCG = row.names(DE_L_CRCG),
  COCG = row.names(DE_L_COCG),
  CRCO = row.names(DE_L_CRCO)
)

venn_R <- list(
  CRCG = row.names(DE_R_CRCG),
  COCG = row.names(DE_R_COCG),
  CRCO = row.names(DE_R_CRCO)
)

# Change category names
# Change the fill color
names(venn_F) <- c("CR/CG","CO/CG","CR/CO")
names(venn_L) <- c("CR/CG","CO/CG","CR/CO")
names(venn_R) <- c("CR/CG","CO/CG","CR/CO")

venn_F_p <- ggvenn(
  venn_F, 
  fill_color = c("#D95F02", "#7570B3","#1B9E77"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = T
)

venn_L_p <- ggvenn(
  venn_L, 
  fill_color = c("#D95F02", "#7570B3", "#1B9E77"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = T
)

venn_R_p <- ggvenn(
  venn_R, 
  fill_color = c("#D95F02","#7570B3","#1B9E77"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = T
)
library(ggpubr)
library("grid")
library(ggthemes) # Load
grid.newpage()
ggarrange(F_CRCG_p2, F_CRCO_p2, F_COCG_p2, venn_F_p,
          L_CRCG_p2, L_CRCO_p2, L_COCG_p2, venn_L_p,
          R_CRCG_p2, R_CRCO_p2, R_COCG_p2, venn_R_p,
          labels = c("A","","","D","B","","","E","C","","","F"),
          nrow = 3, ncol = 4)

# set function for NotIn
'%!in%' <- function(x,y)!('%in%'(x,y))
# Define MatingSystem related genes
# include in CGCR, CGCO but not in CRCO
F.filter1 <- row.names(DE_F_CRCG) %in% row.names(DE_F_COCG)
F.filter2 <- row.names(DE_F_CRCG) %!in% row.names(DE_F_CRCO)
F.coGO_1 <- row.names(DE_F_CRCG) %in% row.names(DE_F_CRCO)  # convergence between CG CO
F.coGO_2 <- row.names(DE_F_CRCG) %!in% row.names(DE_F_COCG)  
F.coGR_1 <- row.names(DE_F_COCG) %in% row.names(DE_F_CRCO)  # convergence between CG CR
F.coGR_2 <- row.names(DE_F_COCG) %!in% row.names(DE_F_CRCG)  

Fkeep <- cbind(F.filter1 & F.filter2)
Fkeep_coGO <- cbind(F.coGO_1 & F.coGO_2)
Fkeep_coGR <- cbind(F.coGR_1 & F.coGR_2)

F.coCRCO_all <- DE_F_CRCG[F.filter1,]
F.coCGCO_all <- DE_F_CRCG[F.coGO_1,]
F.coCGCR_all <- DE_F_COCG[F.coGR_1,]
F.mst <-  DE_F_CRCG[Fkeep,]
F.coGO <- DE_F_CRCG[Fkeep_coGO,]
F.coGR <- DE_F_COCG[Fkeep_coGR,]
head(F.mst)

L.filter1 <- row.names(DE_L_CRCG) %in% row.names(DE_L_COCG)
L.filter2 <- row.names(DE_L_CRCG) %!in% row.names(DE_L_CRCO)
L.coGO_1 <- row.names(DE_L_CRCG) %in% row.names(DE_L_CRCO)  # convergence between CG CO
L.coGO_2 <- row.names(DE_L_CRCG) %!in% row.names(DE_L_COCG)  
L.coGR_1 <- row.names(DE_L_COCG) %in% row.names(DE_L_CRCO)  # convergence between CG CR
L.coGR_2 <- row.names(DE_L_COCG) %!in% row.names(DE_L_CRCG)  
Lkeep <- cbind(L.filter1 & L.filter2)
Lkeep_coGO <- cbind(L.coGO_1 & L.coGO_2)
Lkeep_coGR <- cbind(L.coGR_1 & L.coGR_2)

L.coCRCO_all <- DE_L_CRCG[L.filter1,]
L.coCGCO_all <- DE_L_CRCG[L.coGO_1,]
L.coCGCR_all <- DE_L_COCG[L.coGR_1,]
L.mst <- DE_L_CRCG[Lkeep,]
L.coGO <- DE_L_CRCG[Lkeep_coGO,]
L.coGR <- DE_L_COCG[Lkeep_coGR,]
dim(L.mst)
head(L.mst)

R.filter1 <- row.names(DE_R_CRCG) %in% row.names(DE_R_COCG)
R.filter2 <- row.names(DE_R_CRCG) %!in% row.names(DE_R_CRCO)
R.coGO_1 <- row.names(DE_R_CRCG) %in% row.names(DE_R_CRCO)  # convergence between CG CO
R.coGO_2 <- row.names(DE_R_CRCG) %!in% row.names(DE_R_COCG)  
R.coGR_1 <- row.names(DE_R_COCG) %in% row.names(DE_R_CRCO)  # convergence between CG CR
R.coGR_2 <- row.names(DE_R_COCG) %!in% row.names(DE_R_CRCG)  
Rkeep <- cbind(R.filter1 & R.filter2)
Rkeep_coGO <- cbind(R.coGO_1 & R.coGO_2)
Rkeep_coGR <- cbind(R.coGR_1 & R.coGR_2)

R.coCRCO_all <- DE_R_CRCG[R.filter1,]
R.coCGCO_all <- DE_R_CRCG[R.coGO_1,]
R.coCGCR_all <- DE_R_COCG[R.coGR_1,]
R.mst <- DE_R_CRCG[Rkeep,]
R.coGO <- DE_R_CRCG[Rkeep_coGO,]
R.coGR <- DE_R_COCG[Rkeep_coGR,]
dim(R.mst)
head(R.mst)

# coDEGs
F.coCRCO_all_1 <- F.coCRCO_all[,c("Tissue", "Comparison", "DE")]
F.coCRCO_all_2 <- DE_F_COCG[row.names(DE_F_COCG) %in% row.names(F.coCRCO_all), c("Comparison", "DE")]
coDEGs_CRCO_F <- merge(F.coCRCO_all_1, F.coCRCO_all_2, by=0, all = T)
names(coDEGs_CRCO_F) <- c("Gene", "Tissue", "Comparison_1", "DE_1","Comparison_2", "DE_2")
coDEGs_CRCO_F$Status <- "CRneCO"
coDEGs_CRCO_F[coDEGs_CRCO_F$Gene %in% row.names(F.mst),]$Status <- "CReqCO"
head(coDEGs_CRCO_F)
dim(coDEGs_CRCO_F)
dim(coDEGs_CRCO_F[coDEGs_CRCO_F$Status == "CReqCO",])


L.coCRCO_all_1 <- L.coCRCO_all[,c("Tissue", "Comparison", "DE")]
L.coCRCO_all_2 <- DE_L_COCG[row.names(DE_L_COCG) %in% row.names(L.coCRCO_all), c("Comparison", "DE")]
coDEGs_CRCO_L <- merge(L.coCRCO_all_1, L.coCRCO_all_2, by=0, all = T)
names(coDEGs_CRCO_L) <- c("Gene", "Tissue", "Comparison_1", "DE_1","Comparison_2", "DE_2")
coDEGs_CRCO_L$Status <- "CRneCO"
coDEGs_CRCO_L[coDEGs_CRCO_L$Gene %in% row.names(L.mst),]$Status <- "CReqCO"
head(coDEGs_CRCO_L)
dim(coDEGs_CRCO_L)
dim(coDEGs_CRCO_L[coDEGs_CRCO_L$Status == "CReqCO",])

R.coCRCO_all_1 <- R.coCRCO_all[,c("Tissue", "Comparison", "DE")]
R.coCRCO_all_2 <- DE_R_COCG[row.names(DE_R_COCG) %in% row.names(R.coCRCO_all), c("Comparison", "DE")]
coDEGs_CRCO_R <- merge(R.coCRCO_all_1, R.coCRCO_all_2, by=0, all = T)
names(coDEGs_CRCO_R) <- c("Gene", "Tissue", "Comparison_1", "DE_1","Comparison_2", "DE_2")
coDEGs_CRCO_R$Status <- "CRneCO"
coDEGs_CRCO_R[coDEGs_CRCO_R$Gene %in% row.names(R.mst),]$Status <- "CReqCO"
head(coDEGs_CRCO_R)
dim(coDEGs_CRCO_R)
dim(coDEGs_CRCO_R[coDEGs_CRCO_R$Status == "CReqCO",])

F.coCGCO_all_1 <- F.coCGCO_all[,c("Tissue", "Comparison", "DE")]
F.coCGCO_all_2 <- DE_F_CRCO[row.names(DE_F_CRCO) %in% row.names(F.coCGCO_all), c("Comparison", "DE")]
coDEGs_CGCO_F <- merge(F.coCGCO_all_1, F.coCGCO_all_2, by=0, all = T)
names(coDEGs_CGCO_F) <- c("Gene", "Tissue", "Comparison_1", "DE_1","Comparison_2", "DE_2")
coDEGs_CGCO_F$Status <- "CGneCO"
coDEGs_CGCO_F[coDEGs_CGCO_F$Gene %in% row.names(F.coGO),]$Status <- "CGeqCO"
head(coDEGs_CGCO_F)
dim(coDEGs_CGCO_F)
dim(coDEGs_CGCO_F[coDEGs_CGCO_F$Status == "CGeqCO",])

L.coCGCO_all_1 <- L.coCGCO_all[,c("Tissue", "Comparison", "DE")]
L.coCGCO_all_2 <- DE_L_CRCO[row.names(DE_L_CRCO) %in% row.names(L.coCGCO_all), c("Comparison", "DE")]
coDEGs_CGCO_L <- merge(L.coCGCO_all_1, L.coCGCO_all_2, by=0, all = T)
names(coDEGs_CGCO_L) <- c("Gene", "Tissue", "Comparison_1", "DE_1","Comparison_2", "DE_2")
coDEGs_CGCO_L$Status <- "CGneCO"
coDEGs_CGCO_L[coDEGs_CGCO_L$Gene %in% row.names(L.coGO),]$Status <- "CGeqCO"
head(coDEGs_CGCO_L)
dim(coDEGs_CGCO_L)
dim(coDEGs_CGCO_L[coDEGs_CGCO_L$Status == "CGeqCO",])

R.coCGCO_all_1 <- R.coCGCO_all[,c("Tissue", "Comparison", "DE")]
R.coCGCO_all_2 <- DE_R_CRCO[row.names(DE_R_CRCO) %in% row.names(R.coCGCO_all), c("Comparison", "DE")]
coDEGs_CGCO_R <- merge(R.coCGCO_all_1, R.coCGCO_all_2, by=0, all = T)
names(coDEGs_CGCO_R) <- c("Gene", "Tissue", "Comparison_1", "DE_1","Comparison_2", "DE_2")
coDEGs_CGCO_R$Status <- "CGneCO"
coDEGs_CGCO_R[coDEGs_CGCO_R$Gene %in% row.names(R.coGO),]$Status <- "CGeqCO"
head(coDEGs_CGCO_R)
dim(coDEGs_CGCO_R)
dim(coDEGs_CGCO_R[coDEGs_CGCO_R$Status == "CGeqCO",])


# write table out
# head(F.mst)
write.table(coDEGs_CRCO_F, file = "OutputData/coDEGs_CRCO_F.txt", quote = F, sep = "\t", row.names = F)
write.table(coDEGs_CRCO_L, file = "OutputData/coDEGs_CRCO_L.txt", quote = F, sep = "\t", row.names = F)
write.table(coDEGs_CRCO_R, file = "OutputData/coDEGs_CRCO_R.txt", quote = F, sep = "\t", row.names = F)
write.table(coDEGs_CGCO_F, file = "OutputData/coDEGs_CGCO_F.txt", quote = F, sep = "\t", row.names = F)
write.table(coDEGs_CGCO_L, file = "OutputData/coDEGs_CGCO_L.txt", quote = F, sep = "\t", row.names = F)
write.table(coDEGs_CGCO_R, file = "OutputData/coDEGs_CGCO_R.txt", quote = F, sep = "\t", row.names = F)
# write.table(F.coGR, file = "OutputData/coDEG_CGCR_F.txt", quote = F, sep = "\t", row.names = T)
# write.table(L.coGR, file = "OutputData/coDEG_CGCR_L.txt", quote = F, sep = "\t", row.names = T)
# write.table(R.coGR, file = "OutputData/coDEG_CGCR_R.txt", quote = F, sep = "\t", row.names = T)



head(coDEGs_CRCO_L)


venn_MST <- list(
  MST_F = row.names(F.mst),
  MST_L = row.names(L.mst),
  MST_R = row.names(R.mst)
)

venn_coGO <- list(
  coCGCO_F = row.names(F.coGO),
  coCGCO_L = row.names(L.coGO),
  coCGCO_R = row.names(R.coGO)
)

venn_coGR <- list(
  coCGCR_F = row.names(F.coGR),
  coCGCR_L = row.names(L.coGR),
  coCGCR_R = row.names(R.coGR)
)

venn_MST_p <- ggvenn(
  venn_MST, 
  fill_color = c("#fa990e", "#155800", "#089400"),
  fill_alpha = 0.8,
  stroke_color = "white",
  stroke_size = 1,
  set_name_size = 4,
  show_percentage = T
)

venn_coGO_p <- ggvenn(
  venn_coGO, 
  fill_color = c("#fa990e", "#155800", "#089400"),
  fill_alpha = 0.8,
  stroke_color = "white",
  stroke_size = 1,
  set_name_size = 4,
  show_percentage = T
)

venn_coGR_p <- ggvenn(
  venn_coGR, 
  fill_color = c("#fa990e", "#155800", "#089400"),
  fill_alpha = 0.8,
  stroke_color = "white",
  stroke_size = 1,
  set_name_size = 4,
  show_percentage = T
)

library(grDevices)
# drawing background colour
p_b1 <- ggplot(F.mst) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "#fde6b6")
p_b2 <- ggplot(F.mst) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "#d3e0ce")
p_b3 <- ggplot(F.mst) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "#cbd7b5")

library(ggpubr)
library("grid")
library(ggthemes) # Load
library(cowplot)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2,2)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(ggdraw(p_b1), vp = vplayout(1,1))
print(ggdraw(p_b2), vp = vplayout(1,2))
print(ggdraw(p_b3), vp = vplayout(2,1))
print(venn_F_p,  vp = vplayout(1,1))
print(venn_L_p,  vp = vplayout(1,2))
print(venn_R_p,  vp = vplayout(2,1))
print(venn_MST_p,vp = vplayout(2,2))

#print(venn_coGO_p,vp = vplayout(2,2))

#print(venn_coGR_p,vp = vplayout(3,2))


png("Capsella_PDF/Fig2 Volcano and Venn res600.png",
    width     = 16,
    height    = 10,
    units     = "in",
    res       = 600,
    pointsize = 4
)

grid.newpage()
pushViewport(viewport(layout = grid.layout(3,5)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(ggdraw(p_b1), vp = vplayout(1,1:4))
print(ggdraw(p_b2), vp = vplayout(2,1:4))
print(ggdraw(p_b3), vp = vplayout(3,1:4))
print(ggdraw(F_CRCG_p2), vp = vplayout(1,1))
print(ggdraw(F_CRCO_p2), vp = vplayout(1,2))
print(ggdraw(F_COCG_p2), vp = vplayout(1,3))
print(venn_F_p,  vp = vplayout(1,4))
print(venn_MST_p,vp = vplayout(1,5))
print(ggdraw(L_CRCG_p2), vp = vplayout(2,1))
print(ggdraw(L_CRCO_p2), vp = vplayout(2,2))
print(ggdraw(L_COCG_p2), vp = vplayout(2,3))
print(venn_L_p,  vp = vplayout(2,4))
print(ggdraw(R_CRCG_p2), vp = vplayout(3,1))
print(ggdraw(R_CRCO_p2), vp = vplayout(3,2))
print(ggdraw(R_COCG_p2), vp = vplayout(3,3))
print(venn_R_p,  vp = vplayout(3,4))

dev.off()

# Separate the unique and pleiotropy genes
for (co in c("coCRCO", "coCGCO", "coCGCR")) { # co start
  if(co == "coCRCO"){ # ifelse 1
    co_F <- F.mst
    co_L <- L.mst
    co_R <- R.mst
  } else if (co == "coCGCO") {
    co_F <- F.coGO
    co_L <- L.coGO
    co_R <- R.coGO
  } else {
    co_F <- F.coGR
    co_L <- L.coGR
    co_R <- R.coGR
  } # ifelse 1
  ## in flower
  F.u.1 <- row.names(co_F) %!in% row.names(co_L) # flower not in leaf
  F.u.2 <- row.names(co_F) %!in% row.names(co_R) # flower not in root
  F.u.keep <- cbind(F.u.1 & F.u.2) # flower not in either
  F.unique <- co_F[F.u.keep,] # coCRCO (MST) genes with unique function
  F.pleiotropy <- co_F[!F.u.keep,] # coCRCO (MST) genes with pleiotropy function
  
  assign(paste0("F_", co, ".unique"), F.unique)
  assign(paste0("F_", co, ".pleiotropy"), F.pleiotropy)
  ## in leaf
  L.u.1 <- row.names(co_L) %!in% row.names(co_F) # leaf not in flower
  L.u.2 <- row.names(co_L) %!in% row.names(co_R) # leaf not in root
  L.u.keep <- cbind(L.u.1 & L.u.2) # leaf not in either
  L.unique <- co_L[L.u.keep,] # coCRCO (MST) genes with unique function
  L.pleiotropy <- co_L[!L.u.keep,] # coCRCO (MST) genes with pleiotropy function
  
  assign(paste0("L_", co, ".unique"), L.unique)
  assign(paste0("L_", co, ".pleiotropy"), L.pleiotropy)
  ## in root
  R.u.1 <- row.names(co_R) %!in% row.names(co_F) # root not in flower
  R.u.2 <- row.names(co_R) %!in% row.names(co_L) # root not in leaf
  R.u.keep <- cbind(R.u.1 & R.u.2) # root not in either
  R.unique <- co_R[R.u.keep,] # coCRCO (MST) genes with unique function
  R.pleiotropy <- co_R[!R.u.keep,] # coCRCO (MST) genes with pleiotropy function
  
  assign(paste0("R_", co, ".unique"), R.unique)
  assign(paste0("R_", co, ".pleiotropy"), R.pleiotropy)
  
  # heatmap
  library(ComplexHeatmap)
  
  for (tissue in c("F","L","R")) { # for tissue start
    meanTMM <- read.table(paste0("OutputData/MeanTMM_", tissue, ".txt"))
    keep_u <- row.names(meanTMM) %in% row.names(F.unique)
    keep_p <- row.names(meanTMM) %in% row.names(F.pleiotropy)
    
    # unique
    TMM_U <- meanTMM[keep_u,]
    scaleTMM_U <- scaleRow(TMM_U)
    TMM_U_row_dend = as.dendrogram(hclust(dist(scaleTMM_U))) # generate hierarchical dendrogram for rows
    # Mean_F_row_dend = color_branches(Mean_F_row_dend, k = 4 ) # `color_branches()` returns a dendrogram object
    TMM_U_column_dend = as.dendrogram(hclust(dist(t(TMM_U)))) # generate hierarchical dendrogram for rows - t transpose data
    TMM_U_column_dend = color_branches(TMM_U_column_dend, k = 2) # `color_branches()` returns a dendrogram object
    if(tissue == "F"){ # ifelse 2
      TMM_U_annotation = HeatmapAnnotation(
        Tissues = rep("Flower", 3),
        Species = c("CG","CR","CO"),
        col = list(Tissues=c("Flower" = "#fa990e"),
                   Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
        ),
        annotation_legend_param = list(
          Tissues = list(direction = "horizontal"),
          Species = list(nrow = 1)),
        show_annotation_name = FALSE         # NOT show annotation label
      ) 
    }else if (tissue == "L"){
      TMM_U_annotation = HeatmapAnnotation(
        Tissues = rep("Leaf", 3),
        Species = c("CG","CR","CO"),
        col = list(Tissues=c("Leaf" = "#155800"),
                   Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
        ),
        annotation_legend_param = list(
          Tissues = list(direction = "horizontal"),
          Species = list(nrow = 1)),
        show_annotation_name = FALSE         # NOT show annotation label
      )
    }else{
      TMM_U_annotation = HeatmapAnnotation(
        Tissues = rep("Root", 3),
        Species = c("CG","CR","CO"),
        col = list(Tissues=c("Root" = "#089400"),
                   Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
        ),
        annotation_legend_param = list(
          Tissues = list(direction = "horizontal"),
          Species = list(nrow = 1)),
        show_annotation_name = FALSE         # NOT show annotation label
      )
    } # if else 2
      # Pleiotropy
      TMM_P <- meanTMM[keep_p,]
      scaleTMM_P <- scaleRow(TMM_P)
      TMM_P_row_dend = as.dendrogram(hclust(dist(scaleTMM_P))) # generate hierarchical dendrogram for rows
      # Mean_F_row_dend = color_branches(Mean_F_row_dend, k = 4 ) # `color_branches()` returns a dendrogram object
      TMM_P_column_dend = as.dendrogram(hclust(dist(t(TMM_P)))) # generate hierarchical dendrogram for rows - t transpose data
      TMM_P_column_dend = color_branches(TMM_P_column_dend, k = 2) # `color_branches()` returns a dendrogram object
      if(tissue == "F"){ # ifelse 3
        TMM_P_annotation = HeatmapAnnotation(
          Tissues = rep("Flower", 3),
          Species = c("CG","CR","CO"),
          col = list(Tissues=c("Flower" = "#fa990e"),
                     Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
          ),
          annotation_legend_param = list(
            Tissues = list(direction = "horizontal"),
            Species = list(nrow = 1)),
          show_annotation_name = FALSE         # NOT show annotation label
        ) 
      }else if (tissue == "L"){
        TMM_P_annotation = HeatmapAnnotation(
          Tissues = rep("Leaf", 3),
          Species = c("CG","CR","CO"),
          col = list(Tissues=c("Leaf" = "#155800"),
                     Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
          ),
          annotation_legend_param = list(
            Tissues = list(direction = "horizontal"),
            Species = list(nrow = 1)),
          show_annotation_name = FALSE         # NOT show annotation label
        )
      }else{
        TMM_P_annotation = HeatmapAnnotation(
          Tissues = rep("Root", 3),
          Species = c("CG","CR","CO"),
          col = list(Tissues=c("Root" = "#089400"),
                     Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
          ),
          annotation_legend_param = list(
            Tissues = list(direction = "horizontal"),
            Species = list(nrow = 1)),
          show_annotation_name = FALSE         # NOT show annotation label
        )
      } # ifelse 3
    
    P1 <- draw(Heatmap(scaleTMM_U,                                        # using scaled data
                       name = "TMM", 
                       cluster_rows = TMM_U_row_dend,                         # replace default dendrogram to custom one
                       cluster_columns = TMM_U_column_dend,
                       #column_order = sort(colnames(TMM_F)),             # reorder columns by column names
                       show_row_names = FALSE,                         # hide row names
                       show_column_names = FALSE,                         # hide column names
                       row_title = NULL,                                  # hide cluster order numbers
                       column_title = NULL,
                       row_dend_width = unit(0.5, "cm"),                    # length of dendrogram
                       column_dend_height = unit(1, "cm"),
                       row_split = 2,                                     # split row dendrogram to 3 by hierarchical relationship 
                       column_split = 2,
                       top_annotation = TMM_U_annotation,                      # add column annotation, ie, colour of species and tissues
                       heatmap_legend_param = list(direction = "horizontal"),
                       
    ),
    heatmap_legend_side = "bottom", 
    annotation_legend_side = "bottom",
    #show_annotation_legend = FALSE,                            # hide the annotation legend (tissues and species)
    #show_heatmap_legend = FALSE,                               # hide the heatmap legend (mapping value color red2blue)
    newpage = FALSE                                            # don't create a new page for plot
    )
    
    P2 <- draw(Heatmap(scaleTMM_P,                                        # using scaled data
                       name = "TMM", 
                       cluster_rows = TMM_P_row_dend,                         # replace default dendrogram to custom one
                       cluster_columns = TMM_P_column_dend,
                       #column_order = sort(colnames(TMM_F)),             # reorder columns by column names
                       show_row_names = FALSE,                         # hide row names
                       show_column_names = FALSE,                         # hide column names
                       row_title = NULL,                                  # hide cluster order numbers
                       column_title = NULL,
                       row_dend_width = unit(0.5, "cm"),                    # length of dendrogram
                       column_dend_height = unit(1, "cm"),
                       row_split = 2,                                     # split row dendrogram to 3 by hierarchical relationship 
                       column_split = 2,
                       top_annotation = TMM_P_annotation,                      # add column annotation, ie, colour of species and tissues
                       heatmap_legend_param = list(direction = "horizontal"),
                       
    ),
    heatmap_legend_side = "bottom", 
    annotation_legend_side = "bottom",
    #show_annotation_legend = FALSE,                            # hide the annotation legend (tissues and species)
    #show_heatmap_legend = FALSE,                               # hide the heatmap legend (mapping value color red2blue)
    newpage = FALSE                                            # don't create a new page for plot
    )
    
    assign(paste0("p_", co,"_",tissue,"_U"), P1)
    assign(paste0("p_", co,"_",tissue,"_P"), P2)
    
  } # for tissue end
} # for co end



# drawing background colour
p_b4 <- ggplot(F.mst) +
  geom_rect(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,   fill = "#f0eeef") # light gray

png("Capsella_PDF/UniqueAndPleiotropy.png",
    width     = 16,
    height    = 16,
    units     = "in",
    res       = 600,
    pointsize = 4
)

grid.newpage()
pushViewport(viewport(layout = grid.layout(4,6)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(ggdraw(p_b4), vp = vplayout(2:4,1))
print(ggdraw(p_b4), vp = vplayout(2:4,3))
print(ggdraw(p_b4), vp = vplayout(2:4,5))
print(venn_MST_p,vp = vplayout(1,1:2))
print(venn_coGO_p,vp = vplayout(1,3:4))
print(venn_coGR_p,vp = vplayout(1,5:6))
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(p_coCRCO_F_U, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(p_coCRCO_F_P, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
draw(p_coCRCO_L_U, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 2))
draw(p_coCRCO_L_P, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 1))
draw(p_coCRCO_R_U, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 2))
draw(p_coCRCO_R_P, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
draw(p_coCGCO_F_U, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 4))
draw(p_coCGCO_F_P, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 3))
draw(p_coCGCO_L_U, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 4))
draw(p_coCGCO_L_P, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 3))
draw(p_coCGCO_R_U, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 4))
draw(p_coCGCO_R_P, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 5))
draw(p_coCGCR_F_U, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 6))
draw(p_coCGCR_F_P, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 5))
draw(p_coCGCR_L_U, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 6))
draw(p_coCGCR_L_P, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 5))
draw(p_coCGCR_R_U, newpage = F)
upViewport()
pushViewport(viewport(layout.pos.row = 4, layout.pos.col = 6))
draw(p_coCGCR_R_P, newpage = F)
upViewport()

dev.off()

F_coCRCO_keep_u <- row.names(MeanTMM_F) %in% row.names(F_coCRCO.unique)
F_coCRCO_keep_p <- row.names(MeanTMM_F) %in% row.names(F_coCRCO.pleiotropy)

TMM_F_coCRCO_U <- MeanTMM_F[F_coCRCO_keep_u,]
head(TMM_F_coCRCO_U)
dim(TMM_F_coCRCO_U)
scaleTMM_F_coCRCO_U <- scaleRow(TMM_F_coCRCO_U)
head(scaleTMM_F_coCRCO_U)
TMM_F_coCRCO_U_row_dend = as.dendrogram(hclust(dist(scaleTMM_F_coCRCO_U))) # generate hierarchical dendrogram for rows
# Mean_F_row_dend = color_branches(Mean_F_row_dend, k = 4 ) # `color_branches()` returns a dendrogram object
TMM_F_coCRCO_U_column_dend = as.dendrogram(hclust(dist(t(TMM_F_coCRCO_U)))) # generate hierarchical dendrogram for rows - t transpose data
TMM_F_coCRCO_U_column_dend = color_branches(TMM_F_coCRCO_U_column_dend, k = 2) # `color_branches()` returns a dendrogram object
TMM_F_coCRCO_U_annotation = HeatmapAnnotation(
  Tissues = rep("Flower", 3),
  Species = c("CG","CR","CO"),
  col = list(Tissues=c("Flower" = "#fa990e"),
             Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
  ),
  annotation_legend_param = list(
    Tissues = list(direction = "horizontal"),
    Species = list(nrow = 1)),
  show_annotation_name = FALSE         # NOT show annotation label
)

P1 <- draw(Heatmap(scaleTMM_F_coCRCO_U,                                        # using scaled data
             name = "TMM", 
             cluster_rows = TMM_F_coCRCO_U_row_dend,                         # replace default dendrogram to custom one
             cluster_columns = TMM_F_coCRCO_U_column_dend,
             #column_order = sort(colnames(TMM_F)),             # reorder columns by column names
             show_row_names = FALSE,                         # hide row names
             show_column_names = FALSE,                         # hide column names
             row_title = NULL,                                  # hide cluster order numbers
             column_title = NULL,
             row_dend_width = unit(0.5, "cm"),                    # length of dendrogram
             column_dend_height = unit(1, "cm"),
             row_split = 2,                                     # split row dendrogram to 3 by hierarchical relationship 
             column_split = 2,
             top_annotation = TMM_F_coCRCO_U_annotation,                      # add column annotation, ie, colour of species and tissues
             heatmap_legend_param = list(direction = "horizontal"),
             
),
heatmap_legend_side = "bottom", 
annotation_legend_side = "bottom",
#show_annotation_legend = FALSE,                            # hide the annotation legend (tissues and species)
#show_heatmap_legend = FALSE,                               # hide the heatmap legend (mapping value color red2blue)
newpage = FALSE                                            # don't create a new page for plot
)

# Heatmap for Selfing syndrome genes. 
# run Heatmap for all genes.R first
library(ComplexHeatmap)
MST <- read.csv("OutputData/MatingSystemTransition_gene727.csv", row.names = 1)
head(MST)
head(MeanTMM)
MeanTMM_F <- MeanTMM[,1:3]
head(MeanTMM_F)

keep <- row.names(MeanTMM_F) %in% row.names(MST)

MST_F <- MeanTMM_F[keep,]
head(MST_F)
dim(MST_F)
scaleMST_F <- scaleRow(MST_F)
head(scaleMST_F)
Mean_MST_F_row_dend = as.dendrogram(hclust(dist(scaleMST_F))) # generate hierarchical dendrogram for rows
# Mean_F_row_dend = color_branches(Mean_F_row_dend, k = 4 ) # `color_branches()` returns a dendrogram object
Mean_MST_F_column_dend = as.dendrogram(hclust(dist(t(MST_F)))) # generate hierarchical dendrogram for rows - t transpose data
Mean_MST_F_column_dend = color_branches(Mean_MST_F_column_dend, k = 2) # `color_branches()` returns a dendrogram object
Mean_MST_F_annotation = HeatmapAnnotation(
  Tissues = rep("Flower", 3),
  Species = c("CG","CR","CO"),
  col = list(Tissues=c("Flower" = "#fa990e"),
             Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
  ),
  annotation_legend_param = list(
    Tissues = list(direction = "horizontal"),
    Species = list(nrow = 1)),
  show_annotation_name = FALSE         # NOT show annotation label
)

draw(Heatmap(scaleMST_F,                                        # using scaled data
             name = "TMM", 
             cluster_rows = Mean_MST_F_row_dend,                         # replace default dendrogram to custom one
             cluster_columns = Mean_MST_F_column_dend,
             #column_order = sort(colnames(TMM_F)),             # reorder columns by column names
             show_row_names = FALSE,                         # hide row names
             show_column_names = FALSE,                         # hide column names
             row_title = NULL,                                  # hide cluster order numbers
             column_title = NULL,
             row_dend_width = unit(0.5, "cm"),                    # length of dendrogram
             column_dend_height = unit(1, "cm"),
             row_split = 2,                                     # split row dendrogram to 3 by hierarchical relationship 
             column_split = 2,
             top_annotation = Mean_MST_F_annotation,                      # add column annotation, ie, colour of species and tissues
             heatmap_legend_param = list(direction = "horizontal"),
             
),
heatmap_legend_side = "bottom", 
annotation_legend_side = "bottom",
#show_annotation_legend = FALSE,                            # hide the annotation legend (tissues and species)
#show_heatmap_legend = FALSE,                               # hide the heatmap legend (mapping value color red2blue)
newpage = FALSE                                            # don't create a new page for plot
)


library(ggpubr)
library("grid")
library(ggthemes) # Load

png("Capsella_PDF/Fig2 res600.png",
    width     = 16,
    height    = 10,
    units     = "in",
    res       = 600,
    pointsize = 4
)

grid.newpage()
pushViewport(viewport(layout = grid.layout(3,5)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(ggdraw(p_b1), vp = vplayout(1,1:4))
print(ggdraw(p_b2), vp = vplayout(2,1:4))
print(ggdraw(p_b3), vp = vplayout(3,1:4))
print(ggdraw(F_CGCR_p2), vp = vplayout(1,1))
print(ggdraw(F_CRCO_p2), vp = vplayout(1,2))
print(ggdraw(F_CGCO_p2), vp = vplayout(1,3))
print(venn_F_p,  vp = vplayout(1,4))
print(venn_MST_p,vp = vplayout(1,5))
print(ggdraw(L_CGCR_p2), vp = vplayout(2,1))
print(ggdraw(L_CRCO_p2), vp = vplayout(2,2))
print(ggdraw(L_CGCO_p2), vp = vplayout(2,3))
print(venn_L_p,  vp = vplayout(2,4))
print(ggdraw(R_CGCR_p2), vp = vplayout(3,1))
print(ggdraw(R_CRCO_p2), vp = vplayout(3,2))
print(ggdraw(R_CGCO_p2), vp = vplayout(3,3))
print(venn_R_p,  vp = vplayout(3,4))
pushViewport(viewport(layout.pos.row = 2:3, layout.pos.col = 5))
draw(p_coCGCO_F_P, newpage = F)
draw(Heatmap(scaleMST_F,                                        # using scaled data
             name = "TMM", 
             cluster_rows = Mean_MST_F_row_dend,                         # replace default dendrogram to custom one
             cluster_columns = Mean_MST_F_column_dend,
             #column_order = sort(colnames(TMM_F)),             # reorder columns by column names
             show_row_names = FALSE,                         # hide row names
             show_column_names = FALSE,                         # hide column names
             row_title = NULL,                                  # hide cluster order numbers
             column_title = NULL,
             row_dend_width = unit(0.5, "cm"),                    # length of dendrogram
             column_dend_height = unit(1, "cm"),
             row_split = 2,                                     # split row dendrogram to 3 by hierarchical relationship 
             column_split = 2,
             top_annotation = Mean_MST_F_annotation,                      # add column annotation, ie, colour of species and tissues
             heatmap_legend_param = list(direction = "horizontal"),
             
),
heatmap_legend_side = "bottom", 
annotation_legend_side = "bottom",
#show_annotation_legend = FALSE,                            # hide the annotation legend (tissues and species)
#show_heatmap_legend = FALSE,                               # hide the heatmap legend (mapping value color red2blue)
newpage = FALSE                                            # don't create a new page for plot
)

upViewport()
upViewport()

dev.off()

grid.newpage()
pushViewport(viewport(layout = grid.layout(3,1)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(venn_MST_p,vp = vplayout(1,1))
pushViewport(viewport(layout.pos.row = 2:3, layout.pos.col = 1))
draw(Heatmap(scaleMST_F,                                        # using scaled data
             name = "TMM", 
             cluster_rows = Mean_MST_F_row_dend,                         # replace default dendrogram to custom one
             cluster_columns = Mean_MST_F_column_dend,
             #column_order = sort(colnames(TMM_F)),             # reorder columns by column names
             show_row_names = FALSE,                         # hide row names
             show_column_names = FALSE,                         # hide column names
             row_title = NULL,                                  # hide cluster order numbers
             column_title = NULL,
             row_dend_width = unit(0.5, "cm"),                    # length of dendrogram
             column_dend_height = unit(1, "cm"),
             row_split = 2,                                     # split row dendrogram to 3 by hierarchical relationship 
             column_split = 2,
             top_annotation = Mean_MST_F_annotation,                      # add column annotation, ie, colour of species and tissues
             heatmap_legend_param = list(direction = "horizontal"),
             
),
heatmap_legend_side = "bottom", 
annotation_legend_side = "bottom",
#show_annotation_legend = FALSE,                            # hide the annotation legend (tissues and species)
#show_heatmap_legend = FALSE,                               # hide the heatmap legend (mapping value color red2blue)
newpage = FALSE                                            # don't create a new page for plot
)

upViewport()
upViewport()

dev.off()









# Compare the correlation between CRCG and COCG 
##################################################
# correlation in gene expression changes between CRCG and COCG in Flower
DE_F_CRCG <- read.table("OutputData/DE_F_CRCG_FC2.txt", header = T, sep = "\t", row.names = 1)
DE_F_COCG <- read.table("OutputData/DE_F_COCG_FC2.txt", header = T, sep = "\t", row.names = 1)
head(DE_F_CRCG)
head(DE_F_COCG)

F.CRCG <- DE_F_CRCG
F.CRCG$gene <- rownames(F.CRCG)
F.CRCG <- F.CRCG[,c("gene","logFC")]
rownames(F.CRCG) <- NULL
names(F.CRCG) <- c("gene","CRCGFC")
head(F.CRCG)

F.COCG <- DE_F_COCG
F.COCG$gene <- rownames(F.COCG)
F.COCG <- F.COCG[,c("gene","logFC")]
rownames(F.COCG) <- NULL
names(F.COCG) <- c("gene","COCGFC")
head(F.COCG)

F.CRCO2CG <- merge(F.CRCG, F.COCG, by="gene")
# remove values of these two comparison in opposite direction
F.CRCO2CG <- F.CRCO2CG[!(F.CRCO2CG$CRCGFC > 0 & F.CRCO2CG$COCGFC < 0),] 
F.CRCO2CG <- F.CRCO2CG[!(F.CRCO2CG$COCGFC > 0 & F.CRCO2CG$CRCGFC < 0),] 
head(F.CRCO2CG)
dim(F.CRCO2CG)

library("ggpubr")
PF <- ggscatter(F.CRCO2CG, x = "CRCGFC", y = "COCGFC", 
          color = "firebrick1",
          alpha = 0.5,
          size = 0.7,
          add = "reg.line", 
          add.params = list(color = "black",         # change line colour
                            linetype = "dashed"),   # change line type
          #conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.method = "pearson",
          #xlim = c(-10,10), ylim =c(-10,10),
          title = "Flowers",
          xlab = "log2(CR/CG)", ylab = "log2(CO/CG)")

cor(F.CRCO2CG$COCGFC,F.CRCO2CG$CRCGFC)
relation <- lm(F.CRCO2CG$COCGFC~F.CRCO2CG$CRCGFC)
print(summary(relation))
# correlation in gene expression changes between CRCG and COCG in Leaf
DE_L_CRCG <- read.table("OutputData/DE_L_CRCG_FC2FDR.05.txt", header = T, sep = "\t", row.names = 1)
DE_L_COCG <- read.table("OutputData/DE_L_COCG_FC2FDR.05.txt", header = T, sep = "\t", row.names = 1)
head(DE_L_CRCG)
head(DE_L_COCG)

L.CRCG <- DE_L_CRCG
L.CRCG$gene <- rownames(L.CRCG)
L.CRCG <- L.CRCG[,c("gene","logFC")]
rownames(L.CRCG) <- NULL
names(L.CRCG) <- c("gene","CRCGFC")
head(L.CRCG)

L.COCG <- DE_L_COCG
L.COCG$gene <- rownames(L.COCG)
L.COCG <- L.COCG[,c("gene","logFC")]
rownames(L.COCG) <- NULL
names(L.COCG) <- c("gene","COCGFC")
head(L.COCG)

L.CRCO2CG <- merge(L.CRCG, L.COCG, by="gene")
# remove values of these two comparison in opposite direction
L.CRCO2CG <- L.CRCO2CG[!(L.CRCO2CG$CRCGFC > 0 & L.CRCO2CG$COCGFC < 0),] 
L.CRCO2CG <- L.CRCO2CG[!(L.CRCO2CG$COCGFC > 0 & L.CRCO2CG$CRCGFC < 0),] 
head(L.CRCO2CG)
dim(L.CRCO2CG)

library("ggpubr")
PL <- ggscatter(L.CRCO2CG, x = "CRCGFC", y = "COCGFC", 
          color = "forestgreen",
          alpha = 0.5,
          size = 0.7,
          add = "reg.line", 
          add.params = list(color = "black",         # change line colour
                            linetype = "dashed"),   # change line type
          #conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.method = "pearson",
          #xlim = c(-10,10), ylim =c(-10,10),
          title = "Leaves",
          xlab = "log2(CR/CG)", ylab = "log2(CO/CG)")
# correlation in gene expression changes between CRCG and COCG in Root
DE_R_CRCG <- read.table("OutputData/DE_R_CRCG_FC2FDR.05.txt", header = T, sep = "\t", row.names = 1)
DE_R_COCG <- read.table("OutputData/DE_R_COCG_FC2FDR.05.txt", header = T, sep = "\t", row.names = 1)
head(DE_R_CRCG)
head(DE_R_COCG)

R.CRCG <- DE_R_CRCG
R.CRCG$gene <- rownames(R.CRCG)
R.CRCG <- R.CRCG[,c("gene","logFC")]
rownames(R.CRCG) <- NULL
names(R.CRCG) <- c("gene","CRCGFC")
head(R.CRCG)

R.COCG <- DE_R_COCG
R.COCG$gene <- rownames(R.COCG)
R.COCG <- R.COCG[,c("gene","logFC")]
rownames(R.COCG) <- NULL
names(R.COCG) <- c("gene","COCGFC")
head(R.COCG)

R.CRCO2CG <- merge(R.CRCG, R.COCG, by="gene")
# remove values of these two comparison in opposite direction
R.CRCO2CG <- R.CRCO2CG[!(R.CRCO2CG$CRCGFC > 0 & R.CRCO2CG$COCGFC < 0),] 
R.CRCO2CG <- R.CRCO2CG[!(R.CRCO2CG$COCGFC > 0 & R.CRCO2CG$CRCGFC < 0),] 
head(R.CRCO2CG)
dim(R.CRCO2CG)

library("ggpubr")
PR <- ggscatter(R.CRCO2CG, x = "CRCGFC", y = "COCGFC", 
          color = "#da7702",
          alpha = 0.5,
          size = 0.7,
          add = "reg.line", 
          add.params = list(color = "black",         # change line colour
                            linetype = "dashed"),   # change line type
          conf.int = TRUE, 
          cor.coef = TRUE, 
          cor.method = "pearson",
          #xlim = c(-10,10), ylim =c(-10,10),
          title = "Roots",
          xlab = "log2(CR/CG)", ylab = "log2(CO/CG)")

ggarrange(PF, PL, PR, ncol = 3, nrow = 1)

relation <- lm(R.CRCO2CG$CRCGFC~R.CRCO2CG$COCGFC)

print(summary(relation))
# Consider all genes

head(DE)
# correlation in gene expression changes between CRCG and COCG in Flower
DE_F_CRCG <- read.table("OutputData/DE_F_CRCG.txt", header = T, sep = "\t", row.names = 1)
DE_F_COCG <- read.table("OutputData/DE_F_COCG.txt", header = T, sep = "\t", row.names = 1)

F.CRCG <- DE_F_CRCG
F.CRCG$gene <- rownames(F.CRCG)
F.CRCG <- F.CRCG[,c("gene","logFC")]
rownames(F.CRCG) <- NULL
names(F.CRCG) <- c("gene","CRCGFC")
head(F.CRCG)

F.COCG <- DE_F_COCG
F.COCG$gene <- rownames(F.COCG)
F.COCG <- F.COCG[,c("gene","logFC")]
rownames(F.COCG) <- NULL
names(F.COCG) <- c("gene","COCGFC")
head(F.COCG)

F.CRCO2CG <- merge(F.CRCG, F.COCG, by="gene")

F.CRCO2CG$Label <- "All"
F.CRCO2CG[(F.CRCO2CG$CRCGFC > 1 & F.CRCO2CG$COCGFC > 1),]$Label <- "FC > 2"
F.CRCO2CG[(F.CRCO2CG$COCGFC < -1 & F.CRCO2CG$CRCGFC < -1),]$Label <- "FC > 2"

head(F.CRCO2CG)

with(F.CRCO2CG, plot(CRCGFC, COCGFC, 
                     pch=20,  # dot type
                     cex=0.5, # dot size
                     col = alpha("gray", 0.5),
                     main="Flowers",
                     xlab=expression(log[2]~"("~italic(CR/CG)~")"), 
                     ylab=expression(log[2]~"("~italic(CO/CG)~")")
                     ))
with(subset(F.CRCO2CG,  F.CRCO2CG$Label == "FC > 2"), 
     points(CRCGFC, COCGFC, 
            pch=20, 
            cex=0.5, # dot size
            col=alpha("tomato", 0.5)))
with(subset(F.CRCO2CG),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "gray34",            # Modify color
            lty = "dashed",         # Modify line type
            lwd = 2))
with(subset(F.CRCO2CG, F.CRCO2CG$Label == "FC > 2"),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "red3",            # Modify color
            lty = "dashed",         # Modify line type
            lwd = 2))
text(-8, 10, expression(italic(R)~"= 0.42,"~beta~"= 0.59"), adj = 0)
text(-8, 8, expression(italic(R)~"= 0.87,"~beta~"= 1.10"), adj = 0, col = "red3")


cor(F.CRCO2CG$COCGFC,F.CRCO2CG$CRCGFC)
relation <- lm(F.CRCO2CG$COCGFC~F.CRCO2CG$CRCGFC)
print(summary(relation))
relation <- lm(F.CRCO2CG[F.CRCO2CG$Label == "FC > 2",]$COCGFC~F.CRCO2CG[F.CRCO2CG$Label == "FC > 2",]$CRCGFC)
print(summary(relation))
# correlation in gene expression changes between CRCG and COCG in Leaf
DE_L_CRCG <- read.table("OutputData/DE_L_CRCG.txt", header = T, sep = "\t", row.names = 1)
DE_L_COCG <- read.table("OutputData/DE_L_COCG.txt", header = T, sep = "\t", row.names = 1)
head(DE_L_CRCG)
head(DE_L_COCG)

L.CRCG <- DE_L_CRCG
L.CRCG$gene <- rownames(L.CRCG)
L.CRCG <- L.CRCG[,c("gene","logFC")]
rownames(L.CRCG) <- NULL
names(L.CRCG) <- c("gene","CRCGFC")
head(L.CRCG)

L.COCG <- DE_L_COCG
L.COCG$gene <- rownames(L.COCG)
L.COCG <- L.COCG[,c("gene","logFC")]
rownames(L.COCG) <- NULL
names(L.COCG) <- c("gene","COCGFC")
head(L.COCG)

L.CRCO2CG <- merge(L.CRCG, L.COCG, by="gene")

L.CRCO2CG$Label <- "All"
L.CRCO2CG[(L.CRCO2CG$CRCGFC > 1 & L.CRCO2CG$COCGFC > 1),]$Label <- "FC > 2"
L.CRCO2CG[(L.CRCO2CG$COCGFC < -1 & L.CRCO2CG$CRCGFC < -1),]$Label <- "FC > 2"

head(L.CRCO2CG)

with(L.CRCO2CG, plot(CRCGFC, COCGFC, 
                     pch=20,  # dot type
                     cex=0.5, # dot size
                     col = alpha("gray", 0.5),
                     main="Leaves",
                     xlab=expression(log[2]~"("~italic(CR/CG)~")"), 
                     ylab=expression(log[2]~"("~italic(CO/CG)~")")
))
with(subset(L.CRCO2CG,  L.CRCO2CG$Label == "FC > 2"), 
     points(CRCGFC, COCGFC, 
            pch=20, 
            cex=0.5, # dot size
            col=alpha("#1eb684", 0.5)))
with(subset(L.CRCO2CG),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "gray34",            # Modify color
            lty = "dashed",         # Modify line type light green
            lwd = 2))
with(subset(L.CRCO2CG, L.CRCO2CG$Label == "FC > 2"),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "#046c46",            # Modify color dark green
            lty = "dashed",         # Modify line type
            lwd = 2))
text(-10, 10, expression(italic(R)~"= 0.26,"~beta~"= 0.35"), adj = 0)
text(-10, 8, expression(italic(R)~"= 0.82,"~beta~"= 0.98"), adj = 0, col = "#046c46")


cor(L.CRCO2CG$COCGFC,L.CRCO2CG$CRCGFC)
relation <- lm(L.CRCO2CG$COCGFC~L.CRCO2CG$CRCGFC)
print(summary(relation))
cor(L.CRCO2CG[L.CRCO2CG$Label == "FC > 2",]$COCGFC,L.CRCO2CG[L.CRCO2CG$Label == "FC > 2",]$CRCGFC)
relation <- lm(L.CRCO2CG[L.CRCO2CG$Label == "FC > 2",]$COCGFC~L.CRCO2CG[L.CRCO2CG$Label == "FC > 2",]$CRCGFC)
print(summary(relation))
# correlation in gene expression changes between CRCG and COCG in Root
DE_R_CRCG <- read.table("OutputData/DE_R_CRCG.txt", header = T, sep = "\t", row.names = 1)
DE_R_COCG <- read.table("OutputData/DE_R_COCG.txt", header = T, sep = "\t", row.names = 1)
head(DE_R_CRCG)
head(DE_R_COCG)

R.CRCG <- DE_R_CRCG
R.CRCG$gene <- rownames(R.CRCG)
R.CRCG <- R.CRCG[,c("gene","logFC")]
rownames(R.CRCG) <- NULL
names(R.CRCG) <- c("gene","CRCGFC")
head(R.CRCG)

R.COCG <- DE_R_COCG
R.COCG$gene <- rownames(R.COCG)
R.COCG <- R.COCG[,c("gene","logFC")]
rownames(R.COCG) <- NULL
names(R.COCG) <- c("gene","COCGFC")
head(R.COCG)

R.CRCO2CG <- merge(R.CRCG, R.COCG, by="gene")

R.CRCO2CG$Label <- "All"
R.CRCO2CG[(R.CRCO2CG$CRCGFC > 1 & R.CRCO2CG$COCGFC > 1),]$Label <- "FC > 2"
R.CRCO2CG[(R.CRCO2CG$COCGFC < -1 & R.CRCO2CG$CRCGFC < -1),]$Label <- "FC > 2"

head(R.CRCO2CG)

with(R.CRCO2CG, plot(CRCGFC, COCGFC, 
                     pch=20,  # dot type
                     cex=0.5, # dot size
                     col = alpha("gray", 0.5),
                     main="Leaves",
                     xlab=expression(log[2]~"("~italic(CR/CG)~")"), 
                     ylab=expression(log[2]~"("~italic(CO/CG)~")")
))
with(subset(R.CRCO2CG,  R.CRCO2CG$Label == "FC > 2"), 
     points(CRCGFC, COCGFC, 
            pch=20, 
            cex=0.5, # dot size
            col=alpha("#af8d2c", 0.5))) # Light brown
with(subset(R.CRCO2CG),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "gray34",            # Modify color
            lty = "dashed",         # Modify line type 
            lwd = 2))
with(subset(R.CRCO2CG, R.CRCO2CG$Label == "FC > 2"),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "#856201",            # Modify color dark brown
            lty = "dashed",         # Modify line type
            lwd = 2))
text(-10, 10, expression(italic(R)~"= 0.23,"~beta~"= 0.30"), adj = 0)
text(-10, 8, expression(italic(R)~"= 0.83,"~beta~"= 0.95"), adj = 0, col = "#856201")


cor(R.CRCO2CG$COCGFC,R.CRCO2CG$CRCGFC)
relation <- lm(R.CRCO2CG$COCGFC~R.CRCO2CG$CRCGFC)
print(summary(relation))
cor(R.CRCO2CG[R.CRCO2CG$Label == "FC > 2",]$COCGFC,R.CRCO2CG[R.CRCO2CG$Label == "FC > 2",]$CRCGFC)
relation <- lm(R.CRCO2CG[R.CRCO2CG$Label == "FC > 2",]$COCGFC~R.CRCO2CG[R.CRCO2CG$Label == "FC > 2",]$CRCGFC)
print(summary(relation))

#########################################################################
## Plot all figures together
#opar <- par()
par(mfrow = c(1, 3))
# Flowers
with(F.CRCO2CG, plot(CRCGFC, COCGFC, 
                     pch=20,  # dot type
                     cex=0.5, # dot size
                     col = alpha("gray", 0.5),
                     main="Flowers",
                     xlab=expression(log[2]~"("~italic(CR/CG)~")"), 
                     ylab=expression(log[2]~"("~italic(CO/CG)~")")
))
with(subset(F.CRCO2CG,  F.CRCO2CG$Label == "FC > 2"), 
     points(CRCGFC, COCGFC, 
            pch=20, 
            cex=0.5, # dot size
            col=alpha("tomato", 0.5)))
with(subset(F.CRCO2CG),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "gray34",            # Modify color
            lty = "dashed",         # Modify line type
            lwd = 2))
with(subset(F.CRCO2CG, F.CRCO2CG$Label == "FC > 2"),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "red3",            # Modify color
            lty = "dashed",         # Modify line type
            lwd = 2))
text(-8, 10, expression(italic(R)~"= 0.42,"~beta~"= 0.59"), adj = 0)
text(-8, 8, expression(italic(R)~"= 0.87,"~beta~"= 1.10"), adj = 0, col = "red3")
# Leaves
with(L.CRCO2CG, plot(CRCGFC, COCGFC, 
                     pch=20,  # dot type
                     cex=0.5, # dot size
                     col = alpha("gray", 0.5),
                     main="Leaves",
                     xlab=expression(log[2]~"("~italic(CR/CG)~")"), 
                     ylab=expression(log[2]~"("~italic(CO/CG)~")")
))
with(subset(L.CRCO2CG,  L.CRCO2CG$Label == "FC > 2"), 
     points(CRCGFC, COCGFC, 
            pch=20, 
            cex=0.5, # dot size
            col=alpha("#1eb684", 0.5)))
with(subset(L.CRCO2CG),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "gray34",            # Modify color
            lty = "dashed",         # Modify line type light green
            lwd = 2))
with(subset(L.CRCO2CG, L.CRCO2CG$Label == "FC > 2"),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "#046c46",            # Modify color dark green
            lty = "dashed",         # Modify line type
            lwd = 2))
text(-10, 10, expression(italic(R)~"= 0.26,"~beta~"= 0.35"), adj = 0)
text(-10, 8, expression(italic(R)~"= 0.82,"~beta~"= 0.98"), adj = 0, col = "#046c46")
# Roots
with(R.CRCO2CG, plot(CRCGFC, COCGFC, 
                     pch=20,  # dot type
                     cex=0.5, # dot size
                     col = alpha("gray", 0.5),
                     main="Roots",
                     xlab=expression(log[2]~"("~italic(CR/CG)~")"), 
                     ylab=expression(log[2]~"("~italic(CO/CG)~")")
))
with(subset(R.CRCO2CG,  R.CRCO2CG$Label == "FC > 2"), 
     points(CRCGFC, COCGFC, 
            pch=20, 
            cex=0.5, # dot size
            col=alpha("#af8d2c", 0.5))) # Light brown
with(subset(R.CRCO2CG),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "gray34",            # Modify color
            lty = "dashed",         # Modify line type 
            lwd = 2))
with(subset(R.CRCO2CG, R.CRCO2CG$Label == "FC > 2"),
     abline(lm(COCGFC~CRCGFC),
            cex = 1, pch = 16,
            col = "#856201",            # Modify color dark brown
            lty = "dashed",         # Modify line type
            lwd = 2))
text(-10, 10, expression(italic(R)~"= 0.23,"~beta~"= 0.30"), adj = 0)
text(-10, 8, expression(italic(R)~"= 0.83,"~beta~"= 0.95"), adj = 0, col = "#856201")








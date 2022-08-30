###############################################################
#          Similarity and convergence indices                 #
###############################################################

# To quantify the similarity between each subgenome expression level and the expression level in the parental 
# species, we developed a similarity index (S). For each transcript i and each subgenome j in {CbpCg, CbpCo}, 
# S was computed as the subgenome relative expression deviation from the mean expression level in the parental
# species, u(i) = (E(ico) + E(icg)) / 2:
# S(ij) = (E(ij) - u(i)) / u(i)

library(stringr)
library(ggplot2)
library(ggpubr)

## Set work directory and load files
setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project/")

# This S index calculate the similarity of the expression of both subgenomes of Cbp to the average expression level
# of both parental species. 

parents <- c("CGCO", "CRCO", "CGCR") #In here, parental species could be CG&CO, CR&CO, and CG&CR.
tissue<-c("F","L","R") # set tissue

for(p in parents){ # Loop1 start
  message("Parental species are: ", p)
  for(t in tissue) { # Loop2 start
    Reads <- read.table(paste("InputData/AverageTMM_", t, ".txt", sep = ""), header = T, sep = "\t")
    message("dim TMM ", dim(Reads)[1] )
    CG <- "CG"
    CR <- "CR"
    CO <- "CO"
    Co <- "Co"
    Cg <- "Cg"
    CoColumns <- str_detect(names(Reads), Co) # detect which columns contain the "Co" character -- CbpCo
    Reads$meanCbpCo <- rowMeans(Reads[,CoColumns], na.rm = T) # these two lines must together, otherwise columns number change, logical values also changed
    CgColumns <- str_detect(names(Reads), Cg) 
    Reads$meanCbpCg <- rowMeans(Reads[,CgColumns], na.rm = T)
    
    if(p == "CGCO"){ # Loop3 start
      CGColumns <- str_detect(names(Reads), CG) # detect which columns contain the "CG" 
      COColumns <- str_detect(names(Reads), CO)
      Reads$meanParents <- rowMeans(Reads[, CGColumns | COColumns], na.rm = T)
    }else{
      if(p == "CRCO"){
        CRColumns <- str_detect(names(Reads), CR)
        COColumns <- str_detect(names(Reads), CO)
        Reads$meanParents <- rowMeans(Reads[, CRColumns | COColumns], na.rm = T)
      }else{
        CGColumns <- str_detect(names(Reads), CG) # detect which columns contain the "CG" 
        CRColumns <- str_detect(names(Reads), CR)
        Reads$meanParents <- rowMeans(Reads[, CGColumns | CRColumns], na.rm = T)
      }
    }# Loop3 end
    
    Reads$S_CbpCg <- (Reads$meanCbpCg - Reads$meanParents) / Reads$meanParents  # Subgenome Cg
    Reads$S_CbpCo <- (Reads$meanCbpCo - Reads$meanParents) / Reads$meanParents  # Subgenome Cg 

    # Loop4 for S value orientation
    if(p == "CGCO"){ # Loop4 start
      CGColumns <- str_detect(names(Reads), CG) # detect which columns contain the "CG" 
      COColumns <- str_detect(names(Reads), CO)
      orientate <- Reads[,CGColumns] > Reads[,COColumns]
    }else{
      if(p == "CRCO"){
        CRColumns <- str_detect(names(Reads), CR) 
        COColumns <- str_detect(names(Reads), CO)
        orientate <- Reads[,CRColumns] > Reads[,COColumns]
      }else{
        CGColumns <- str_detect(names(Reads), CG) 
        CRColumns <- str_detect(names(Reads), CR)
        orientate <- Reads[,CGColumns] > Reads[,CRColumns]
      }
    } # Loop4 end
    
    Reads[orientate,]$S_CbpCg <- Reads[orientate,]$S_CbpCg * (-1)
    Reads[orientate,]$S_CbpCo <- Reads[orientate,]$S_CbpCo * (-1)
    Reads$deltaS <- Reads$S_CbpCg + Reads$S_CbpCo
    
    # assign data
    assign(paste0("S_", t), Reads)
    
    # DE genes in different comparisons
    for (com in c("CGCO", "CGCR", "CRCO")){ # Loop5 start
      DE <- read.table(paste0("OutputData/DE_", t, "_", com, "_FDR.05.txt"), header = T, sep = "\t")
      
      # keep genes with CRCO FDR ≤ 0.05
      keep <- row.names(Reads) %in% row.names(DE)
      Reads_DE <- Reads[keep,]
      
      # assign data
      assign(paste0("S_", t, "_DE_",com), Reads_DE)
    } # loop5 end
    
  } # Loop2 end
  
  # Combine S data to a new dataframe for plot
  for (com in c("All", "CGCO", "CGCR", "CRCO")){ # Loop6 start
  S_case <- rep(c(rep(com,3)),3) # 3 tissues
  S_index <- c(rep(c("delta_S","S_Co", "S_Cg"), 3)) # 4 cases, 3 tissues (F,L,R)
  S_tissue <- c(rep("Flower", 3), rep("Leaf", 3), rep("Root", 3))
  
  if(com == "All"){
    S_value <- c(median(S_F$S_CbpCg + S_F$S_CbpCo), # Flower all
                 median(S_F$S_CbpCo),
                 median(S_F$S_CbpCg),
                 median(S_L$S_CbpCg + S_L$S_CbpCo), # Leaf all
                 median(S_L$S_CbpCo),
                 median(S_L$S_CbpCg),
                 median(S_R$S_CbpCg + S_R$S_CbpCo), # Root all
                 median(S_R$S_CbpCo),
                 median(S_R$S_CbpCg)
    )
  }else{
    if(com == "CGCO"){
      S_value <- c(median(S_F_DE_CGCO$S_CbpCg + S_F_DE_CGCO$S_CbpCo), # Flower all
                   median(S_F_DE_CGCO$S_CbpCo),
                   median(S_F_DE_CGCO$S_CbpCg),
                   median(S_L_DE_CGCO$S_CbpCg + S_L_DE_CGCO$S_CbpCo), # Leaf all
                   median(S_L_DE_CGCO$S_CbpCo),
                   median(S_L_DE_CGCO$S_CbpCg),
                   median(S_R_DE_CGCO$S_CbpCg + S_R_DE_CGCO$S_CbpCo), # Root all
                   median(S_R_DE_CGCO$S_CbpCo),
                   median(S_R_DE_CGCO$S_CbpCg)
      )
    }else{
      if (com == "CGCR"){
        S_value <- c(median(S_F_DE_CGCR$S_CbpCg + S_F_DE_CGCR$S_CbpCo), # Flower all
                     median(S_F_DE_CGCR$S_CbpCo),
                     median(S_F_DE_CGCR$S_CbpCg),
                     median(S_L_DE_CGCR$S_CbpCg + S_L_DE_CGCR$S_CbpCo), # Leaf all
                     median(S_L_DE_CGCR$S_CbpCo),
                     median(S_L_DE_CGCR$S_CbpCg),
                     median(S_R_DE_CGCR$S_CbpCg + S_R_DE_CGCR$S_CbpCo), # Root all
                     median(S_R_DE_CGCR$S_CbpCo),
                     median(S_R_DE_CGCR$S_CbpCg)
        )
      }else{
        S_value <- c(median(S_F_DE_CRCO$S_CbpCg + S_F_DE_CRCO$S_CbpCo), # Flower all
                     median(S_F_DE_CRCO$S_CbpCo),
                     median(S_F_DE_CRCO$S_CbpCg),
                     median(S_L_DE_CRCO$S_CbpCg + S_L_DE_CRCO$S_CbpCo), # Leaf all
                     median(S_L_DE_CRCO$S_CbpCo),
                     median(S_L_DE_CRCO$S_CbpCg),
                     median(S_R_DE_CRCO$S_CbpCg + S_R_DE_CRCO$S_CbpCo), # Root all
                     median(S_R_DE_CRCO$S_CbpCo),
                     median(S_R_DE_CRCO$S_CbpCg)
        )
      }
    }
  }
  
  # Combine all vectors to a dataframe
  S_plot <- cbind(S_case, S_index, S_tissue, S_value)
  # change the format of the dataframe
  S_plot <- print.data.frame(data.frame(S_plot), quote=FALSE)
  str(S_plot)
  S_plot$S_index <- factor(S_plot$S_index, levels = c("S_Co","S_Cg","delta_S"))
  S_plot$S_case <- factor(S_plot$S_case)
  S_plot$S_tissue <- factor(S_plot$S_tissue, levels = c("Flower","Leaf","Root"))
  S_plot$S_value <- as.numeric(S_plot$S_value)
  str(S_plot)
  S_plot$Comparison <- p
  names(S_plot) <- c("Case","Index","Tissue","Value","Comparison")
  write.table(S_plot, file = paste0("OutputData/ggplot2PlotData_for_S_index_", p,".txt"), col.names = T, row.names = F, sep = "\t")
  
  # assign data to a new name
  assign(paste0(p,"_DE_", com, "_S_plot") , S_plot)
  }# Loop6 end
  
} # Loop1 end


# Plot

ggplot(CGCO_DE_CGCO_S_plot, aes(x=Tissue, y=Value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=Index, shape=Index),size = 3, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  scale_y_continuous(limits=c(-0.15, 0.15)) +
  #scale_x_discrete(breaks = c('Flower', 'Leaf', 'Root'), 
  #                 labels = c('Flowers\nN=738', 'Leaves\nN=1108', 'Roots\nN=1715')) +
  facet_grid( ~ Case) + 
  labs(x ="Tissues", y = expression(median * " " * S[i])) +
  theme_bw()


ggplot(CGCO_DE_CGCO_S_plot, aes(x=Tissue, y=Value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=Index, shape=Index),size = 3, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  scale_y_continuous(limits=c(-0.3, 0.3)) +
  #scale_x_discrete(breaks = c('Flower', 'Leaf', 'Root'), 
  #                 labels = c('Flowers\nN=738', 'Leaves\nN=1108', 'Roots\nN=1715')) +
  facet_grid( ~ Case) + 
  labs(x ="Tissues", y = expression(median * " " * S[i])) +
  theme_bw()

head(S_plot)
p1 <- ggplot(CGCO_DE_CGCO_S_plot, aes(x=Tissue, y=Value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=Index, shape=Index),size = 4, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  scale_y_continuous(limits=c(-0.4, 0.4)) +
  #scale_x_discrete(breaks = c('Flower', 'Leaf', 'Root'), 
  #                 labels = c('Flowers\nN=738', 'Leaves\nN=1108', 'Roots\nN=1715')) +
  #facet_grid( ~ Comparison) + 
  labs(x ="Tissues", y = expression(median * " " * S[i])) +
  theme_bw()

p2 <- ggplot(CRCO_DE_CGCO_S_plot, aes(x=Tissue, y=Value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=Index, shape=Index),size = 4, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  scale_y_continuous(limits=c(-0.4, 0.4)) +
  #scale_x_discrete(breaks = c('Flower', 'Leaf', 'Root'), 
  #                 labels = c('Flowers\nN=738', 'Leaves\nN=1108', 'Roots\nN=1715')) +
  #facet_grid( ~ Comparison) + 
  labs(x ="Tissues", y = expression(median * " " * S[i])) +
  theme_bw()

p3 <- ggplot(CGCR_DE_CGCO_S_plot, aes(x=Tissue, y=Value)) +
  geom_hline(yintercept=0, linetype="dashed", color = "grey31") +
  geom_point(aes(fill=Index, shape=Index),size = 4, stroke = 1) + # stroke: making symbol bold
  scale_shape_manual(values=c(24, 25, 23))+
  scale_fill_manual(values=c('#56B4E9',"tomato", '#E69F00'))+
  scale_y_continuous(limits=c(-0.4, 0.4)) +
  #scale_x_discrete(breaks = c('Flower', 'Leaf', 'Root'), 
  #                 labels = c('Flowers\nN=738', 'Leaves\nN=1108', 'Roots\nN=1715')) +
  #facet_grid( ~ Comparison) + 
  labs(x ="Tissues", y = expression(median * " " * S[i])) +
  theme_bw()

# grid.newpage()
ggarrange(p1, p2, p3, ncol=3, nrow=1, common.legend = TRUE, legend="right") # set the common legend







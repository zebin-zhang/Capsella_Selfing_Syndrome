setwd("~/Desktop/My_Computer/Uppsala_University/Capsella_Project/")
# read PAC file from josephsetal
josephsetal <- read.table(file="connectivity/josephsetal_genedata.txt", header = T, sep = " ")
head(josephsetal)
# read PAC to gene file
PAC2gene <- read.table(file="connectivity/pac-to-carubv.txt", header = F, sep = "\t")
names(PAC2gene) <- c("PAC", "gene")
head(PAC2gene)
Cg_connect <- merge(PAC2gene, josephsetal, by="PAC") 

row.names(Cg_connect) <- Cg_connect$gene
head(Cg_connect)
row.names(Cg_connect) <- NULL
# write.table(Cg_connect, file = "connectivity/josephsetal_connectivity.txt", sep = "\t", quote = F, row.names = F)
# Load the coDEG files of pairwise comparison
F.mst <- read.table(file = "OutputData/coDEGs_CRCO_F.txt", sep = "\t", row.names = 1, header = T)
L.mst <- read.table(file = "OutputData/coDEGs_CRCO_L.txt", sep = "\t", row.names = 1, header = T)
R.mst <- read.table(file = "OutputData/coDEGs_CRCO_R.txt", sep = "\t", row.names = 1, header = T)
F.coGO <- read.table(file = "OutputData/coDEGs_CGCO_F.txt", sep = "\t", row.names = 1, header = T)
L.coGO <- read.table(file = "OutputData/coDEGs_CGCO_L.txt", sep = "\t", row.names = 1, header = T)
R.coGO <- read.table(file = "OutputData/coDEGs_CGCO_R.txt", sep = "\t", row.names = 1, header = T)
# F.coGR <- read.table(file = "OutputData/coDEG_CGCR_F.txt", sep = "\t", row.names = 1)
# L.coGR <- read.table(file = "OutputData/coDEG_CGCR_L.txt", sep = "\t", row.names = 1)
# R.coGR <- read.table(file = "OutputData/coDEG_CGCR_R.txt", sep = "\t", row.names = 1)

head(F.coGO)
head(F.mst)
'%!in%' <- function(x,y)!('%in%'(x,y))
# Distangle the unique and pleiotropy genes
for (co in c("coCRCO", "coCGCO")) { # co start
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
  co_F$gene <- row.names(co_F)
  co_F$DEGs <- co
  F.u.1 <- row.names(co_F) %!in% row.names(co_L) # flower not in leaf
  F.u.2 <- row.names(co_F) %!in% row.names(co_R) # flower not in root
  F.u.keep <- cbind(F.u.1 & F.u.2) # flower not in either
  row.names(co_F) <- NULL
  co_F$Status <- "pleiotropy"
  co_F[F.u.keep,]$Status <- "unique"
  assign(paste0("F_", co, ".UniPlei"), co_F)
  ## in leaf
  co_L$gene <- row.names(co_L)
  co_L$DEGs <- co
  L.u.1 <- row.names(co_L) %!in% row.names(co_F) # leaf not in flower
  L.u.2 <- row.names(co_L) %!in% row.names(co_R) # leaf not in root
  L.u.keep <- cbind(L.u.1 & L.u.2) # leaf not in either
  row.names(co_L) <- NULL
  co_L$Status <- "pleiotropy"
  co_L[L.u.keep,]$Status <- "unique"
  assign(paste0("L_", co, ".UniPlei"), co_L)
  ## in root
  co_R$gene <- row.names(co_R)
  co_R$DEGs <- co
  R.u.1 <- row.names(co_R) %!in% row.names(co_F) # root not in flower
  R.u.2 <- row.names(co_R) %!in% row.names(co_L) # root not in leaf
  R.u.keep <- cbind(R.u.1 & R.u.2) # root not in either
  row.names(co_R) <- NULL
  co_R$Status <- "pleiotropy"
  co_R[R.u.keep,]$Status <- "unique"
  assign(paste0("R_", co, ".UniPlei"), co_R)
}

# Merge coDEGs from different tissues and comparisons
head(F_coCRCO.UniPlei)
coDEGs <- rbind(F_coCRCO.UniPlei,
                L_coCRCO.UniPlei,
                R_coCRCO.UniPlei,
                F_coCGCO.UniPlei,
                L_coCGCO.UniPlei,
                R_coCGCO.UniPlei)
                #F_coCGCR.UniPlei,
                #L_coCGCR.UniPlei,
                #R_coCGCR.UniPlei)
head(coDEGs)

# coDEGs <- coDEGs[,c("gene","Tissue","DEGs","Status")]
# write.table(coDEGs, file = "OutputData/coDEGs_Unique_Pleiotropy_FLR.txt", sep = "\t", row.names = F, quote = F)
head(coDEGs)
# Merge coDEGs with connectivity
# coDEGs_CRCO <- rbind(coDEGs_CRCO_F,coDEGs_CRCO_L,coDEGs_CRCO_R)
# coDEGs_CRCO$coDEGs <- "CRCO"
# names(coDEGs_CRCO)[1] <- "gene"
# head(coDEGs_CRCO)
# 
# coDEGs_CGCO <- rbind(coDEGs_CGCO_F,coDEGs_CGCO_L,coDEGs_CGCO_R)
# coDEGs_CGCO$coDEGs <- "CGCO"
# names(coDEGs_CGCO)[1] <- "gene"
# head(coDEGs_CGCO)
# 
# coDEGs <- rbind(coDEGs_CRCO, coDEGs_CGCO)
# dim(coDEGs)
coDEGs_connect <- merge(coDEGs, Cg_connect, by = "gene", all = T)
head(coDEGs_connect)
# write.table(coDEGs_connect, file = "OutputData/coDEGs_connectivity_FLR.txt", sep = "\t", row.names = F, quote = F)

coDEGs_connect_T <- na.omit(coDEGs_connect) # Ture, without any NA value

head(coDEGs_connect_T)
dim(coDEGs_connect_T)
# pick random value of connectivity from non-coDEGs

random <- coDEGs_connect_T
random$DEGs <- "random"
head(random)
library(data.table)
set.seed(23)

Cg_connect <- data.table(Cg_connect) # important
random_connect <- Cg_connect[sample(.N, 1163)] # random select 4008 rows
head(random_connect)

random$sumConnec <- random_connect$sumConnec
dim(random)

coDEGs_connect_plot <- rbind(coDEGs_connect_T, random)
head(coDEGs_connect_plot)
coDEGs_connect_plot$ColorFill <- "random"
coDEGs_connect_plot[coDEGs_connect_plot$Tissue == "Flower" & coDEGs_connect_plot$DEGs == "coCRCO",]$ColorFill <- "F_CRCO"
coDEGs_connect_plot[coDEGs_connect_plot$Tissue == "Flower" & coDEGs_connect_plot$DEGs == "coCGCO",]$ColorFill <- "F_CGCO"
coDEGs_connect_plot[coDEGs_connect_plot$Tissue == "Leaf" & coDEGs_connect_plot$DEGs == "coCRCO",]$ColorFill <- "L_CRCO"
coDEGs_connect_plot[coDEGs_connect_plot$Tissue == "Leaf" & coDEGs_connect_plot$DEGs == "coCGCO",]$ColorFill <- "L_CGCO"
coDEGs_connect_plot[coDEGs_connect_plot$Tissue == "Root" & coDEGs_connect_plot$DEGs == "coCRCO",]$ColorFill <- "R_CRCO"
coDEGs_connect_plot[coDEGs_connect_plot$Tissue == "Root" & coDEGs_connect_plot$DEGs == "coCGCO",]$ColorFill <- "R_CGCO"

coDEGs_connect_plot$ColorFill <- factor(coDEGs_connect_plot$ColorFill, 
                                        levels = c("F_CRCO", "F_CGCO", "L_CRCO", "L_CGCO", "R_CRCO", "R_CGCO", "random"))
### plot
library(ggplot2)
coDEGs_connect_plot$Tissue <- factor(coDEGs_connect_plot$Tissue, levels = c("Flower","Leaf","Root"))
coDEGs_connect_plot$DEGs <- factor(coDEGs_connect_plot$DEGs, levels = c("coCRCO","coCGCO","random"))

mycomparison <- c("coCRCO", "coCGCO")
library(ggpubr)
compare_means(sumConnec ~ DEGs, data = coDEGs_connect_plot, 
              group.by = "Tissue", method = "wilcox",
              paired = F)


head(coDEGs_connect_plot)

ggplot(coDEGs_connect_plot,
       aes(Tissue, sumConnec)) +
  geom_boxplot(aes(fill = ColorFill, color = ColorFill), notch = TRUE) +
  facet_wrap(.~DE_1, ncol = 1) +
  scale_color_manual(values=c("black","#fa990e", "black", "#155800", "black", "#089400","black"), 
                     name = "coDEGs") +
  scale_fill_manual(values=c("#fa990e", "gray",  "#155800", "gray" ,"#089400","gray", "white"), 
                    name = "coDEGs") +
  #stat_compare_means(aes(group = coDEGs), method = "wilcox.test", label = "p.signif") +
  labs(#title = "Up-regulated coDEGs",
       x = "Tissues",
       y = "Connectivity") +
  theme_classic() 
 # theme(legend.position = "none") 
 

  
  #stat_compare_means(method = "wilcoxon") +
  #


p1 <- ggplot(coDEGs_connect_plot[coDEGs_connect_plot$DEGs == "coCRCO" | coDEGs_connect_plot$DEGs == "random",],
       aes(Tissue, sumConnec)) +
  geom_boxplot(aes(fill = Tissue), notch = TRUE) +
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(method = "wilcoxon") +
  scale_fill_manual(values=c("#fa990e", "#155800", "#089400","gray")) +
  labs(title = "coDEGs of CR & CO",
       x = "Tissues",
       y = "Connectivity") +
  theme_classic() +
  theme(legend.position = "none") 

p2 <- ggplot(coDEGs_connect_plot[coDEGs_connect_plot$DEGs == "coCGCO" | coDEGs_connect_plot$DEGs == "random",],
       aes(Tissue, sumConnec)) +
  geom_boxplot(aes(fill = Tissue), notch = TRUE) +
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(method = "wilcoxon") +
  scale_fill_manual(values=c("#fa990e", "#155800", "#089400","gray")) +
  labs(title = "coDEGs of CG & CO",
       x = "Tissues",
       y = "Connectivity") +
  theme_classic() +
  theme(legend.position = "none") 

p3 <- ggplot(coDEGs_connect_plot[coDEGs_connect_plot$DEGs == "coCGCR" | coDEGs_connect_plot$DEGs == "random",],
       aes(Tissue, sumConnec)) +
  geom_boxplot(aes(fill = Tissue), notch = TRUE) +
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(method = "wilcoxon") +
  scale_fill_manual(values=c("#fa990e", "#155800", "#089400","gray")) +
  labs(title = "coDEGs of CG & CR",
       x = "Tissues",
       y = "Connectivity") +
  theme_classic() +
  theme(legend.position = "none") 

library("grid")
library(ggthemes) # Load
grid.newpage()
pushViewport(viewport(layout = grid.layout(2,3)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p1, vp = vplayout(2,1))
print(p2, vp = vplayout(2,2))
print(p3, vp = vplayout(2,3))
print(venn_MST_p,vp = vplayout(1,1)) # from R code of EdgeR pairwise comparison
print(venn_coGO_p,vp = vplayout(1,2))
print(venn_coGR_p,vp = vplayout(1,3))










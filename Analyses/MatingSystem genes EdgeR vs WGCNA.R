# Mating system genes
setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project/OutputData")

# set function for notin
'%!in%' <- function(x,y)!('%in%'(x,y))

# Co expressed gene list of CR & CO from Wozniak etal Plant cell 2020
coCRCO <- read.csv("../InputData/Wozniak_etal_list_coDEG_all.csv")
head(coCRCO)
F.DE.CGCR <- read.table("DE_F_CGCR_FDR.05.txt", header = T)
F.DE.CGCO <- read.table("DE_F_CGCO_FDR.05.txt", header = T)
F.DE.CRCO <- read.table("DE_F_CRCO_FDR.05.txt", header = T)
head(F.DE.CRCO)
dim(F.DE.CGCO)
dim(F.DE.CRCO)

library("ggvenn")
MyVenn <- list(
  CG2CR = gsub(".g", "", row.names(F.DE.CGCR)),
  CG2CO = gsub(".g", "", row.names(F.DE.CGCO)),
  CR2CO = gsub(".g", "", row.names(F.DE.CRCO)),
  coCRCO = coCRCO$id
)


ggvenn(
  MyVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F
)

# Define MatingSystem related genes
# include in CGCR, CGCO but not in CRCO
F.filter1 <- row.names(F.DE.CGCR) %in% row.names(F.DE.CGCO)
F.filter2 <- row.names(F.DE.CGCR) %!in% row.names(F.DE.CRCO)
Fkeep <- cbind(F.filter1 & F.filter2)
F.mst <- F.DE.CGCR[Fkeep,]
dim(F.mst)
head(F.mst)
# write.csv(F.mst, "F_MatingSystemTransition_gene.csv", row.names = T, col.names = T, sep = "\t", quote = F)

L.DE.CGCR <- read.table("DE_L_CGCR_FDR.05.txt", header = T)
L.DE.CGCO <- read.table("DE_L_CGCO_FDR.05.txt", header = T)
L.DE.CRCO <- read.table("DE_L_CRCO_FDR.05.txt", header = T)
# Define MatingSystem related genes
# include in CGCR, CGCO but not in CRCO
L.filter1 <- row.names(L.DE.CGCR) %in% row.names(L.DE.CGCO)
L.filter2 <- row.names(L.DE.CGCR) %!in% row.names(L.DE.CRCO)
Lkeep <- cbind(L.filter1 & L.filter2)
L.mst <- L.DE.CGCR[Lkeep,]
dim(L.mst)
head(L.mst)
#write.csv(L.mst, "L_MatingSystemTransition_gene.csv", row.names = T, col.names = T, sep = "\t", quote = F)

R.DE.CGCR <- read.table("DE_R_CGCR_FDR.05.txt", header = T)
R.DE.CGCO <- read.table("DE_R_CGCO_FDR.05.txt", header = T)
R.DE.CRCO <- read.table("DE_R_CRCO_FDR.05.txt", header = T)
# Define MatingSystem related genes
# include in CGCR, CGCO but not in CRCO
R.filter1 <- row.names(R.DE.CGCR) %in% row.names(R.DE.CGCO)
R.filter2 <- row.names(R.DE.CGCR) %!in% row.names(R.DE.CRCO)
Rkeep <- cbind(R.filter1 & R.filter2)
R.mst <- R.DE.CGCR[Rkeep,]
dim(R.mst)
head(R.mst)
#write.csv(R.mst, "R_MatingSystemTransition_gene.csv", row.names = T, col.names = T, sep = "\t", quote = F)

dim(F.mst)
dim(L.mst)
dim(R.mst)
MSTfilter1 <- row.names(F.mst) %!in% row.names(L.mst)
MSTfilter2 <- row.names(F.mst) %!in% row.names(R.mst)
MSTkeep <- cbind(MSTfilter1 & MSTfilter2)
MST <- F.mst[MSTkeep,]
MST$Ortholog.gene.name <- row.names(MST)
head(MST)
dim(MST)
# cr2at <- read.table("../InputData/AT_CR_mart_export.txt", header=TRUE, sep = "\t")
# head(cr2at)
# dim(cr2at)

# MST <- merge(MST, cr2at, by = "Ortholog.gene.name")
# head(MST)
# dim(MST)

#write.csv(MST, "MatingSystemTransition_gene727.csv", row.names = T, col.names = T, sep = "\t", quote = F)


GreyGene <- read.table("WGCNA_gene_grey.txt", header = F)

TanjaLimma <- read.csv("/Users/zebinzhang/Desktop/My_Computer/Bioinformatics/R/TanjaNGData/limma_Slotte_etal_DE.csv", row.names = 1)
head(TanjaLimma)
TanjaPermuatation <- read.csv("/Users/zebinzhang/Desktop/My_Computer/Bioinformatics/R/TanjaNGData/premuatation_Slotte_etal_DE.csv", row.names = 1)
head(TanjaPermuatation)
Wozniak_etal <- read.csv("../InputData/Wozniak_etal_list_coDEG_all.csv")
head(Wozniak_etal)

Pink <- read.table("/Users/zebinzhang/Desktop/My_Computer/Bioinformatics/R/WGCNA/Capsella-Data/F_CRgene-pink.txt", header = F)
Purple <- read.table("/Users/zebinzhang/Desktop/My_Computer/Bioinformatics/R/WGCNA/Capsella-Data/F_CRgene-purple.txt", header = F)
Blue <- read.table("/Users/zebinzhang/Desktop/My_Computer/Bioinformatics/R/WGCNA/Capsella-Data/F_CRgene-blue.txt", header = F)

Lightgreen <- read.table("F_CRgene-lightgreen.txt", header = F)
Turquoise <- read.table("F_CRgene-turquoise.txt", header = F)
Tan <- read.table("F_CRgene-tan.txt", header = F)

library("ggvenn")
MyVenn <- list(
  EdgeR = gsub(".g", "", row.names(MST)),
  WGCNA_Pink = gsub(".g", "", Pink$V1),
  WGCNA_Purple = gsub(".g", "", Purple$V1),
  WGCNA_Blue = gsub(".g", "", Blue$V1)
)

MyVenn <- list(
  EdgeR = gsub(".g", "", row.names(MST)),
  WGCNA_Lightgreen = gsub(".g", "", Lightgreen$V1),
  WGCNA_Turquoise = gsub(".g", "", Turquoise$V1),
  WGCNA_Tan = gsub(".g", "", Tan$V1)
)
# without Cbp 
# ggvenn(
#   MyVenn, 
#   fill_color = c("#e48826", "#fd7c84", "#690fad", "#006beb"),
#   stroke_size = 0.5, set_name_size = 4,
#   show_percentage = F
# )

pdf("~/Desktop/My_Computer/Uppsala_University/Capsella_Project/Capsella_PDF/MST_GenesOverlap_EdgeR_WGCNA.pdf", width = 6.5, height = 6.5)

ggvenn(
  MyVenn, 
  fill_color = c("red", "green", "darkturquoise", "tan"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F
)

dev.off()

# without Cbp
#pinkKeep <- row.names(MST) %in% Pink$V1
#purpleKeep <- row.names(MST) %in% Purple$V1
#blueKeep <- row.names(MST) %in% Blue$V1

#wgcnaKeep <- cbind(pinkKeep | purpleKeep | blueKeep)

# with Cbp
Lightgreenkeep <- row.names(MST) %in% Lightgreen$V1
Turquoisekeep <- row.names(MST) %in% Turquoise$V1
Tankeep <- row.names(MST) %in% Tan$V1

wgcnaKeep <- cbind(Lightgreenkeep | Turquoisekeep | Tankeep)
RealMST <- MST[wgcnaKeep,]
head(RealMST)
dim(RealMST)

#write.csv(RealMST, "MST_genes_EdgeR_WGCNA_withCbp_482.csv", row.names = T, col.names = T, sep = "\t", quote = F)

MyVenn <- list(
  MST = gsub(".g", "", row.names(RealMST)), # MST genes detected in both EdgeR & WGCNA
  Slotte_etal_Limma = row.names(TanjaLimma),
  Slotte_etal_Permuatation = row.names(TanjaPermuatation),
  Wozniak_etal = Wozniak_etal$id
)


ggvenn(
  MyVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F
)

head(F.DE.CGCR)
head(GreyGene)


MyVenn <- list(
  Flowers = gsub(".g", "", row.names(F.mst)),
  Leaves = gsub(".g", "", row.names(L.mst)),
  Roots = gsub(".g", "", row.names(R.mst))
)


ggvenn(
  MyVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F
)


head(F.mst)




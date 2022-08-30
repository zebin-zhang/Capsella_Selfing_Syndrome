setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project/OutputData")

# Gene Ontology Analysis
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   +     install.packages("BiocManager")
# BiocManager::install("topGO")
# BiocManager::install("GOstats")

library(topGO)
library(GO.db)
library(GOstats)

# Mating system genes from EdgeR
# set function for notin
'%!in%' <- function(x,y)!('%in%'(x,y))
MST <- read.csv("MST_genes_EdgeR_WGCNA_WithCbp_428.csv", row.names = 1)
#ST$Ortholog.gene.name <- row.names(MST)
head(MST)
dim(MST)
# write.csv(MST, "MatingSystemTransition_gene.csv", row.names = T, col.names = T, sep = "\t", quote = F)

All <- read.table("TMM_CGCRCO_F.txt", header = T, sep = "\t")
head(All)
dim(All)
########################################################################
# # Method 1
# # Construct local R library for Capsella rubella genome wide annotation
# 
# # Biological Sequence Retrieval
# # BiocManager::install("biomartr")
# library(biomartr)
# 
# # check genome availability for Capsella rubella
# is.genome.available(db = "refseq", organism = "Capsella rubella", details = TRUE)
  # is.genome.available(db = "refseq", organism = "Arabidopsis lyrata", details = TRUE)
# # Making Capsella rubella packages for genome wide annotation
# library(AnnotationForge)
# # makeOrgPackageFromNCBI(version = "0.1",
# #                        author = "Zebin Zhang <zebin.zhang@ebc.uu.se>",
# #                        maintainer = "Zebin Zhang <zebin.zhang@ebc.uu.se>",
# #                        outputDir = ".",
# #                        tax_id = "81985",
# #                        genus = "Capsella",
# #                        species = "rubella")

# makeOrgPackageFromNCBI(version = "0.1",
#    author = "Zebin Zhang <zebin.zhang@ebc.uu.se>",
#    maintainer = "Zebin Zhang <zebin.zhang@ebc.uu.se>",
#    outputDir = ".",
#    tax_id = "81972",
#    genus = "Arabidopsis",
#    species = "lyrata")
# 
# # install.packages("./org.Crubella.eg.db", repos=NULL, type="source")
# library(org.Crubella.eg.db)
# 
# # convert gene name to Entrez
# MappingFile <- read.table("../InputData/Cru_id_mapping.txt", header = T, sep = "\t")
# head(MappingFile)
# MSTEntrezKeep <- MappingFile$Gene_id %in% row.names(MST)
# MSTEntrez <- MappingFile[MSTEntrezKeep,]
# head(MSTEntrez)
# dim(MSTEntrez)
# # remove unmatched / double matched rows
# MSTEntrez <- MSTEntrez[!(MSTEntrez$EntrezGene == "" | MSTEntrez$EntrezGene == "17881288, 17882126"), ]
# MSTEntrez <- MSTEntrez$EntrezGene
# head(MSTEntrez)
# 
# AllEntrezKeep <- MappingFile$Gene_id %in% row.names(All)
# AllEntrez <- MappingFile[AllEntrezKeep,]
# head(AllEntrez)
# dim(AllEntrez)
# AllEntrez <- AllEntrez$EntrezGene
# 
# # Perform Gene Ontology Analysis
# MyGO <- new("GOHyperGParams",
#             geneIds = MSTEntrez,
#             universeGeneIds = AllEntrez,
#             annotation = "org.Crubella.eg.db",
#             ontology = "BP",
#             pvalueCutoff = 0.05,
#             conditional = FALSE,
#             testDirection = "over")
# 
# MyGO_BP <- hyperGTest(MyGO)
# summary(MyGO_BP)[1:35,]
########################################################################

# since the annotation pathways are not good I try another method according to the Kryvokhyzha et al plos genetics 2019
# Method 2
# First download and load annotation file from http://plantregmap.gao-lab.org/download.php#go-annotation

# Step 1: Read mappings file
geneID2GO <- readMappings(file = "../Annotation/Cru_gene2go")
str(head(geneID2GO))
str(geneID2GO)
summary(geneID2GO)
# Step2: Create the background set
head(All)
AllID <- row.names(All)

# Step3: set interesting genes
head(MST)
mstID <- row.names(MST)

# Step4: Create the sample set
# This is the set of ID's in which we will look for enriched GO terms. 
# This set can be any subset of the background, and the set that is chosen affects the interpretation of the analysis.
geneList <- factor(as.integer(AllID %in% mstID))
names(geneList) <- AllID
str(geneList)

# Step5: Create GOData object
GO_BP <- new("topGOdata", 
              ontology = "BP", 
              allGenes = geneList,
              nodeSize = 5,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO)


# Step6: Run analysis
# This runs a Fisher's exact test
resultFis_BP <- runTest(GO_BP, algorithm = "classic", statistic = "fisher")

score(resultFis_BP)
# Step7: Look at the results
# This creates a plot, but the table allRes is perhaps more useful.
allRes_BP <- GenTable(GO_BP, 
                   classic = resultFis_BP,
                   orderBy = "weight", 
                   ranksOf = "classic", 
                   topNodes = 110)

GOBP_edger <- allRes_BP[allRes_BP$classic < 0.05,]
# output data
write.table(GOBP_edger, file = "482_MST_EdgeR_WGCNA_GO_BP_p0.05.txt", sep = "\t", row.names = F)

par(cex = 0.3) # smart way to increase the font size in plot
showSigOfNodes(GO_BP, score(resultFis_BP), firstSigNodes = 5, useInfo = 'all')


GO_MF <- new("topGOdata", 
             ontology = "MF", 
             allGenes = geneList,
             nodeSize = 5,
             annot = annFUN.gene2GO, 
             gene2GO = geneID2GO)

resultFis_MF <- runTest(GO_MF, algorithm = "classic", statistic = "fisher")

allRes_MF <- GenTable(GO_MF, 
                      classic = resultFis_MF,
                      orderBy = "weight", 
                      ranksOf = "classic", 
                      topNodes = 50)

GOMF_edger <- allRes_MF[allRes_MF$classic < 0.05,]

write.table(GOMF_edger, file = "482_MST_EdgeR_WGCNA_GO_MF_p0.05.txt", sep = "\t", row.names = F)

par(cex = 0.3)
showSigOfNodes(GO_MF, score(resultFis_MF), firstSigNodes = 5, useInfo = 'all')


######## WGCNA grey MST genes
# Step 1: Read mappings file
geneID2GO <- readMappings(file = "../Annotation/Cru_gene2go")
str(head(geneID2GO))

# Step2: Create the background set
head(All)
AllID <- row.names(All)

# Step3: set interesting genes
head(GreyWGCNA)
wgcnaID <- GreyWGCNA$V1

# Step4: Create the sample set
# This is the set of ID's in which we will look for enriched GO terms. 
# This set can be any subset of the background, and the set that is chosen affects the interpretation of the analysis.
geneList <- factor(as.integer(AllID %in% wgcnaID))
names(geneList) <- AllID

# Step5: Create GOData object
GO_BP <- new("topGOdata", 
             ontology = "BP", 
             allGenes = geneList,
             nodeSize = 5,
             annot = annFUN.gene2GO, 
             gene2GO = geneID2GO)


# Step6: Run analysis
# This runs a Fisher's exact test
resultFis_BP <- runTest(GO_BP, algorithm = "classic", statistic = "fisher")

score(resultFis_BP)
# Step7: Look at the results
# This creates a plot, but the table allRes is perhaps more useful.
allRes_BP <- GenTable(GO_BP, 
                      classic = resultFis_BP,
                      orderBy = "weight", 
                      ranksOf = "classic", 
                      topNodes = 50)

GOBP_wgcna <- allRes_BP[allRes_BP$classic < 0.05,]$GO.ID

# output data
write.table(allRes_BP, file = "WGCNAgrey_GO_BP_TOP50.txt", sep = "\t", row.names = F)

par(cex = 0.3) # smart way to increase the font size in plot
showSigOfNodes(GO_BP, score(resultFis_BP), firstSigNodes = 5, useInfo = 'all')


GO_MF <- new("topGOdata", 
             ontology = "MF", 
             allGenes = geneList,
             nodeSize = 5,
             annot = annFUN.gene2GO, 
             gene2GO = geneID2GO)

resultFis_MF <- runTest(GO_MF, algorithm = "classic", statistic = "fisher")

allRes_MF <- GenTable(GO_MF, 
                      classic = resultFis_MF,
                      orderBy = "weight", 
                      ranksOf = "classic", 
                      topNodes = 50)
GOMF_wgcna <- allRes_MF[allRes_MF$classic < 0.05,]$GO.ID


write.table(allRes_MF, file = "WGCNAgrey_GO_MF_TOP50.txt", sep = "\t", row.names = F)

par(cex = 0.3)
showSigOfNodes(GO_MF, score(resultFis_MF), firstSigNodes = 5, useInfo = 'all')


library("ggvenn")
BPVenn <- list(
  EdgeR = GOBP_edger,
  WGCNA = GOBP_wgcna
)


ggvenn(
  BPVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F
)

BPOverlap <- GOBP_edger %in% GOBP_wgcna

BPGO <- GOBP_edger[BPOverlap]
## MF
MFVenn <- list(
  EdgeR = GOMF_edger,
  WGCNA = GOMF_wgcna
)


ggvenn(
  MFVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F
)

MFOverlap <- GOMF_edger %in% GOMF_wgcna

MFGO <- GOMF_edger[MFOverlap]


## GO analysis on overlapping genes of 482 MST genes and genes from two other references (Slotte et al & Wozniak et al)
TanjaLimma <- read.csv("/Users/zebinzhang/Desktop/My_Computer/Bioinformatics/R/TanjaNGData/limma_Slotte_etal_DE.csv", row.names = 1)
head(TanjaLimma)
TanjaPremuatation <- read.csv("/Users/zebinzhang/Desktop/My_Computer/Bioinformatics/R/TanjaNGData/premuatation_Slotte_etal_DE.csv", row.names = 1)
head(TanjaPremuatation)
Wozniak_etal <- read.csv("../InputData/Wozniak_etal_list_coDEG_all.csv")
head(Wozniak_etal)
head(MST)
dim(MST)


library("ggvenn")

MyVenn <- list(
  MST = gsub(".g", "", row.names(MST)), # MST genes detected in both EdgeR & WGCNA
  Slotte_etal_Limma = row.names(TanjaLimma),
  Slotte_etal_Premuatation = row.names(TanjaPremuatation),
  Wozniak_etal = Wozniak_etal$id
)


ggvenn(
  MyVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F
)

head(MST)
dim(MST)
MST$Slotte <- "-"
MST$Woznia <- "-"

overlap_1 <- gsub(".g", "", row.names(MST)) %in% row.names(TanjaLimma)
overlap_2 <- gsub(".g", "", row.names(MST)) %in% row.names(TanjaPremuatation)
overlap_3 <- gsub(".g", "", row.names(MST)) %in% Wozniak_etal$id
overlap_keep <- cbind(overlap_1 | overlap_2 | overlap_3)

MST[overlap_1,]$Slotte <- "+"
MST[overlap_2,]$Slotte <- "+"
MST[overlap_3,]$Woznia <- "+"
head(MST)
# subset 218 overlapping genes.a
Gene_218 <- MST[overlap_keep,]
head(Gene_218)
dim(Gene_218)

write.table(MST, file = "482_MST_genes_overlap_Slotte_Wozniak.txt", sep = "\t", row.names = T, quote = F)
write.table(Gene_218, file = "218_MST_genes_overlap_Slotte_Wozniak.txt", sep = "\t", row.names = T, quote = F)

# Step 1: Read mappings file
geneID2GO <- readMappings(file = "../Annotation/Cru_gene2go")
str(head(geneID2GO))

# Step2: Create the background set
head(All)
AllID <- row.names(All)

# Step3: set interesting genes
head(Gene_218)
Gene_218_ID <- row.names(Gene_218)

# Step4: Create the sample set
# This is the set of ID's in which we will look for enriched GO terms. 
# This set can be any subset of the background, and the set that is chosen affects the interpretation of the analysis.
geneList <- factor(as.integer(AllID %in% Gene_218_ID))
names(geneList) <- AllID

# Step5: Create GOData object
GO_BP <- new("topGOdata", 
             ontology = "BP", 
             allGenes = geneList,
             nodeSize = 5,
             annot = annFUN.gene2GO, 
             gene2GO = geneID2GO)


# Step6: Run analysis
# This runs a Fisher's exact test
resultFis_BP <- runTest(GO_BP, algorithm = "classic", statistic = "fisher")

score(resultFis_BP)
# Step7: Look at the results
# This creates a plot, but the table allRes is perhaps more useful.
allRes_BP <- GenTable(GO_BP, 
                      classic = resultFis_BP,
                      orderBy = "weight", 
                      ranksOf = "classic", 
                      topNodes = 50)

GOBP_gene_218 <- allRes_BP[allRes_BP$classic < 0.05,]$GO.ID

# output data
write.table(allRes_BP, file = "Gene_218_GO_BP_TOP50.txt", sep = "\t", row.names = F)

par(cex = 0.3) # smart way to increase the font size in plot
showSigOfNodes(GO_BP, score(resultFis_BP), firstSigNodes = 5, useInfo = 'all')


GO_MF <- new("topGOdata", 
             ontology = "MF", 
             allGenes = geneList,
             nodeSize = 5,
             annot = annFUN.gene2GO, 
             gene2GO = geneID2GO)

resultFis_MF <- runTest(GO_MF, algorithm = "classic", statistic = "fisher")

allRes_MF <- GenTable(GO_MF, 
                      classic = resultFis_MF,
                      orderBy = "weight", 
                      ranksOf = "classic", 
                      topNodes = 50)
GOMF_gene_218 <- allRes_MF[allRes_MF$classic < 0.05,]$GO.ID


write.table(allRes_MF, file = "Gene_218_GO_MF_TOP50.txt", sep = "\t", row.names = F)

par(cex = 0.3)
showSigOfNodes(GO_MF, score(resultFis_MF), firstSigNodes = 5, useInfo = 'all')


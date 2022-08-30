#####################################################################################
#
#              Heatmap for all gene expressions across CG,CR,CO
#
#####################################################################################

# R Packages/functions for drawing heatmaps
# 
# There are a multiple numbers of R packages and functions for drawing interactive and static heatmaps, including:
#   
# heatmap() [R base function, stats package]: Draws a simple heatmap
# heatmap.2() [gplots R package]: Draws an enhanced heatmap compared to the R base function.
# pheatmap() [pheatmap R package]: Draws pretty heatmaps and provides more control to change the appearance of heatmaps.
# d3heatmap() [d3heatmap R package]: Draws an interactive/clickable heatmap
# Heatmap() [ComplexHeatmap R/Bioconductor package]: Draws, annotates and arranges complex heatmaps (very useful for genomic data analysis)

setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project")

# install packages
# library(devtools)
# install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library("gplots")
My_data <- read.table("InputData/Capsella_Diploids_ReadsCount_TMM_Filtered.txt", header=TRUE)
head(My_data)
dim(My_data)
TMM <- My_data[ ,2:37] # only kept columns of CGCRCO 
TMM <- TMM[ ,c(5:8,1:4,9:12,17:20,13:16,21:24,29:32,25:28,33:36)]
head(TMM)
names(TMM)

TMM <- as.matrix(TMM)

MyColors <- cbind(Species=rep(c(rep("#2176ff",4), rep("#33a1fd",4), rep("#f79824",4)), 3),
                  Tissues=c(rep("#ef476f",12), rep("#06d6a0",12), rep("#ffd166",12))
)

FlowerColors<-cbind(Species=c(rep("#2176ff",4), rep("#33a1fd",4), rep("#f79824",4)),
                    Tissues=rep("#ef476f",12)
                    )

LeafColors<-cbind(Species=c(rep("#2176ff",4), rep("#33a1fd",4), rep("#f79824",4)),
                    Tissues=rep("#06d6a0",12)
)

RootColors<-cbind(Species=c(rep("#2176ff",4), rep("#33a1fd",4), rep("#f79824",4)),
                    Tissues=rep("#ffd166",12)
)


#####################################################################
# Method 2 Complex Heatmap
scaleRow <- function(x) {
  rm <- rowMeans(x, na.rm = T)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd, na.rm = T)
  x <- sweep(x, 1, sx, "/")
  return(round(x, 6))
}
scaleTMM <- scaleRow(TMM)
head(scaleTMM)

#install.packages("magick")
library(magick)
library(dendextend)
row_dend = as.dendrogram(hclust(dist(scaleTMM))) # generate hierarchical dendrogram for rows
# row_dend = color_branches(row_dend ) # `color_branches()` returns a dendrogram object
column_dend = as.dendrogram(hclust(dist(t(TMM)))) # generate hierarchical dendrogram for rows - t transpose data
# column_dend = color_branches(column_dend, k = 3) # `color_branches()` returns a dendrogram object
#  Different method try
png("~/Desktop/ComplexHeatmap.png", width=1000, height=1000)
Heatmap(scaleTMM, 
        name = "TMM", 
        cluster_rows = row_dend,
        cluster_columns = column_dend,
        column_order = sort(colnames(TMM)),
        row_dend_width = unit(4, "cm"),
        column_dend_height = unit(4, "cm"),
        row_split = 3,
        column_split = 3
        )
dev.off()

# Heatmap for flower
TMM_F <- TMM[,1:12]
head(TMM_F)

scaleTMM_F <- scaleRow(TMM_F)
head(scaleTMM_F)
F_row_dend = as.dendrogram(hclust(dist(scaleTMM_F))) # generate hierarchical dendrogram for rows
#F_row_dend = color_branches(F_row_dend, k = 3 ) # `color_branches()` returns a dendrogram object
F_column_dend = as.dendrogram(hclust(dist(t(TMM_F)))) # generate hierarchical dendrogram for rows - t transpose data
F_column_dend = color_branches(F_column_dend, k = 2) # `color_branches()` returns a dendrogram object

# Add annotation in different colours

MyColors <- cbind(Species=rep(c(rep("#2176ff",4), rep("#33a1fd",4), rep("#f79824",4)), 3),
                  Tissues=c(rep("#ef476f",12), rep("#06d6a0",12), rep("#ffd166",12))
)

F_annotation = HeatmapAnnotation(
  Tissues = rep("Flower", 12),
  Species = c(rep("CG",4), rep("CR",4), rep("CO",4)),
  col = list(Tissues=c("Flower" = "#fa990e", "Leaf" = "#155800", "Root" = "#089400"),
             Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
  ),
  # annotation_name_side = "left",       # set annotation label to left
  show_annotation_name = FALSE         # NOT show annotation label
)

Heatmap(scaleTMM_F, 
        name = "TMM", 
        cluster_rows = F_row_dend,
        cluster_columns = F_column_dend,
        column_order = sort(colnames(TMM_F)),
        row_dend_width = unit(4, "cm"),
        column_dend_height = unit(4, "cm"),
        row_split = 3,
        column_split = 2,
        top_annotation = F_annotation
) 


# Heatmap for LEAF
TMM_L <- TMM[,13:24]
head(TMM_L)

scaleTMM_L <- scaleRow(TMM_L)
head(scaleTMM_L)
L_row_dend = as.dendrogram(hclust(dist(scaleTMM_L))) # generate hierarchical dendrogram for rows
#L_row_dend = color_branches(L_row_dend, k = 3 ) # `color_branches()` returns a dendrogram object
L_column_dend = as.dendrogram(hclust(dist(t(TMM_L)))) # generate hierarchical dendrogram for rows - t transpose data
L_column_dend = color_branches(L_column_dend, k = 2) # `color_branches()` returns a dendrogram object

L_annotation = HeatmapAnnotation(
  Tissues = rep("Leaf", 12),
  Species = c(rep("CG",4), rep("CR",4), rep("CO",4)),
  col = list(Tissues=c("Flower" = "#fa990e", "Leaf" = "#155800", "Root" = "#089400"),
             Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
  ),
  # annotation_name_side = "left",       # set annotation label to left
  show_annotation_name = FALSE         # NOT show annotation label
)

Heatmap(scaleTMM_L, 
        name = "TMM", 
        cluster_rows = L_row_dend,
        cluster_columns = L_column_dend,
        column_order = sort(colnames(TMM_L)),
        row_dend_width = unit(4, "cm"),
        column_dend_height = unit(4, "cm"),
        row_split = 3,
        column_split = 2,
        top_annotation = L_annotation
)

# Heatmap for Root
TMM_R <- TMM[,25:36]
head(TMM_R)

scaleTMM_R <- scaleRow(TMM_R)
head(scaleTMM_R)
R_row_dend = as.dendrogram(hclust(dist(scaleTMM_R))) # generate hierarchical dendrogram for rows
#R_row_dend = color_branches(R_row_dend, k = 3 ) # `color_branches()` returns a dendrogram object
R_column_dend = as.dendrogram(hclust(dist(t(TMM_R)))) # generate hierarchical dendrogram for rows - t transpose data
R_column_dend = color_branches(R_column_dend, k = 2) # `color_branches()` returns a dendrogram object

R_annotation = HeatmapAnnotation(
  Tissues = rep("Root", 12),
  Species = c(rep("CG",4), rep("CR",4), rep("CO",4)),
  col = list(Tissues=c("Flower" = "#fa990e", "Leaf" = "#155800", "Root" = "#089400"),
             Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
  ),
  # annotation_name_side = "left",       # set annotation label to left
  show_annotation_name = FALSE         # NOT show annotation label
)

Heatmap(scaleTMM_R, 
        name = "TMM", 
        cluster_rows = R_row_dend,
        cluster_columns = R_column_dend,
        column_order = sort(colnames(TMM_R)),
        row_dend_width = unit(4, "cm"),
        column_dend_height = unit(4, "cm"),
        row_split = 3,
        column_split = 2,
        top_annotation = R_annotation
)


# Only plot legend
library(circlize)
# col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
# lgd1 = Legend(col_fun = col_fun, title = "TMM")
# lgd2 = Legend(at = c("Flower","Leaf","Root"), 
#               legend_gp = gpar(fill = c("#ef476f","#06d6a0","#ffd166")),
#               title = "Tissues")
# lgd3 = Legend(at = c("CG","CR","CO"), 
#               legend_gp = gpar(fill = c("#2176ff","#33a1fd","#f79824")),
#               title = "Species")
# pd = packLegend(lgd1, lgd2, lgd3)
# 
# draw(pd,  x = unit(0.1, "cm"), y = unit(0.5, "npc"), just="left") 

# Plot all figures
#pdf("Capsella_PDF/TMM Heatmap by tissues 10x6.pdf", width=10, height=6)

library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun, title = "TMM")
lgd2 = Legend(at = c("Flower","Leaf","Root"), 
              legend_gp = gpar(fill = c("#fa990e","#155800","#089400")),
              title = "Tissues")
lgd3 = Legend(at = c("CG","CR","CO"), 
              legend_gp = gpar(fill = c("#f90052","darkviolet","#006beb")),
              title = "Species")
pd = packLegend(lgd1, lgd2, lgd3)

draw(pd,  x = unit(0.1, "cm"), y = unit(0.5, "npc"), just="left") 

png("Capsella_PDF/TMM Heatmap by tissues 10x6.png", width=2000, height=1500, res = 300) # resolution=300
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 10)))
# Draw flower
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:3)) # set the region for plot
draw(Heatmap(scaleTMM_F,                                        # using scaled data
             name = "TMM", 
             cluster_rows = F_row_dend,                         # replace default dendrogram to custom one
             cluster_columns = F_column_dend,
             #column_order = sort(colnames(TMM_F)),             # reorder columns by column names
             show_row_names = FALSE,                         # hide row names
             show_column_names = FALSE,
             row_title = NULL,
             column_title = NULL,
             row_dend_width = unit(0.5, "cm"),
             column_dend_height = unit(1, "cm"),
             row_split = 3,
             column_split = 2,
             top_annotation = F_annotation                      # add column annotation, ie, colour of species and tissues
             ),
     show_annotation_legend = FALSE,                            # hide the annotation legend (tissues and species)
     show_heatmap_legend = FALSE,                               # hide the heatmap legend (mapping value color red2blue)
     newpage = FALSE)                                           # don't create a new page for plot

upViewport()

# Draw leaf
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4:6))
draw(Heatmap(scaleTMM_L, 
             name = "TMM", 
             cluster_rows = L_row_dend,
             cluster_columns = L_column_dend,
             #column_order = sort(colnames(TMM_L)),
             show_row_names = FALSE,                         # hide row names
             show_column_names = FALSE,
             row_title = NULL,
             column_title = NULL,
             row_dend_width = unit(0.5, "cm"),
             column_dend_height = unit(1, "cm"),
             row_split = 3,
             column_split = 2,
             top_annotation = L_annotation
             ),
     show_heatmap_legend = FALSE,
     show_annotation_legend = FALSE,
     #annotation_legend_side = "bottom",
     newpage = FALSE
)
upViewport()
# Draw Root
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 7:9))
draw(Heatmap(scaleTMM_R, 
             name = "TMM", 
             cluster_rows = R_row_dend,
             cluster_columns = R_column_dend,
             #column_order = sort(colnames(TMM_R)),
             show_row_names = FALSE,                         # hide row names
             show_column_names = FALSE,
             row_title = NULL,
             column_title = NULL,
             row_dend_width = unit(0.5, "cm"),
             column_dend_height = unit(1, "cm"),
             row_split = 3,
             column_split = 2,
             top_annotation = R_annotation
),
show_heatmap_legend = FALSE,
show_annotation_legend = FALSE,
#annotation_legend_side = "bottom",
newpage = FALSE
)
upViewport()
# Draw Legend
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 10))

draw(pd,  x = unit(0.01, "cm"), 
     y = unit(0.5, "npc"), 
     just="left") 

upViewport()

upViewport()

dev.off()

#######################
# use average value of each species 
TMM <- read.table("InputData/Diploids_individual_TMM_FLR.txt", header=TRUE)
head(TMM)
TMM <- TMM[ ,c(5:8,1:4,9:12,17:20,13:16,21:24,29:32,25:28,33:36)]
head(TMM)
names(TMM)
##### define popultaion
CG_F <- c("CG1_F","CG2_F","CG3_F","CG4_F")
CR_F <- c("CR1_F","CR2_F","CR3_F","CR4_F")
CO_F <- c("CO1_F","CO2_F","CO3_F","CO4_F")
CG_L <- c("CG1_L","CG2_L","CG3_L","CG4_L")
CR_L <- c("CR1_L","CR2_L","CR3_L","CR4_L")
CO_L <- c("CO1_L","CO2_L","CO3_L","CO4_L")
CG_R <- c("CG1_R","CG2_R","CG3_R","CG4_R")
CR_R <- c("CR1_R","CR2_R","CR3_R","CR4_R")
CO_R <- c("CO1_R","CO2_R","CO3_R","CO4_R")

head(TMM[,CR_F])
# Use the MEAN value of each sample to represents as population value
# Flower of Cr, Cg, Co
TMM$CG_F <- rowMeans(TMM[,CG_F], na.rm = T)
TMM$CR_F <- rowMeans(TMM[,CR_F], na.rm = T)
TMM$CO_F <- rowMeans(TMM[,CO_F], na.rm = T)
# Leaf of Cr, Cg, Co
TMM$CG_L <- rowMeans(TMM[,CG_L], na.rm = T)
TMM$CR_L <- rowMeans(TMM[,CR_L], na.rm = T)
TMM$CO_L <- rowMeans(TMM[,CO_L], na.rm = T)
# Root of Cr, Cg, Co
TMM$CG_R <- rowMeans(TMM[,CG_R], na.rm = T)
TMM$CR_R <- rowMeans(TMM[,CR_R], na.rm = T)
TMM$CO_R <- rowMeans(TMM[,CO_R], na.rm = T)

# Split mean value to a new dataframe
MeanTMM <- TMM[,37:45]
head(MeanTMM)
MeanTMM <- as.matrix(MeanTMM)

scaleRow <- function(x) {
  rm <- rowMeans(x, na.rm = T)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd, na.rm = T)
  x <- sweep(x, 1, sx, "/")
  return(round(x, 6))
}
scaleMeanTMM <- scaleRow(MeanTMM)
head(scaleMeanTMM)

#install.packages("magick")
library(magick)
library(dendextend)
row_dend_Mean = as.dendrogram(hclust(dist(scaleMeanTMM))) # generate hierarchical dendrogram for rows
row_dend_Mean = color_branches(row_dend_Mean, k = 3) # `color_branches()` returns a dendrogram object
column_dend_Mean = as.dendrogram(hclust(dist(t(MeanTMM)))) # generate hierarchical dendrogram for rows - t transpose data
column_dend_Mean = color_branches(column_dend_Mean, k = 3) # `color_branches()` returns a dendrogram object

library(ComplexHeatmap)

#  Different method try
# png("~/Desktop/ComplexHeatmap.png", width=1000, height=1000)
Heatmap(scaleMeanTMM, 
        name = "MeanTMM", 
        cluster_rows = row_dend_Mean,
        cluster_columns = column_dend_Mean,
        #column_order = sort(colnames(meanTMM)),
        show_row_names = FALSE,                         # hide row names
        #show_column_names = FALSE,                         # hide column names
        row_dend_width = unit(4, "cm"),
        column_dend_height = unit(4, "cm"),
        row_split = 3,
        column_split = 3
)
# dev.off()

# Heatmap for flower
head(MeanTMM)
MeanTMM_F <- MeanTMM[,1:3]
head(MeanTMM_F)
# write.table(MeanTMM_F, file = "OutputData/MeanTMM_F.txt", sep = "\t", quote = F)
scaleMeanTMM_F <- scaleRow(MeanTMM_F)
head(scaleMeanTMM_F)
Mean_F_row_dend = as.dendrogram(hclust(dist(scaleMeanTMM_F))) # generate hierarchical dendrogram for rows
# Mean_F_row_dend = color_branches(Mean_F_row_dend, k = 4 ) # `color_branches()` returns a dendrogram object
Mean_F_column_dend = as.dendrogram(hclust(dist(t(MeanTMM_F)))) # generate hierarchical dendrogram for rows - t transpose data
Mean_F_column_dend = color_branches(Mean_F_column_dend, k = 2) # `color_branches()` returns a dendrogram object
Mean_F_annotation = HeatmapAnnotation(
  Tissues = rep("Flower", 3),
  Species = c("CG","CR","CO"),
  col = list(Tissues=c("Flower" = "#fa990e", "Leaf" = "#155800", "Root" = "#089400"),
             Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
  ),
  # annotation_name_side = "left",       # set annotation label to left
  show_annotation_name = FALSE         # NOT show annotation label
)
# Leaf
MeanTMM_L <- MeanTMM[,4:6]
head(MeanTMM_L)
# write.table(MeanTMM_L, file = "OutputData/MeanTMM_L.txt", sep = "\t", quote = F)

scaleMeanTMM_L <- scaleRow(MeanTMM_L)
head(scaleMeanTMM_L)
Mean_L_row_dend = as.dendrogram(hclust(dist(scaleMeanTMM_L))) # generate hierarchical dendrogram for rows
#Mean_L_row_dend = color_branches(Mean_L_row_dend, k = 4 ) # `color_branches()` returns a dendrogram object
Mean_L_column_dend = as.dendrogram(hclust(dist(t(MeanTMM_L)))) # generate hierarchical dendrogram for rows - t transpose data
Mean_L_column_dend = color_branches(Mean_L_column_dend, k = 2) # `color_branches()` returns a dendrogram object
Mean_L_annotation = HeatmapAnnotation(
  Tissues = rep("Leaf", 3),
  Species = c("CG","CR","CO"),
  col = list(Tissues=c("Flower" = "#fa990e", "Leaf" = "#155800", "Root" = "#089400"),
             Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
  ),
  # annotation_name_side = "left",       # set annotation label to left
  show_annotation_name = FALSE         # NOT show annotation label
)

# Root
MeanTMM_R <- MeanTMM[,7:9]
head(MeanTMM_R)
# write.table(MeanTMM_R, file = "OutputData/MeanTMM_R.txt", sep = "\t", quote = F)

scaleMeanTMM_R <- scaleRow(MeanTMM_R)
head(scaleMeanTMM_R)
Mean_R_row_dend = as.dendrogram(hclust(dist(scaleMeanTMM_R))) # generate hierarchical dendrogram for rows
# Mean_R_row_dend = color_branches(Mean_R_row_dend, k = 4 ) # `color_branches()` returns a dendrogram object
Mean_R_column_dend = as.dendrogram(hclust(dist(t(MeanTMM_R)))) # generate hierarchical dendrogram for rows - t transpose data
Mean_R_column_dend = color_branches(Mean_R_column_dend, k = 2) # `color_branches()` returns a dendrogram object
Mean_R_annotation = HeatmapAnnotation(
  Tissues = rep("Root", 3),
  Species = c("CG","CR","CO"),
  col = list(Tissues=c("Flower" = "#fa990e", "Leaf" = "#155800", "Root" = "#089400"),
             Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
  ),
  # annotation_name_side = "left",       # set annotation label to left
  show_annotation_name = FALSE         # NOT show annotation label
)
# Only plot legend
library(circlize)
col_fun = colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
lgd1 = Legend(col_fun = col_fun, title = "TMM")
lgd2 = Legend(at = c("Flower","Leaf","Root"), 
              legend_gp = gpar(fill = c("#fa990e","#155800","#089400")),
              title = "Tissues")
lgd3 = Legend(at = c("CG","CR","CO"), 
              legend_gp = gpar(fill = c("#f90052","darkviolet","#006beb")),
              title = "Species")
pd = packLegend(lgd1, lgd2, lgd3)

draw(pd,  x = unit(0.1, "cm"), y = unit(0.5, "npc"), just="left") 

# Plot all figures
# pdf("Capsella_PDF/Fig 1 Heatmap 6x6.pdf", width=6.1, height=6)
png("Capsella_PDF/Fig1 Newheatmap 7x5.png",
  width     = 7,
  height    = 5,
  units     = "in",
  res       = 1200,
  pointsize = 4
)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 1, nc = 10)))
# Draw flower
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1:3)) # set the region for plot

draw(Heatmap(scaleMeanTMM_R, 
             name = "TMM", 
             cluster_rows = Mean_R_row_dend,
             cluster_columns = Mean_R_column_dend,
             show_row_names = FALSE,                         # hide row names
             show_column_names = FALSE,
             row_title = NULL,
             column_title = NULL,
             row_dend_width = unit(0.5, "cm"),
             column_dend_height = unit(1, "cm"),
             row_split = 5,
             column_split = 2,
             top_annotation = Mean_R_annotation
),
show_heatmap_legend = FALSE,
show_annotation_legend = FALSE,
#annotation_legend_side = "bottom",
newpage = FALSE
)

upViewport()
# Draw leaf
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4:6))
draw(Heatmap(scaleMeanTMM_L, 
             name = "TMM", 
             cluster_rows = Mean_L_row_dend,
             cluster_columns = Mean_L_column_dend,
             show_row_names = FALSE,                         # hide row names
             show_column_names = FALSE,
             row_title = NULL,
             column_title = NULL,
             row_dend_width = unit(0.5, "cm"),
             column_dend_height = unit(1, "cm"),
             row_split = 5,
             column_split = 2,
             top_annotation = Mean_L_annotation
),
show_heatmap_legend = FALSE,
show_annotation_legend = FALSE,
#annotation_legend_side = "bottom",
newpage = FALSE
)
upViewport()
# Draw Root
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 7:9))

draw(Heatmap(scaleMeanTMM_F,                                        # using scaled data
             name = "TMM", 
             cluster_rows = Mean_F_row_dend,                         # replace default dendrogram to custom one
             cluster_columns = Mean_F_column_dend,
             #column_order = sort(colnames(TMM_F)),             # reorder columns by column names
             show_row_names = FALSE,                         # hide row names
             show_column_names = FALSE,                         # hide column names
             row_title = NULL,                                  # hide cluster order numbers
             column_title = NULL,
             row_dend_width = unit(0.5, "cm"),                    # length of dendrogram
             column_dend_height = unit(1, "cm"),
             row_split = 5,                                     # split row dendrogram to 3 by hierarchical relationship 
             column_split = 2,
             top_annotation = Mean_F_annotation                      # add column annotation, ie, colour of species and tissues
),
show_annotation_legend = FALSE,                            # hide the annotation legend (tissues and species)
show_heatmap_legend = FALSE,                               # hide the heatmap legend (mapping value color red2blue)
newpage = FALSE)                                           # don't create a new page for plot

upViewport()
# Draw Legend
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 10))

draw(pd,  x = unit(0.01, "cm"), 
     y = unit(0.5, "npc"), 
     just="left") 

upViewport()

upViewport()

dev.off()

# Heatmap for Selfing syndrome genes. 
MST <- read.csv("OutputData/MST_genes_EdgeR_WGCNA_withCbp_482.csv", row.names = 1)
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
  col = list(Tissues=c("Flower" = "#fa990e", "Leaf" = "#155800", "Root" = "#089400"),
             Species=c("CG" = "#f90052", "CR" = "darkviolet", "CO" = "#006beb")
  ),
  # annotation_name_side = "left",       # set annotation label to left
  show_annotation_name = FALSE         # NOT show annotation label
)

png("Capsella_PDF/MSTgenes_482_Newheatmap_4x6.png",
    width     = 4,
    height    = 6,
    units     = "in",
    res       = 1200,
    pointsize = 4
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
             top_annotation = Mean_MST_F_annotation                      # add column annotation, ie, colour of species and tissues
),
#show_annotation_legend = FALSE,                            # hide the annotation legend (tissues and species)
#show_heatmap_legend = FALSE,                               # hide the heatmap legend (mapping value color red2blue)
newpage = FALSE
)                                           # don't create a new page for plot

dev.off()







#######################
# Another methods
library(heatmap3)
png("~/Desktop/heatmap3_2.png", width=1000, height=1000)

heatmap3(TMM, 
         scale = "row",                     # Z scale the data by row
         showRowDendro = T,                 # NOT show the row dendrogram
         labRow = F,                        # NOT show the row names
         #labCol = F,                      # NOT show the row names
         ColSideColors=MyColors,
         col=colorRampPalette(c("#4B0082","#4B0082","white", "#228B22", "#228B22","darkgreen","darkgreen"))(100)
)
dev.off()

# The dendrogram was generated by raw TMM data WITHOUT Z-score transformation, however the mapping colour was performed
# by Z transformated data.

png("~/Desktop/Flower heatmap3.png", width=1000, height=1000)

heatmap3(TMM[,1:12], 
         scale = "row",                     # Z scale the data by row
         showRowDendro = T,                 # NOT show the row dendrogram
         labRow = F,                        # NOT show the row names
         labCol = F,                      # NOT show the row names
         ColSideColors=FlowerColors,
         col=colorRampPalette(c("#4B0082","#4B0082","white","#228B22","#228B22"))(1024)
)
dev.off()

png("~/Desktop/Leaf heatmap3.png", width=1000, height=1000)

heatmap3(TMM[,13:24], 
         scale = "row",                     # Z scale the data by row
         showRowDendro = F,                 # NOT show the row dendrogram
         labRow = F,                        # NOT show the row names
         labCol = F,                      # NOT show the row names
         ColSideColors=LeafColors,
         col=colorRampPalette(c("blue","gray34", "yellow"))(1024)
)
dev.off()

png("~/Desktop/Root heatmap3.png", width=1000, height=1000)

heatmap3(TMM[,25:36], 
         scale = "row",                     # Z scale the data by row
         showRowDendro = F,                 # NOT show the row dendrogram
         labRow = F,                        # NOT show the row names
         labCol = F,                      # NOT show the row names
         ColSideColors=RootColors,
         col=colorRampPalette(c("blue","gray34", "yellow"))(1024)
)
dev.off()

p = heatmap.2(cpm, dendrogram="both", scale="row")

p # Outputs all the data in the list; lots of output to the console

str(p) # Struture of p; also lots of output to the console

names(p) # Names of all the list elements

p$rowInd # Ordering of the data rows

p$carpet # The heatmap values

heatmap2_value <- t(p$carpet)
# sording dataframe by CG1_F
heatmap2_value <- as.data.frame(heatmap2_value)
heatmap2_value <- heatmap2_value[,c(
  "CG1_F","CG2_F","CG3_F","CG4_F",
  "CR1_F","CR2_F","CR3_F","CR4_F",
  "CO1_F","CO2_F","CO3_F","CO4_F",
  "CG1_L","CG2_L","CG3_L","CG4_L",
  "CR1_L","CR2_L","CR3_L","CR4_L",
  "CO1_L","CO2_L","CO3_L","CO4_L",
  "CG1_R","CG2_R","CG3_R","CG4_R",
  "CR1_R","CR2_R","CR3_R","CR4_R",
  "CO1_R","CO2_R","CO3_R","CO4_R"
  )]

sort_heatmap2_value <- heatmap2_value[order(heatmap2_value$CG1_F),]
sort_heatmap2_value <- round(sort_heatmap2_value, 6)
rownames(sort_heatmap2_value) <- NULL
colnames(sort_heatmap2_value) <- NULL
head(sort_heatmap2_value)
str(sort_heatmap2_value)

manual_value <- as.data.frame(scale1)
sort_manual_value <- manual_value[order(manual_value$CG1_F),]
rownames(sort_manual_value) <- NULL
colnames(sort_manual_value) <- NULL
head(sort_manual_value)
str(sort_manual_value)

all((sort_manual_value[] == sort_heatmap2_value) == TRUE)

require(dplyr)
anti_join(sort_manual_value,sort_heatmap2_value)


install.packages("compare")
library(compare)
comparison <- compare(sort_manual_value, sort_heatmap2_value, allowAll=TRUE)
head(comparison)
comparison$tM

heatmap.2(cpm, 
          scale = "row",           #scale by row
          col = bluered(100), 
          trace = "none", 
          dendrogram="column",
          density.info = "none")

heatmap.2(as.matrix(p$carpet), 
          scale = "none",           #scale by row
          col = bluered(100), 
          trace = "none", 
          dendrogram="both",
          density.info = "none")

# heatmap.2 default scaling
scaleRow <- function(x) {
  rm <- rowMeans(x, na.rm = T)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd, na.rm = T)
  x <- sweep(x, 1, sx, "/")
  return(round(x, 6))
}
scale1 <- scaleRow(cpm)
head(scale1)

# Manual scaling
heatmap_scale_rows = function (x) 
{
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m)/s)
}

scale2 <- heatmap_scale_rows(cpm)

# compare two scales
all((scale1 == round(scale2, 6)) == TRUE)

  
heatmap.2(scale1, 
          scale = "none",           # without scale
          col = bluered(100), 
          trace = "none", 
          dendrogram="column",
          density.info = "none")

heatmap.2(scale2,                  
          scale = "none",           # without scale
          col = bluered(100), 
          trace = "none", 
          dendrogram="column",
          density.info = "none")

# heatmap.2(scale1,
#           Rowv=F,
#           Colv=F,
#           dendrogram="none") 

library("ggplot2")
# install.packages("ggdendro")
library("ggdendro")
library("reshape2")

cpm_gg <- My_data[ ,2:37] # only kept columns of CGCRCO 
cpm_gg <- cpm_gg[ ,c(5:8,1:4,9:12,17:20,13:16,21:24,29:32,25:28,33:36)]
row.names(cpm_gg) <- My_data$Gene
cpm_gg <- as.matrix(cpm_gg)
head(cpm_gg)

head(t(cpm_gg)[,1:10])
dim(t(cpm_gg))
str(t(cpm_gg))

species_scaled <- scale(t(cpm_gg))
head(species_scaled[,1:10])

species.dendro <- as.dendrogram(hclust(d = dist(x = species_scaled)))
ggdendrogram(data = species.dendro)




cpm_gg[,1:36] <- scale1
cpm_gg <- cpm_gg[,c(37,1:36)]
head(cpm_gg)
gg <- melt(cpm_gg)
head(gg)

pdf("~/Desktop/ggplot2 heatmap.pdf", width=8, height=10)
ggplot(gg, aes(variable, gene, fill= value)) + 
  geom_tile()+
  scale_fill_gradient2()
dev.off()




############################################################
#     ggplot2 method

# Compute the correlation matrix
cormat <- round(cor(cpm),2)
head(cormat)
# Create the correlation heatmap with ggplot2
library(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)

library(ggplot2)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# Get the lower and upper triangles of the correlation matrix
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
# usage
upper_tri <- get_upper_tri(cormat)
upper_tri

# Finished correlation matrix heatmap
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", mid = "blue", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()





heatmap.2(scale_cpm, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")

all.equal(scale1,scale_cpm, check.attributes = FALSE)

library("pheatmap")


head(cpm)
head(scale_cpm)
pheatmap(scale_cpm)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("made4")
library(made4)
heatplot(scale_cpm,scale="row",margins=c(4,7))


# Different clustering methods
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)

Heatmap(t(scale(t(cpm[,1:12]))), 
        name = " ",
        #cluster_columns=FALSE, # turn off column clustering
        #show_row_dend = FALSE, # hide row dendrogram
        clustering_distance_rows ="euclidean",
        clustering_method_rows = "ward.D",
        #cluster_columns = agnes,
        row_dend_width = unit(3, "cm"), # Make the column dendrogram wider
        column_dend_height = unit(4, "cm") # Make the column dendrogram wider
)
# Z score scaling
library(matrixStats)
# create a function of rowscaling
# https://www.r-bloggers.com/2016/02/a-faster-scale-function/
scale_matrix_rows = function(x,
                    center = TRUE,
                    scale = TRUE,
                    add_attr = F,
                    rows = NULL,
                    cols = NULL) {
  if (!is.null(rows) && !is.null(cols)) {
    x <- x[rows, cols, drop = FALSE]
  } else if (!is.null(rows)) {
    x <- x[rows, , drop = FALSE]
  } else if (!is.null(cols)) {
    x <- x[, cols, drop = FALSE]
  }
  ################
  # Get the column means
  ################
  cm = rowMeans(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
  if (scale) {
    csd = rowSds(x, center = cm)
  } else {
    # just divide by 1 if not
    csd = rep(1, length = length(cm))
  }
  if (!center) {
    # just subtract 0
    cm = rep(0, length = length(cm))
  }
  x = (x - cm) / csd
  if (add_attr) {
    if (center) {
      attr(x, "scaled:center") <- cm
    }
    if (scale) {
      attr(x, "scaled:scale") <- csd
    }
  }
  return(x)
}


Z_cpm <- scale_matrix_rows(
  cpm,
  center = TRUE,
  scale = TRUE,
  add_attr = FALSE,
  rows = NULL,
  cols = NULL
)

head(Z_cpm)

Z_cpm_matrix <- as.matrix(Z_cpm) #  convert to matrix

# Keep Sample info
id_info <- data.frame(id = colnames(Z_cpm) )
id_info

# Make heatmap# Different distance calculation methods
# "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
# euclidean is the default

# Different clustering methods
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC)
library(cluster)
# pdf("~/Desktop/My_Computer/Uppsala_University/Capsella_Project/Capsella_PDF/heatmap of all genes 9x10.pdf", width=8, height=10)
# grid.newpage()
Heatmap(Z_cpm_matrix, 
        name = " ",
        #cluster_columns=FALSE, # turn off column clustering
        show_row_dend = FALSE, # hide row dendrogram
        #clustering_distance_rows ="kendall",
        #clustering_method_rows = "ward.D",
        #cluster_columns = agnes,
        row_dend_width = unit(3, "cm"), # Make the column dendrogram wider
        column_dend_height = unit(4, "cm") # Make the column dendrogram wider
        )
dev.off()

##  method 2
# install.packages("gplots")
install.packages("heatmap.plus")
# install.packages("RColorBrewer")

library("gplots")
library("heatmap.plus")
library("RColorBrewer")




if (!require("NMF")) {
  install.packages("NMF", dependencies = TRUE)
  library(NMF)
}

pdf("~/Desktop/NMF heatmap.pdf", width=8, height=10)
aheatmap(cpm, 
         color = rev(brewer.pal(9,"RdBu")), 
         scale="row",
         annColors = "Set2")
dev.off()


head(scale1)


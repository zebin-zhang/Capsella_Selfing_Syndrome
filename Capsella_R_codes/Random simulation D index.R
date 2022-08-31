#  D index pairwise comparison among Cg, Cr, and Co.
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

#####################################################################################################
# Random Simulation
listOfDataFrames <- vector(mode = "list", length = 1000)  # Create a new list for all simulated data

for (i in 1:1000) {                    # set 10,000 simulations
  set.seed(i)                           # set seed for each simulation for obtain certain values 
  if (i == 1){
    Value <- rnorm(17396, 0, 0.2)       # generate randome values with mean=0 and sd=0.2
    Random <- data.frame(Value)
    Random$Rep <- paste0("rep",i)
    dim(Random)
    listOfDataFrames[[i]] <- data.frame(Random)
    Median <- median(Value)            # Calculate the median value for each simulation
  }else{
    Value <- rnorm(17396, 0, 0.2)
    App <- data.frame(Value)
    App$Rep <- paste0("rep",i)
    listOfDataFrames[[i]] <- data.frame(App)
    Median[i] <- median(Value)
  }
}

Random <- do.call("rbind", listOfDataFrames) # combine all simulation by do.call

head(Random)
dim(Random)
str(Random)
Random$Rep <- factor(Random$Rep)      # as factor
Random_sd <- sd(Random$Value)
greyScale <- rep(gray.colors(10, start = 0, end = 0.7), 100)     # scale grey color to 10 scales from dark to light, repeat 100 times.

library(ggplot2)

# Plot simulation result
p1 <- ggplot(Random, aes(x=Value, color = Rep))+
  geom_density() + 
  scale_color_manual(values=greyScale) +   
  xlim(-1,1) +
  theme_classic() +                                 # set theme to classic
  theme(legend.position = "none") +                   # remove legend
  labs( #title = "Fuel economy declines as weight increases",
    #subtitle = "(1973-74)",
    #caption = "Data from the 1974 Motor Trend US magazine.",
    tag = "D",
    x = " ",
    y = "Density"
    #colour = "Gears"
  )


# Plot D index distribution in F L R.
# create new dataframe for F, L R
Value <- CPMnorm.Pop.Comp$F.D  # From D index pairwise comparison across Cg Cr Co.R
D_F <- data.frame(Value)
D_F$Rep <- "Flower"

Value <- CPMnorm.Pop.Comp$L.D
D_L <- data.frame(Value)
D_L$Rep <- "Leaf"

Value <- CPMnorm.Pop.Comp$R.D
D_R <- data.frame(Value)
D_R$Rep <- "Root"

D <- rbind(D_F, D_L, D_R)
str(D)
D$Rep <- factor(D$Rep, levels = c("Flower","Leaf","Root"))
D$Class <- "D = 0"
Flower <- D$Rep == "Flower"
F_sd <- sd(D[Flower,]$Value)
Leaf <- D$Rep == "Leaf"
L_sd <- sd(D[Leaf,]$Value)
Root <- D$Rep == "Root"
R_sd <- sd(D[Root,]$Value)

head(D)

D[Flower & D$Value > F_sd,]$Class <- "D > 0"
D[Flower & D$Value < -F_sd,]$Class <- "D < 0"
D[Leaf & D$Value > L_sd,]$Class <- "D > 0"
D[Leaf & D$Value < -L_sd,]$Class <- "D < 0"
D[Root & D$Value > R_sd,]$Class <- "D > 0"
D[Root & D$Value < -R_sd,]$Class <- "D < 0"


F_density <- with(density(D[Flower,]$Value), data.frame(x, y)) # calculate the density for flower and convert to dataframe

adjustcolor("firebrick1", alpha.f = 0.3) # generate transparented color

p2 <- ggplot(data = F_density, mapping = aes(x = x, y = y)) +
  geom_line(size=1 )+
  geom_area(mapping = aes(x = ifelse(x > F_sd , x, 0)), fill = "#ffcecc", size=1) + #light red
  geom_area(mapping = aes(x = ifelse(x > -F_sd & x < F_sd, x, 0)), fill = "white") + 
  geom_area(mapping = aes(x = ifelse(x < -F_sd , x, -F_sd)), fill = "gray94") + 
  geom_vline(aes(xintercept=0),
             color="gray34", linetype="dashed", size=1) +
  geom_vline(aes(xintercept = median(D[Flower,]$Value)),
             color="firebrick1", linetype="dashed", size=1) + #dark red color
  xlim(-1,1) +
  ylim(0,max(F_density$y)) +
  annotate(geom="text", x=-1, y=2, label=paste("median =", round(median(D[Flower,]$Value),2)) , color = "black",hjust=0) + # add annotation for median  fontface =2, bold
  annotate(geom="text", x=-1, y=1.8, label= expression(paste(sigma, " = 0.26")), color = "black", fontface =2, hjust=0) + # standard divation
  theme_classic() +
  labs(tag = "A",
       x = "D index",
       y = "Density"
  )


L_density <- with(density(D[Leaf,]$Value), data.frame(x, y)) # calculate the density for flower and convert to dataframe

adjustcolor("forestgreen", alpha.f = 0.3) # generate transparented color

p3 <- ggplot(data = L_density, mapping = aes(x = x, y = y)) +
  geom_line(size=1)+
  geom_area(mapping = aes(x = ifelse(x > L_sd , x, 0)), fill = "#d5f3d3") + # light green
  geom_area(mapping = aes(x = ifelse(x > -L_sd & x < L_sd, x, 0)), fill = "white") + 
  geom_area( mapping = aes(x = ifelse(x < -L_sd , x, -L_sd)), fill = "gray94") + 
  geom_vline(aes(xintercept=0),
             color="gray34", linetype="dashed", size=1) +
  geom_vline(aes(xintercept = median(D[Leaf,]$Value)),
             color="forestgreen", linetype="dashed", size=1) + # dark green
  xlim(-1,1) +
  ylim(0,max(F_density$y)) +
  annotate(geom="text", x=-1, y=2, label=paste("median =", round(median(D[Leaf,]$Value),2)) , color = "black",  hjust=0) + # add annotation for median
  annotate(geom="text", x=-1, y=1.8, label= expression(paste(sigma, " = 0.31")), color = "black", fontface =2, hjust=0) + # standard divation
  theme_classic() +
  labs(tag = "B",
       x = "D index",
       y = "Density"
  )

adjustcolor("purple", alpha.f = 0.3)

R_density <- with(density(D[Root,]$Value), data.frame(x, y)) # calculate the density for flower and convert to dataframe

p4 <- ggplot(data = R_density, mapping = aes(x = x, y = y)) +
  geom_line( size=1)+
  geom_area(mapping = aes(x = ifelse(x > R_sd , x, 0)), fill = "#f4d8b3") + 
  geom_area(mapping = aes(x = ifelse(x > -R_sd & x < R_sd, x, 0)), fill = "white") + 
  geom_area(mapping = aes(x = ifelse(x < -R_sd , x, -R_sd)), fill = "gray94") + 
  geom_vline(aes(xintercept=0),
             color="gray34", linetype="dashed", size=1) +
  geom_vline(aes(xintercept = median(D[Root,]$Value)),
             color="#da7702", linetype="dashed", size=1) +
  xlim(-1,1) +
  ylim(0,max(F_density$y)) +
  annotate(geom="text", x=-1, y=2, label=paste("median =", round(median(D[Root,]$Value),2)) , color = "black",  hjust=0) + # add annotation for median
  annotate(geom="text", x=-1, y=1.8, label= expression(paste(sigma, " = 0.30")), color = "black", fontface =2, hjust=0) + # standard divation
  theme_classic() +
  labs(tag = "C",
       x = "D index",
       y = "Density"
  )




# p2 <- ggplot(D, aes(x=Value, color = Rep, fill = Rep))+
#   geom_density(alpha=0.2, size=1) +                         # set transparent rate
#   scale_color_manual(values=c("firebrick1","forestgreen","purple4")) +       
#   scale_fill_manual(values= c("sienna1", "lightseagreen",  "orchid1")) +  
#   xlim(-1,1) +
#   theme_classic() +                                 # set theme to classic
#   theme(legend.position = "none") +                  # remove legend
#   labs(tag = "B",
#     x = "D index",
#     y = "Density"
#   )

D_median <- c(median(D_F$Value), median(D_L$Value), median(D_R$Value))
D_median <- c(D_median, Median) 
D_median <- data.frame(D_median)
D_median$Rep <- c("Flower","Leaf","Root", rep("simulation",i))

for(i in 1:1000){
  D_median[3+i,]$Rep <- paste0("S",i)
}
D_median$y <- c(2,2,2, rep(1,i))
head(D_median)
str(D_median)
D_median1 <- D_median
D_median1$y <- 0
head(D_median1)
D_median <- rbind(D_median1, D_median)

greyScale <- rep(gray.colors(10, start = 0, end = 0.7), 1000)
colorScale <- c("firebrick1","forestgreen","#da7702",greyScale)

p5 <- ggplot(D_median,aes(x=D_median,y=y, group=Rep))+
  geom_line(aes(color=Rep), size=1) +
  scale_color_manual(values=colorScale) +  
  theme_classic() +
  theme(legend.position = "none",
        panel.grid.minor.y  = element_blank(),
        #axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(), 
        #axis.line.y = element_blank(),  # remove y axis line
        panel.background = element_blank()) +
  labs(tag = "E",
       x = "Median D index")
  
#####################################################################################################
F_sd2 <- sd(CPMnorm.Pop.Comp$F.D) * 2
L_sd2 <- sd(CPMnorm.Pop.Comp$L.D) * 2
R_sd2 <- sd(CPMnorm.Pop.Comp$R.D) * 2

# Statistic how many genes with D value bigger than 0 + 2sd

dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$F.D > F_sd2,])[1]
dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$L.D > L_sd2,])[1]
dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$R.D > L_sd2,])[1]
message("In Flower, The number of genes with higher differences between outcrosser and selfer than differnces within selfer are: ", dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$F.D > F_sd,])[1])
message("In Leaf, The number of genes with higher differences between outcrosser and selfer than differnces within selfer are: ", dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$L.D > L_sd,])[1])
message("In Root, The number of genes with higher differences between outcrosser and selfer than differnces within selfer are: ", dim(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$R.D > L_sd,])[1])

venn_D_FLR <- list(
  Flower = row.names(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$F.D > F_sd2,]),
  Leaf = row.names(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$L.D > L_sd2,]),
  Root = row.names(CPMnorm.Pop.Comp[CPMnorm.Pop.Comp$R.D > L_sd2,])
)

names(venn_D_FLR) <- c("Flowers","Leaves","Roots")

p6 <- ggvenn(
  venn_D_FLR, 
  fill_color = c("firebrick1","forestgreen","#da7702"),
  stroke_size = 0.3, set_name_size = 4,
  show_percentage = F
) +
  labs(tag = "F")

library("grid")
# install.packages("ggthemes") # Install 
library(ggthemes) # Load

pdf("~/Desktop/My_Computer/Uppsala_University/Capsella_Project/Capsella_PDF/D index with 1000 simulations1 8x5.pdf", width=8, height=5)
# Cbp in all populations
grid.newpage()
pushViewport(viewport(layout = grid.layout(4,6)))
vplayout <- function(x,y){viewport(layout.pos.row = x, layout.pos.col = y)}  
print(p2, vp = vplayout(1:2,1:2))
print(p3, vp = vplayout(1:2,3:4))
print(p4, vp = vplayout(1:2,5:6))
print(p1, vp = vplayout(3:4,1:2))
print(p5, vp = vplayout(4,3:4))
print(p6, vp = vplayout(3:4,5:6))

dev.off()





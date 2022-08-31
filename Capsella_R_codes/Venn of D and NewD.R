setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project/OutputData")

# set function for notin
'%!in%' <- function(x,y)!('%in%'(x,y))
library("ggvenn")
library(ggplot2)

F_DE_CGCO_D <- read.table("DE_F_CGCO_FDR.05_D_index.txt", header = T)
head(F_DE_CGCO_D)
F_DE_CGCO_D_Pos <- F_DE_CGCO_D[F_DE_CGCO_D$Value > 0, ]
F_DE_CGCO_D_Neg <- F_DE_CGCO_D[F_DE_CGCO_D$Value < 0, ]

F_DE_CGCR_New_D <- read.table("DE_F_CGCR_FDR.05_New_D_index.txt", header = T)
head(F_DE_CGCR_New_D)
F_DE_CGCR_New_D_Pos <- F_DE_CGCR_New_D[F_DE_CGCR_New_D$NewIndex > 0, ]
F_DE_CGCR_New_D_Neg <- F_DE_CGCR_New_D[F_DE_CGCR_New_D$NewIndex < 0, ]

Merge_F <- merge(F_DE_CGCO_D,F_DE_CGCR_New_D,by="row.names",all = T)
head(Merge_F)
Merge_F <- Merge_F[,c("Row.names","Value","NewIndex")]
names(Merge_F) <- c("Gene","D","NewD")
#write.table(Merge_F, "D_NewD_F.txt", append = F, sep = "\t")

L_DE_CGCO_D <- read.table("DE_L_CGCO_FDR.05_D_index.txt", header = T)
head(L_DE_CGCO_D)
L_DE_CGCO_D_Pos <- L_DE_CGCO_D[L_DE_CGCO_D$Value > 0, ]
L_DE_CGCO_D_Neg <- L_DE_CGCO_D[L_DE_CGCO_D$Value < 0, ]

L_DE_CGCR_New_D <- read.table("DE_L_CGCR_FDR.05_New_D_index.txt", header = T)
head(L_DE_CGCR_New_D)
L_DE_CGCR_New_D_Pos <- L_DE_CGCR_New_D[L_DE_CGCR_New_D$NewIndex > 0, ]
L_DE_CGCR_New_D_Neg <- L_DE_CGCR_New_D[L_DE_CGCR_New_D$NewIndex < 0, ]

Merge_L <- merge(L_DE_CGCO_D,L_DE_CGCR_New_D,by="row.names",all = T)
head(Merge_L)
Merge_L <- Merge_L[,c("Row.names","Value","NewIndex")]
names(Merge_L) <- c("Gene","D","NewD")

R_DE_CGCO_D <- read.table("DE_R_CGCO_FDR.05_D_index.txt", header = T)
head(R_DE_CGCO_D)
R_DE_CGCO_D_Pos <- R_DE_CGCO_D[R_DE_CGCO_D$Value > 0, ]
R_DE_CGCO_D_Neg <- R_DE_CGCO_D[R_DE_CGCO_D$Value < 0, ]

R_DE_CGCR_New_D <- read.table("DE_R_CGCR_FDR.05_New_D_index.txt", header = T)
head(R_DE_CGCR_New_D)
R_DE_CGCR_New_D_Pos <- R_DE_CGCR_New_D[R_DE_CGCR_New_D$NewIndex > 0, ]
R_DE_CGCR_New_D_Neg <- R_DE_CGCR_New_D[R_DE_CGCR_New_D$NewIndex < 0, ]

Merge_R <- merge(R_DE_CGCO_D,R_DE_CGCR_New_D,by="row.names",all = T)
head(Merge_R)
Merge_R <- Merge_R[,c("Row.names","Value","NewIndex")]
names(Merge_R) <- c("Gene","D","NewD")

F_NegVenn <- list(
  NegCGCO = gsub(".g", "", row.names(F_DE_CGCO_D_Neg)),
  NegCGCR = gsub(".g", "", row.names(F_DE_CGCR_New_D_Neg))
)

p1 <- ggvenn(
  F_NegVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = T
)

F_PosVenn <- list(
  PosCGCO = gsub(".g", "", row.names(F_DE_CGCO_D_Pos)),
  PosCGCR = gsub(".g", "", row.names(F_DE_CGCR_New_D_Pos))
)

p2 <- ggvenn(
  F_PosVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = T
)

L_NegVenn <- list(
  NegCGCO = gsub(".g", "", row.names(L_DE_CGCO_D_Neg)),
  NegCGCR = gsub(".g", "", row.names(L_DE_CGCR_New_D_Neg))
)

p3 <- ggvenn(
  L_NegVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = T
)

L_PosVenn <- list(
  PosCGCO = gsub(".g", "", row.names(L_DE_CGCO_D_Pos)),
  PosCGCR = gsub(".g", "", row.names(L_DE_CGCR_New_D_Pos))
)

p4 <- ggvenn(
  L_PosVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = T
)

R_NegVenn <- list(
  NegCGCO = gsub(".g", "", row.names(R_DE_CGCO_D_Neg)),
  NegCGCR = gsub(".g", "", row.names(R_DE_CGCR_New_D_Neg))
)

p5 <- ggvenn(
  R_NegVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = T
)

R_PosVenn <- list(
  PosCGCO = gsub(".g", "", row.names(R_DE_CGCO_D_Pos)),
  PosCGCR = gsub(".g", "", row.names(R_DE_CGCR_New_D_Pos))
)

p6 <- ggvenn(
  R_PosVenn, 
  #fill_color = c("#1B9E77", "#D95F02", "#7570B3"),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = T
)

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 2, nrow = 3)




for (tissue in c("F","L","R")) {
  Percentage <- c()  # Create a new vector for all simulated data
  Percentage2 <- c()
  mystart <- -1.0
  myend <- mystart + 0.1
  
  if (tissue == "F"){
    mytissue <- "Flower"
    D_Neg <- F_DE_CGCO_D_Neg
    NewD_Neg <- F_DE_CGCR_New_D_Neg
    D_Pos <- F_DE_CGCO_D_Pos
    NewD_Pos <- F_DE_CGCR_New_D_Pos
  }else{
    if (tissue == "L"){
      mytissue <- "Leaf"
      D_Neg <- L_DE_CGCO_D_Neg
      NewD_Neg <- L_DE_CGCR_New_D_Neg
      D_Pos <- L_DE_CGCO_D_Pos
      NewD_Pos <- L_DE_CGCR_New_D_Pos
    }else{
      mytissue <- "Root"
      D_Neg <- R_DE_CGCO_D_Neg
      NewD_Neg <- R_DE_CGCR_New_D_Neg
      D_Pos <- R_DE_CGCO_D_Pos
      NewD_Pos <- R_DE_CGCR_New_D_Pos
    }
  }
  
  for (i in 1:20) {
    if(myend <= 0){
      D <- D_Neg
      NewD <- NewD_Neg
      Dsub <- D[D$Value >= mystart & D$Value < myend,]
      overlap <- row.names(Dsub) %in% row.names(NewD)
      Keep <- Dsub[overlap,]
      percent <- dim(Keep)[1] / (dim(Dsub)[1] + dim(NewD)[1] - dim(Keep)[1])
      percent1 <- dim(Keep)[1] / dim(Dsub)[1] 
      Percentage[i] <- percent * 100
      Percentage2[i] <- percent1 * 100
    }else{
      D <- D_Pos
      NewD <- NewD_Pos
      Dsub <- D[D$Value >= mystart & D$Value < myend,]
      overlap <- row.names(Dsub) %in% row.names(NewD)
      Keep <- Dsub[overlap,]
      percent <- dim(Keep)[1] / (dim(Dsub)[1] + dim(NewD)[1] - dim(Keep)[1])
      percent1 <- dim(Keep)[1] / dim(Dsub)[1]
      Percentage[i] <- percent * 100
      Percentage2[i] <- percent1 * 100
    }
    mystart <- mystart + 0.1
    myend <- myend + 0.1
  }
  
  Range <- seq(from = -0.95, to = 0.95, by =0.1)
  Tissue <- rep(mytissue, 20)
  myData <- data.frame(Range, Percentage, Percentage2, Tissue)
  
  assign(paste0("Percentage_", tissue), myData)
}

Percentage_FLR <- rbind(Percentage_F, Percentage_L, Percentage_R)

ggplot(Percentage_FLR, aes(Range, Percentage, colour = Tissue)) +
  geom_point() + geom_line()
  #geom_step()
ggplot(Percentage_FLR, aes(Range, Percentage2, colour = Tissue)) +
  geom_point() + geom_line()


MaxOverlap <- rbind(Percentage_F, Percentage_L, Percentage_R)
MaxOverlap$Direction <- rep(c(rep("Negative", 10), rep("Positive", 10)), 3)
head(MaxOverlap)

p7 <- ggplot(data = MaxOverlap, aes(x = Range, y =Percentage2, col = Tissue, fill = Tissue, shape = Direction)) +
  geom_rect(aes(xmin=-1.0,xmax=0,ymin=-Inf,ymax=Inf), color = NA,  fill="#f4f4f2") +
  #geom_segment(aes(x = 0, y = 0, xend = 0, yend = 100 ), color = "grey", linetype = "dashed") + # Flowers
  geom_point(size = 3 , stroke = 1) +
  scale_shape_manual(values=c(24, 25))+
  scale_color_manual(values=rep("black", 3)) +
  scale_fill_manual(values=c("#fa990e", "#155800", "#089400")) +
  scale_x_continuous(limits = c(-1, 1), breaks=seq(-1, 1.0, 0.2)) +
  scale_y_continuous(limits = c(0, 100),labels = function(x) paste0(x, "%")) +
  theme_classic() +
  theme(legend.position='none') + 
  labs(title = "Two D indices overlap",
       x = "D index ranges",
       y = "Percentage"
  ) 

ggplot(data = MaxOverlap) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 100 ), color = "grey", linetype = "dashed") + # Flowers
  geom_point(aes(x =  Range, y = Flower)) +
  geom_point(aes(x =  Range, y = Leaf)) +
  geom_point(aes(x =  Range, y = Root)) +
  geom_line(mapping = aes(x = ifelse(Range < 0, Range, NA), y = Leaf), color = "#155800", size=1, arrow = arrow(angle = 15, ends = "first", type = "closed"))+
  geom_line(mapping = aes(x = ifelse(Range < 0, Range, NA), y = Root), color = "#089400", size=1, arrow = arrow(angle = 15, ends = "first", type = "closed")) +
  geom_line(mapping = aes(x = ifelse(Range < 0, Range, NA), y = Flower), color = "#fa990e",size=1,arrow = arrow(angle = 15, ends = "first", type = "closed") )+
  geom_line(mapping = aes(x = ifelse(Range > 0, Range, NA), y = Leaf), linetype = "dashed", color = "#155800", size=1, arrow = arrow(angle = 15, ends = "last", type = "closed"))+
  geom_line(mapping = aes(x = ifelse(Range > 0, Range, NA), y = Root), linetype = "dashed", color = "#089400", size=1, arrow = arrow(angle = 15, ends = "last", type = "closed")) +
  geom_line(mapping = aes(x = ifelse(Range > 0, Range, NA), y = Flower), linetype = "dashed", color = "#fa990e",size=1,arrow = arrow(angle = 15, ends = "last", type = "closed") )+
  xlim(-1,1) + 
  ylim(0,100) + scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_classic() +
  labs(title = "Two D indices overlap (maximum)",
       x = "D index ranges",
       y = "Percentage"
  ) 

OverallOverlap <- cbind(Percentage_F, Percentage_L, Percentage_R)
OverallOverlap <- OverallOverlap[,c(1,2,6,10)]
names(OverallOverlap) <- c("Range","Flower","Leaf","Root")

p8 <- ggplot(data = OverallOverlap) +
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 100 ), color = "grey", linetype = "dashed") + # Flowers
  geom_point(aes(x =  Range, y = Flower)) +
  geom_point(aes(x =  Range, y = Leaf)) +
  geom_point(aes(x =  Range, y = Root)) +
  geom_line(mapping = aes(x = ifelse(Range < 0, Range, NA), y = Leaf), color = "#155800", size=1, arrow = arrow(angle = 15, ends = "first", type = "closed"))+
  geom_line(mapping = aes(x = ifelse(Range < 0, Range, NA), y = Root), color = "#089400", size=1, arrow = arrow(angle = 15, ends = "first", type = "closed")) +
  geom_line(mapping = aes(x = ifelse(Range < 0, Range, NA), y = Flower), color = "#fa990e",size=1,arrow = arrow(angle = 15, ends = "first", type = "closed") )+
  geom_line(mapping = aes(x = ifelse(Range > 0, Range, NA), y = Leaf), linetype = "dashed", color = "#155800", size=1, arrow = arrow(angle = 15, ends = "last", type = "closed"))+
  geom_line(mapping = aes(x = ifelse(Range > 0, Range, NA), y = Root), linetype = "dashed", color = "#089400", size=1, arrow = arrow(angle = 15, ends = "last", type = "closed")) +
  geom_line(mapping = aes(x = ifelse(Range > 0, Range, NA), y = Flower), linetype = "dashed", color = "#fa990e",size=1,arrow = arrow(angle = 15, ends = "last", type = "closed") )+
  xlim(-1,1) + 
  ylim(0,100) + scale_y_continuous(labels = function(x) paste0(x, "%")) +
  theme_classic() +
  labs(title = "Overlap (overall)",
       x = "D index ranges",
       y = "Two D indices overlap"
  ) 

ggarrange(p7, p8)

# Deviation from 0

for (tissue in c("F","L","R")) {
  if(tissue == "F"){
    dataf <- Merge_F
  }else{
    if(tissue == "L"){
      dataf <- Merge_L
    }else{
      dataf <- Merge_R
    }
  }
  
  neg<-NULL
  pos<-NULL
  neg2<-NULL
  pos2<-NULL
  k<-0
  for( i in seq(0,0.9,0.1)){
    k<-k+1
    right<-dataf$Gene[which(dataf$D < -i)]
    left<-dataf$Gene[which(dataf$NewD < -i)]
    
    right2<-dataf$Gene[which(dataf$D > i)]
    left2<-dataf$Gene[which(dataf$NewD > i)]
    
    
    neg[k]<-length(which(right %in% left)) / min(length(right),length(left))
    pos[k]<-length(which(right2 %in% left2)) / min(length(right2),length(left2))
    
    
    neg2[k]<-length(which(right %in% left)) / (length(right)+length(left)-length(which(right %in% left)))
    pos2[k]<-length(which(right2 %in% left2)) /(length(right2)+length(left2)-length(which(right2 %in% left2)))
    
  }
  MyDeviation <- rep(seq(0,0.9,0.1), 2)
  Overlap1 <- c(neg, pos) *100
  Overlap2 <- c(neg2, pos2) *100
  if(tissue == "F"){
    tissues <- rep("Flower",20)
  }else{
    if(tissue == "L"){
      tissues <- rep("Leaf",20)
    }else{
      tissues <- rep("Root",20)
    }
  }
  MyDirection <- c(rep("Negative", 10), rep("Positive", 10))
  
  MaxD <- data.frame(MyDeviation, Overlap1, MyDirection, tissues)
  names(MaxD) <- c("Deviation", "Overlap", "Direction", "Tissue")
  assign(paste0("MaxD_", tissue), MaxD)
  
  OverallD <- data.frame(MyDeviation, Overlap2, MyDirection, tissues)
  names(OverallD) <- c("Deviation", "Overlap", "Direction", "Tissue")
  assign(paste0("OverallD_", tissue), OverallD)
}

MaxD_FLR <- rbind(MaxD_F, MaxD_L, MaxD_R)
OverallD_FLR <- rbind(OverallD_F, OverallD_L, OverallD_R)

p9 <- ggplot(data = MaxD_FLR, aes(x = Deviation, y =Overlap, col = Tissue, fill = Tissue, shape = Direction)) +
  geom_point(size = 3 , stroke = 1) +
  scale_shape_manual(values=c(24, 25))+
  scale_color_manual(values=rep("black", 3)) +
  scale_fill_manual(values=c("#fa990e", "#155800", "#089400")) +
  scale_x_continuous(limits = c(0, 1), breaks=seq(0.0, 1.0, 0.2)) +
  scale_y_continuous(limits = c(0, 100),labels = function(x) paste0(x, "%")) +
  theme_bw() +
  theme(legend.position='none') + 
  labs(title = "Two D indices overlap (overall)",
       x = "Deviation from 0",
       y = "Percentage"
  ) 

p10 <- ggplot(data = OverallD_FLR, aes(x = Deviation, y =Overlap, col = Tissue, fill = Tissue, shape = Direction)) +
  geom_point(size = 3 , stroke = 1) +
  scale_shape_manual(values=c(24, 25))+
  scale_color_manual(values=rep("black", 3)) +
  scale_fill_manual(values=c("#fa990e", "#155800", "#089400")) +
  scale_x_continuous(limits = c(0, 1), breaks=seq(0.0, 1.0, 0.2)) +
  scale_y_continuous(limits = c(0, 100),labels = function(x) paste0(x, "%")) +
  theme_bw() +
  theme(legend.position='none') + 
  labs(title = "Two D indices overlap (overall)",
       x = "Deviation from 0",
       y = "Two D indices overlap"
  ) 

ggarrange(p9,p10,common.legend =T)


####
install.packages("VennDiagram")
library(VennDiagram)
grid.newpage()
draw.pairwise.venn(area1 = 638,
                   area2 = 377,
                   cross.area = 220)


draw.pairwise.venn(area1 = 2605,
                   area2 = 522,
                   cross.area = 126)


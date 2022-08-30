setwd("/Users/zebinzhang/Desktop/My_Computer/Uppsala_University/Capsella_Project/")
require(edgeR)
dataf<-read.table("InputData/Cbp_diploids_scaled_FRL_masked_RNA_ASE.csv",header = T)
head(dataf)
dim(dataf)
nrow(na.omit(dataf))
trans<-dataf$gene
dataf<-dataf[,-c(1,2:5,14:17,26:29)] # remove CR
rownames(dataf)<-trans

CG_F<-rowMeans(dataf[,c(1:4)],na.rm=T)
CO_F<-rowMeans(dataf[,c(5:8)],na.rm=T)
CG_L<-rowMeans(dataf[,c(9:12)],na.rm=T)
CO_L<-rowMeans(dataf[,c(13:16)],na.rm=T)
CG_R<-rowMeans(dataf[,c(17:20)],na.rm=T)
CO_R<-rowMeans(dataf[,c(21:24)],na.rm=T)
Cbp1_F_Co<-rowMeans(dataf[,seq(25,56,2)],na.rm=T)
Cbp1_F_Cg<-rowMeans(dataf[,seq(26,56,2)],na.rm=T)
Cbp1_L_Co<-rowMeans(dataf[,seq(57,88,2)],na.rm=T)
Cbp1_L_Cg<-rowMeans(dataf[,seq(58,88,2)],na.rm=T)
Cbp1_R_Co<-rowMeans(dataf[,seq(89,120,2)],na.rm=T)
Cbp1_R_Cg<-rowMeans(dataf[,seq(90,120,2)],na.rm=T)


dataf<-cbind.data.frame(CG_F,CG_L,CG_R,CO_F,CO_L,CO_R,Cbp1_F_Co,Cbp1_L_Co,Cbp1_R_Co,Cbp1_F_Cg,Cbp1_L_Cg,Cbp1_R_Cg)
rownames(dataf)<-trans


dataCGCbp<-read.table("Pascal_Result/DEgenes_CGvsCbp_F.csv",header=T)
dataCGCO<-read.table("Pascal_Result/DEgenes_CGvsCO_F.csv",header=T)
dataCOCbp<-read.table("Pascal_Result/DEgenes_COvsCbp_F.csv",header=T)
head(dataCGCbp)
CGCbp<-rownames(dataCGCbp)[which(dataCGCbp$FDR>=0.01)]
CGCO<-rownames(dataCGCO)[which(dataCGCO$FDR>=0.05)]
COCbp<-rownames(dataCOCbp)[which(dataCOCbp$FDR>=0.01)]
liste_F<-CGCO

dataCGCbp<-read.table("Pascal_Result/DEgenes_CGvsCbp_L.csv",header=T)
dataCGCO<-read.table("Pascal_Result/DEgenes_CGvsCO_L.csv",header=T)
dataCOCbp<-read.table("Pascal_Result/DEgenes_COvsCbp_L.csv",header=T)
head(dataCGCbp)
CGCbp<-rownames(dataCGCbp)[which(dataCGCbp$FDR>=0.01)]
CGCO<-rownames(dataCGCO)[which(dataCGCO$FDR>=0.05)]
COCbp<-rownames(dataCOCbp)[which(dataCOCbp$FDR>=0.01)]
liste_L<-CGCO

dataCGCbp<-read.table("Pascal_Result/DEgenes_CGvsCbp_R.csv",header=T)
dataCGCO<-read.table("Pascal_Result/DEgenes_CGvsCO_R.csv",header=T)
dataCOCbp<-read.table("Pascal_Result/DEgenes_COvsCbp_R.csv",header=T)
head(dataCGCbp)
CGCbp<-rownames(dataCGCbp)[which(dataCGCbp$FDR>=0.01)]
CGCO<-rownames(dataCGCO)[which(dataCGCO$FDR>=0.05)]
COCbp<-rownames(dataCOCbp)[which(dataCOCbp$FDR>=0.01)]
liste_R<-CGCO


domh1<-function(x){
  mean<-c(mean(x[which(F_d_phased_t$genome==0)],na.rm = T),mean(x[which(F_d_phased_t$genome==1)],na.rm = T),mean(x[which(F_d_phased_t$genome==2)],na.rm = T))
  mod2<-lm(mean[c(1,3)]~c(0,2))
  pred2<-mean(mod2$fitted.values)
  obs2<-mean[2]
  
  sign<-mod2$coefficients[2]
  if(sign<0){
    coeff2<-(pred2-obs2)/pred2} else {
      coeff2<-(obs2-pred2)/pred2
    }
}

domh2<-function(x){
  mean<-c(mean(x[which(F_d_phased_t$genome==0)],na.rm = T),mean(x[which(F_d_phased_t$genome==3)],na.rm = T),mean(x[which(F_d_phased_t$genome==2)],na.rm = T))
  mod2<-lm(mean[c(1,3)]~c(0,2))
  pred2<-mean(mod2$fitted.values)
  obs2<-mean[2]
  
  sign<-mod2$coefficients[2]
  if(sign<0){
    coeff2<-(pred2-obs2)/pred2} else {
      coeff2<-(obs2-pred2)/pred2
    }
}

tissue<-c("F","L","R")
ref<-c()
par(mfrow=c(3,2),oma=c(0,0,2,0))
for( i in tissue){print(i)
  if(i=="F"){
    F_d_phased<-dataf[,which(sub(".*_","",colnames(dataf))==i |
                               sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                               sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]} else {
                                 if(i=="L"){F_d_phased<-dataf[,which(sub(".*_","",colnames(dataf))==i |
                                                                       sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                       sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]}else{
                                                                         F_d_phased<-dataf[,which(sub(".*_","",colnames(dataf))==i |
                                                                                                    sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                                                    sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]
                                                                       }
                               }
  
  genome<-c(rep(0,1),rep(2,1),rep(1,2))
  F_d_phased<-na.omit(F_d_phased)
  #F_d_phased<-as.data.frame(t(F_d_phased))
  d<-DGEList(counts=F_d_phased,group=genome)
  tmmy<-calcNormFactors(d,method = "TMM")
  d<- estimateCommonDisp(tmmy)
  F_d_phased<-d$pseudo.counts
  genome2<-genome<-c(0,2,3,1)
  F_d_phased_t<-t.default(F_d_phased)
  F_d_phased_t<-cbind.data.frame(genome,F_d_phased_t)
  
  testh1<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh1(x))
  #test2<-t(test)
  #colnames(test2)<-c("coeff","coeff2")
  #hist(test,breaks=100,freq = F)
  hist(testh1[which(testh1 < 1 & testh1 > -1)],ylim=c(0,2),xlim=c(-1,1),breaks=50,freq = F, xlab = "Relative deviance index",main = paste(i,"_Cbp-Cg"," (N = ",length(testh1),")",sep=""))
  abline(v=0, col="red2",lwd=2,lty=2)
  abline(v=median(testh1,na.rm=T),col="green",lwd=2,lty=3)
  coeff_plus<-mean(testh1[which(testh1 <= 1 & testh1 > 0)])
  coeff_moins<-mean(testh1[which(testh1 <0 & testh1 >= -1)])
  ref2<-length(which(testh1<0))/length(testh1)
  text(-1,1.5,lab=paste(round(ref2,2)),pos=4)#," / ",round(coeff_moins,2),sep=""),pos=4)
  text(1,1.5,lab=paste(round(1-ref2,2)),pos=2)#," / ",round(coeff_plus,2),sep=""),pos=2)
  
  #theo<-rnorm(length(test),0,sd(test))
  #lines(density(theo),col="red",lwd=2)
  ref[length(ref)+1]<-length(which(testh1<0))/length(testh1)
  
  ################################ H2_flower ##################
  #abs(mean[3]-mean[1])
  
  testh2<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh2(x))
  #test2<-t(test)
  #colnames(test2)<-c("coeff","coeff2")
  #hist(test,breaks=100,freq = F)
  
  hist(testh2[which(testh2 < 1 & testh2 > -1)],ylim = c(0,2),xlim=c(-1,1),breaks=50,freq = F, xlab = "Relative deviance index",main = paste(i,"_Cbp-Co"," (N = ",length(testh2),")",sep=""))
  abline(v=0, col="red2",lwd=2,lty=2)
  abline(v=median(testh2,na.rm=T),col="green",lwd=2,lty=3)
  coeff_plus<-mean(testh2[which(testh2 <= 1 & testh2 > 0)])
  coeff_moins<-mean(testh2[which(testh2 <0 & testh2 >= -1)])
  ref2<-length(which(testh2<0))/length(testh2)
  text(-1,1.5,lab=paste(round(ref2,2)),pos=4)#," / ",round(coeff_moins,2),sep=""),pos=4)
  text(1,1.5,lab=paste(round(1-ref2,2)),pos=2)#," / ",round(coeff_plus,2),sep=""),pos=2)
  #theo<-rnorm(length(test),0,sd(test))
  #lines(density(theo),col="red",lwd=2)
  #length(test[which(test <= 1 & test >= -1)])/length(test)
  ref[length(ref)+1]<-length(which(testh2<0))/length(testh2)
  #length(which(test<0 & test >= -1))/length(test[which(test <= 1 & test >= -1)])
}
mtext(side = c(3),text= "All genes",font=2,outer = T)



###################### DE gene #####################
tissue<-c("F","L","R")
ref<-c()
par(mfrow=c(3,2),oma=c(0,0,2,0))
for( i in tissue){print(i)
  if(i=="F"){
    F_d_phased<-dataf[-which(rownames(dataf)%in%liste_F),which(sub(".*_","",colnames(dataf))==i |
                                                                 sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                 sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]} else {
                                                                   if(i=="L"){F_d_phased<-dataf[-which(rownames(dataf)%in%liste_L),which(sub(".*_","",colnames(dataf))==i |
                                                                                                                                           sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                                                                                           sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]}else{
                                                                                                                                             F_d_phased<-dataf[-which(rownames(dataf)%in%liste_R),which(sub(".*_","",colnames(dataf))==i |
                                                                                                                                                                                                          sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                                                                                                                                                          sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]
                                                                                                                                           }
                                                                 }
  
  genome<-c(rep(0,1),rep(2,1),rep(1,2))
  F_d_phased<-na.omit(F_d_phased)
  #F_d_phased<-as.data.frame(t(F_d_phased))
  d<-DGEList(counts=F_d_phased,group=genome)
  tmmy<-calcNormFactors(d,method = "TMM")
  d<- estimateCommonDisp(tmmy)
  F_d_phased<-d$pseudo.counts
  genome2<-genome<-c(0,2,3,1)
  F_d_phased_t<-t.default(F_d_phased)
  F_d_phased_t<-cbind.data.frame(genome,F_d_phased_t)
  
  testh1<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh1(x))
  #test2<-t(test)
  #colnames(test2)<-c("coeff","coeff2")
  #hist(test,breaks=100,freq = F)
  hist(testh1[which(testh1 < 1 & testh1 > -1)],ylim=c(0,2),xlim=c(-1,1),breaks=50,freq = F, xlab = "Relative deviance index",main = paste(i,"_Cbp-Cg"," (N = ",length(testh1),")",sep=""))
  abline(v=0, col="red2",lwd=2,lty=2)
  abline(v=median(testh1,na.rm=T),col="green",lwd=2,lty=3)
  coeff_plus<-mean(testh1[which(testh1 <= 1 & testh1 > 0)])
  coeff_moins<-mean(testh1[which(testh1 <0 & testh1 >= -1)])
  ref2<-length(which(testh1<0))/length(testh1)
  text(-1,1.5,lab=paste(round(ref2,2)),pos=4)#," / ",round(coeff_moins,2),sep=""),pos=4)
  text(1,1.5,lab=paste(round(1-ref2,2)),pos=2)#," / ",round(coeff_plus,2),sep=""),pos=2)
  
  #theo<-rnorm(length(test),0,sd(test))
  #lines(density(theo),col="red",lwd=2)
  ref[length(ref)+1]<-length(which(testh1<0))/length(testh1)
  
  ################################ H2_flower ##################
  #abs(mean[3]-mean[1])
  
  testh2<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh2(x))
  #test2<-t(test)
  #colnames(test2)<-c("coeff","coeff2")
  #hist(test,breaks=100,freq = F)
  
  hist(testh2[which(testh2 < 1 & testh2 > -1)],ylim = c(0,2),xlim=c(-1,1),breaks=50,freq = F, xlab = "Relative deviance index",main = paste(i,"_Cbp-Co"," (N = ",length(testh2),")",sep=""))
  abline(v=0, col="red2",lwd=2,lty=2)
  abline(v=median(testh2,na.rm=T),col="green",lwd=2,lty=3)
  coeff_plus<-mean(testh2[which(testh2 <= 1 & testh2 > 0)])
  coeff_moins<-mean(testh2[which(testh2 <0 & testh2 >= -1)])
  ref2<-length(which(testh2<0))/length(testh2)
  text(-1,1.5,lab=paste(round(ref2,2)),pos=4)#," / ",round(coeff_moins,2),sep=""),pos=4)
  text(1,1.5,lab=paste(round(1-ref2,2)),pos=2)#," / ",round(coeff_plus,2),sep=""),pos=2)
  #theo<-rnorm(length(test),0,sd(test))
  #lines(density(theo),col="red",lwd=2)
  #length(test[which(test <= 1 & test >= -1)])/length(test)
  ref[length(ref)+1]<-length(which(testh2<0))/length(testh2)
  #length(which(test<0 & test >= -1))/length(test[which(test <= 1 & test >= -1)])
}
mtext(side = c(3),text= "Differentially expressed genes",font=2,outer = T)




######################################################################################################
### FILTRAGE BY POP

d <- read.table("~/Documents/Post-doc/Marion_distib/Cbp_diploids_scaled_FRL_masked_RNA_ASE.csv",header=TRUE, row.names="gene")
annot <- read.csv('~/Documents/Post-doc/papier4/first_results/Dima/28genomes_annot_R.csv', header = T)

# Define populations
CR_F <- c(annot$populations=="CR" & annot$tissue=="flower")
CG_F <- c(annot$populations=="CG" & annot$tissue=="flower")
CO_F <- c(annot$populations=="CO" & annot$tissue=="flower")

CR_L <- c(annot$populations=="CR" & annot$tissue=="leaf")
CG_L <- c(annot$populations=="CG" & annot$tissue=="leaf")
CO_L <- c(annot$populations=="CO" & annot$tissue=="leaf")

CR_R <- c(annot$populations=="CR" & annot$tissue=="root")
CG_R <- c(annot$populations=="CG" & annot$tissue=="root")
CO_R <- c(annot$populations=="CO" & annot$tissue=="root")

Cbp_Co <- c("ASI_Co", "EUR_Co", "ME_Co", "CASI_Co")
Cbp_Cg <- c("ASI_Cg", "EUR_Cg", "ME_Cg", "CASI_Cg")

Cbp_Co_F <- c(annot$populations %in% Cbp_Co & annot$tissue=="flower")
Cbp_Co_L <- c(annot$populations %in% Cbp_Co & annot$tissue=="leaf")
Cbp_Co_R <- c(annot$populations %in% Cbp_Co & annot$tissue=="root")
Cbp_Cg_F <- c(annot$populations %in% Cbp_Cg & annot$tissue=="flower")
Cbp_Cg_L <- c(annot$populations %in% Cbp_Cg & annot$tissue=="leaf")
Cbp_Cg_R <- c(annot$populations %in% Cbp_Cg & annot$tissue=="root")

ASI_Co_F <- c(annot$populations=="ASI_Co" & annot$tissue=="flower")
EUR_Co_F <- c(annot$populations=="EUR_Co" & annot$tissue=="flower")
ME_Co_F <- c(annot$populations=="ME_Co" & annot$tissue=="flower")
CASI_Co_F <- c(annot$populations=="CASI_Co" & annot$tissue=="flower")

ASI_Cg_F <- c(annot$populations=="ASI_Cg" & annot$tissue=="flower")
EUR_Cg_F <- c(annot$populations=="EUR_Cg" & annot$tissue=="flower")
ME_Cg_F <- c(annot$populations=="ME_Cg" & annot$tissue=="flower")
CASI_Cg_F <- c(annot$populations=="CASI_Cg" & annot$tissue=="flower")

ASI_Co_L <- c(annot$populations=="ASI_Co" & annot$tissue=="leaf")
EUR_Co_L <- c(annot$populations=="EUR_Co" & annot$tissue=="leaf")
ME_Co_L <- c(annot$populations=="ME_Co" & annot$tissue=="leaf")
CASI_Co_L <- c(annot$populations=="CASI_Co" & annot$tissue=="leaf")

ASI_Cg_L <- c(annot$populations=="ASI_Cg" & annot$tissue=="leaf")
EUR_Cg_L <- c(annot$populations=="EUR_Cg" & annot$tissue=="leaf")
ME_Cg_L <- c(annot$populations=="ME_Cg" & annot$tissue=="leaf")
CASI_Cg_L <- c(annot$populations=="CASI_Cg" & annot$tissue=="leaf")

ASI_Co_R <- c(annot$populations=="ASI_Co" & annot$tissue=="root")
EUR_Co_R <- c(annot$populations=="EUR_Co" & annot$tissue=="root")
ME_Co_R <- c(annot$populations=="ME_Co" & annot$tissue=="root")
CASI_Co_R <- c(annot$populations=="CASI_Co" & annot$tissue=="root")

ASI_Cg_R <- c(annot$populations=="ASI_Cg" & annot$tissue=="root")
EUR_Cg_R <- c(annot$populations=="EUR_Cg" & annot$tissue=="root")
ME_Cg_R <- c(annot$populations=="ME_Cg" & annot$tissue=="root")
CASI_Cg_R <- c(annot$populations=="CASI_Cg" & annot$tissue=="root")


#Define the filtering criteria and filter
popMis <- 1
popKeep <- cbind( (rowSums(!is.na(d[,ASI_Co_F]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,EUR_Co_F]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,ME_Co_F]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,CASI_Co_F]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,ASI_Cg_F]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,EUR_Cg_F]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,ME_Cg_F]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,CASI_Cg_F]) >= 1) >= popMis) &
                    
                    (rowSums(!is.na(d[,ASI_Co_L]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,EUR_Co_L]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,ME_Co_L]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,CASI_Co_L]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,ASI_Cg_L]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,EUR_Cg_L]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,ME_Cg_L]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,CASI_Cg_L]) >= 1) >= popMis) &
                    
                    (rowSums(!is.na(d[,ASI_Co_R]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,EUR_Co_R]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,ME_Co_R]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,CASI_Co_R]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,ASI_Cg_R]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,EUR_Cg_R]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,ME_Cg_R]) >= 1) >= popMis) &
                    (rowSums(!is.na(d[,CASI_Cg_R]) >= 1) >= popMis))
dfpop <- d[popKeep,]
#message(paste("Number of genes before", dim(d)[1], "\nNumber of genes after", dim(dfpop)[1]))


trans<-rownames(dfpop)
dataf<-dfpop[,-c(1:4,13:16,25:28)]
rownames(dataf)<-trans

CG_F<-rowMeans(dataf[,c(1:4)],na.rm=T)
CO_F<-rowMeans(dataf[,c(5:8)],na.rm=T)
CG_L<-rowMeans(dataf[,c(9:12)],na.rm=T)
CO_L<-rowMeans(dataf[,c(13:16)],na.rm=T)
CG_R<-rowMeans(dataf[,c(17:20)],na.rm=T)
CO_R<-rowMeans(dataf[,c(21:24)],na.rm=T)
Cbp1_F_Co<-rowMeans(dataf[,seq(25,56,2)],na.rm=T)
Cbp1_F_Cg<-rowMeans(dataf[,seq(26,56,2)],na.rm=T)
Cbp1_L_Co<-rowMeans(dataf[,seq(57,88,2)],na.rm=T)
Cbp1_L_Cg<-rowMeans(dataf[,seq(58,88,2)],na.rm=T)
Cbp1_R_Co<-rowMeans(dataf[,seq(89,120,2)],na.rm=T)
Cbp1_R_Cg<-rowMeans(dataf[,seq(90,120,2)],na.rm=T)


dataf<-cbind.data.frame(CG_F,CG_L,CG_R,CO_F,CO_L,CO_R,Cbp1_F_Co,Cbp1_L_Co,Cbp1_R_Co,Cbp1_F_Cg,Cbp1_L_Cg,Cbp1_R_Cg)
rownames(dataf)<-trans

domh1<-function(x){
  mean<-c(mean(x[which(F_d_phased_t$genome==0)],na.rm = T),mean(x[which(F_d_phased_t$genome==1)],na.rm = T),mean(x[which(F_d_phased_t$genome==2)],na.rm = T))
  mod2<-lm(mean[c(1,3)]~c(0,2))
  pred2<-mean(mod2$fitted.values)
  obs2<-mean[2]
  
  sign<-mod2$coefficients[2]
  if(sign<0){
    coeff2<-(pred2-obs2)/pred2} else {
      coeff2<-(obs2-pred2)/pred2
    }
}

domh2<-function(x){
  mean<-c(mean(x[which(F_d_phased_t$genome==0)],na.rm = T),mean(x[which(F_d_phased_t$genome==3)],na.rm = T),mean(x[which(F_d_phased_t$genome==2)],na.rm = T))
  mod2<-lm(mean[c(1,3)]~c(0,2))
  pred2<-mean(mod2$fitted.values)
  obs2<-mean[2]
  
  sign<-mod2$coefficients[2]
  if(sign<0){
    coeff2<-(pred2-obs2)/pred2} else {
      coeff2<-(obs2-pred2)/pred2
    }
}

tissue<-c("F","L","R")
ref<-c()
cor2<-c()
par(mfrow=c(3,3),oma=c(0,0,2,0),mar=c(4,4,4,1))
for( i in tissue){print(i)
  if(i=="F"){
    F_d_phased<-dataf[,which(sub(".*_","",colnames(dataf))==i |
                               sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                               sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]} else {
                                 if(i=="L"){F_d_phased<-dataf[,which(sub(".*_","",colnames(dataf))==i |
                                                                       sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                       sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]}else{
                                                                         F_d_phased<-dataf[,which(sub(".*_","",colnames(dataf))==i |
                                                                                                    sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                                                    sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]
                                                                       }
                               }
  
  genome<-c(rep(0,1),rep(2,1),rep(1,2))
  F_d_phased<-na.omit(F_d_phased)
  #F_d_phased<-as.data.frame(t(F_d_phased))
  d<-DGEList(counts=F_d_phased,group=genome)
  tmmy<-calcNormFactors(d,method = "TMM")
  d<- estimateCommonDisp(tmmy)
  relation<-as.factor(rownames(d$samples))
  design<-model.matrix(~relation,datat=d)
  
  dGLM<-estimateGLMCommonDisp(d,design)
  
  
  F_d_phased<-d$pseudo.counts
  genome2<-genome<-c(0,2,3,1)
  F_d_phased_t<-t.default(F_d_phased)
  F_d_phased_t<-cbind.data.frame(genome,F_d_phased_t)
  if(i=="F"){titre<-"Flowers"}else{if(i=="L"){titre<-"Leaves"}else{titre<-"Roots"}}
  
  testh1<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh1(x))
  #test2<-t(test)
  #colnames(test2)<-c("coeff","coeff2")
  #hist(test,breaks=100,freq = F)
  hist(testh1[which(testh1 < 1 & testh1 > -1)],ylim=c(0,2.5),xlim=c(-1,1),breaks=50,freq = F, xlab = expression(italic(S[Cg])),main = paste(titre," (",length(which(testh1!=0)),")",sep=""))
  abline(v=0, col="red2",lwd=2,lty=2)
  abline(v=median(testh1,na.rm=T),col="green",lwd=2,lty=3)
  coeff_plus<-mean(testh1[which(testh1 <= 1 & testh1 > 0)])
  coeff_moins<-mean(testh1[which(testh1 <0 & testh1 >= -1)])
  ref2<-length(which(testh1<0))/length(testh1)
  b_test<-binom.test(c(length(which(testh1<0)),length(which(testh1>0))),p=0.5)
  if(b_test$p.value>=0.05){out<-paste(round(ref2,2),"ns",sep=" ")}else{if(b_test$p.value<0.001){out<-paste(round(ref2,2),"***",sep=" ")}else{if(b_test$p.value<0.01){out<-paste(round(ref2,2),"**",sep=" ")}else{out<-paste(round(ref2,2),"*",sep=" ")}}}
  text(-1,1.5,lab=out,pos=4)#," / ",round(coeff_moins,2),sep=""),pos=4)
  text(1,1.5,lab=paste(round(1-ref2,2)),pos=2)#," / ",round(coeff_plus,2),sep=""),pos=2)
  #theo<-rnorm(length(test),0,sd(test))
  #lines(density(theo),col="red",lwd=2)
  ref[length(ref)+1]<-median(testh1,na.rm=T)
  
  ################################ H2_flower ##################
  #abs(mean[3]-mean[1])
  
  testh2<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh2(x))
  #test2<-t(test)
  #colnames(test2)<-c("coeff","coeff2")
  #hist(test,breaks=100,freq = F)
  
  hist(testh2[which(testh2 < 1 & testh2 > -1)],ylim = c(0,2.5),xlim=c(-1,1),breaks=50,freq = F, xlab = expression(italic(S[Co])),main = paste(titre," (",length(which(testh2!=0)),")",sep=""))
  abline(v=0, col="red2",lwd=2,lty=2)
  abline(v=median(testh2,na.rm=T),col="green",lwd=2,lty=3)
  coeff_plus<-mean(testh2[which(testh2 <= 1 & testh2 > 0)])
  coeff_moins<-mean(testh2[which(testh2 <0 & testh2 >= -1)])
  ref2<-length(which(testh2<0))/length(testh2)
  b_test<-binom.test(c(length(which(testh2<0)),length(which(testh2>0))),p=0.5)
  if(b_test$p.value>=0.05){out<-paste(round(ref2,2),"ns",sep=" ")}else{if(b_test$p.value<0.001){out<-paste(round(ref2,2),"***",sep=" ")}else{if(b_test$p.value<0.01){out<-paste(round(ref2,2),"**",sep=" ")}else{out<-paste(round(ref2,2),"*",sep=" ")}}}
  text(-1,1.5,lab=out,pos=4)#," / ",round(coeff_moins,2),sep=""),pos=4)
  text(1,1.5,lab=paste(round(1-ref2,2)),pos=2)#," / ",round(coeff_plus,2),sep=""),pos=2)
  #theo<-rnorm(length(test),0,sd(test))
  #lines(density(theo),col="red",lwd=2)
  #length(test[which(test <= 1 & test >= -1)])/length(test)
  ref[length(ref)+1]<-median(testh2,na.rm=T)
  cuth1<-names(head(sort(abs(testh1),decreasing  =T),50))
  cuth2<-names(head(sort(abs(testh2),decreasing  =T),50))
  cut<-unique(c(cuth1,cuth2))
  h1<-testh1[-which(names(testh1)%in%cut)]
  h2<-testh2[-which(names(testh2)%in%cut)]
  h1<-h1[order(names(h1))]
  h2<-h2[order(names(h2))]
  cor<-cor.test(h1,h2,method="spearman")
  cor<-c(cor$estimate,cor$parameter,cor$statistic,cor$p.value)
  cor2<-rbind.data.frame(cor2,cor)
  
  plot(h2~h1,xlab=expression(italic(S[Cg])),asp=1,ylab=expression(italic(S[Co])),main=paste(titre," (",length(h1),")",sep=""))
  abline(lm(h2~h1),lty=2,col="red")
  #length(which(test<0 & test >= -1))/length(test[which(test <= 1 & test >= -1)])
  #print(wilcox.test(testh1,testh2,paired = F))
  
}
mtext(side = c(3),text= "All genes",font=2,outer = T)
cor2


###################### DE gene #####################
tissue<-c("F","L","R")
ref<-c()
cor2<-c()
new2<-c()
par(mfrow=c(3,3),oma=c(0,0,2,0),mar=c(4,4,4,1))
for( i in tissue){print(i)
  if(i=="F"){
    F_d_phased<-dataf[-which(rownames(dataf)%in%liste_F),which(sub(".*_","",colnames(dataf))==i |
                                                                 sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                 sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]} else {
                                                                   if(i=="L"){F_d_phased<-dataf[-which(rownames(dataf)%in%liste_L),which(sub(".*_","",colnames(dataf))==i |
                                                                                                                                           sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                                                                                           sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]}else{
                                                                                                                                             F_d_phased<-dataf[-which(rownames(dataf)%in%liste_R),which(sub(".*_","",colnames(dataf))==i |
                                                                                                                                                                                                          sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                                                                                                                                                          sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]
                                                                                                                                           }
                                                                 }
  
  genome<-c(rep(0,1),rep(2,1),rep(1,2))
  F_d_phased<-na.omit(F_d_phased)
  #F_d_phased<-as.data.frame(t(F_d_phased))
  d<-DGEList(counts=F_d_phased,group=genome)
  tmmy<-calcNormFactors(d,method = "TMM")
  d<- estimateCommonDisp(tmmy)
  F_d_phased<-d$pseudo.counts
  genome2<-genome<-c(0,2,3,1)
  F_d_phased_t<-t.default(F_d_phased)
  F_d_phased_t<-cbind.data.frame(genome,F_d_phased_t)
  if(i=="F"){titre<-"Flowers"}else{if(i=="L"){titre<-"Leaves"}else{titre<-"Roots"}}
  
  testh1<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh1(x))
  #test2<-t(test)
  #colnames(test2)<-c("coeff","coeff2")
  #hist(test,breaks=100,freq = F)
  hist(testh1[which(testh1 < 1 & testh1 > -1)],ylim=c(0,2.5),xlim=c(-1,1),breaks=50,freq = F, xlab = expression(italic(S[Cg])),main = paste(titre," (",length(which(testh1!=0)),")",sep=""))
  abline(v=0, col="red2",lwd=2,lty=2)
  abline(v=median(testh1,na.rm=T),col="green",lwd=2,lty=3)
  coeff_plus<-mean(testh1[which(testh1 <= 1 & testh1 > 0)])
  coeff_moins<-mean(testh1[which(testh1 <0 & testh1 >= -1)])
  ref2<-length(which(testh1<0))/length(testh1)
  b_test<-binom.test(c(length(which(testh1<0)),length(which(testh1>0))),p=0.5)
  if(b_test$p.value>=0.05){out<-paste(round(ref2,2),"ns",sep=" ")}else{if(b_test$p.value<0.001){out<-paste(round(ref2,2),"***",sep=" ")}else{if(b_test$p.value<0.01){out<-paste(round(ref2,2),"**",sep=" ")}else{out<-paste(round(ref2,2),"*",sep=" ")}}}
  text(-1,1.5,lab=out,pos=4)#," / ",round(coeff_moins,2),sep=""),pos=4)
  text(1,1.5,lab=paste(round(1-ref2,2)),pos=2)#," / ",round(coeff_plus,2),sep=""),pos=2)
  
  #theo<-rnorm(length(test),0,sd(test))
  #lines(density(theo),col="red",lwd=2)
  ref[length(ref)+1]<-median(testh1,na.rm=T)
  
  ################################ H2_flower ##################
  #abs(mean[3]-mean[1])
  
  testh2<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh2(x))
  #test2<-t(test)
  #colnames(test2)<-c("coeff","coeff2")
  #hist(test,breaks=100,freq = F)
  
  hist(testh2[which(testh2 < 1 & testh2 > -1)],ylim = c(0,2.5),xlim=c(-1,1),breaks=50,freq = F, xlab = expression(italic(S[Co])),main = paste(titre," (",length(which(testh1!=0)),")",sep=""))
  abline(v=0, col="red2",lwd=2,lty=2)
  abline(v=median(testh2,na.rm=T),col="green",lwd=2,lty=3)
  coeff_plus<-mean(testh2[which(testh2 <= 1 & testh2 > 0)])
  coeff_moins<-mean(testh2[which(testh2 <0 & testh2 >= -1)])
  ref2<-length(which(testh2<0))/length(testh2)
  b_test<-binom.test(c(length(which(testh2<0)),length(which(testh2>0))),p=0.5)
  if(b_test$p.value>=0.05){out<-paste(round(ref2,2),"ns",sep=" ")}else{if(b_test$p.value<0.001){out<-paste(round(ref2,2),"***",sep=" ")}else{if(b_test$p.value<0.01){out<-paste(round(ref2,2),"**",sep=" ")}else{out<-paste(round(ref2,2),"*",sep=" ")}}}
  text(-1,1.5,lab=out,pos=4)#," / ",round(coeff_moins,2),sep=""),pos=4)
  text(1,1.5,lab=paste(round(1-ref2,2)),pos=2)#," / ",round(coeff_plus,2),sep=""),pos=2)
  #theo<-rnorm(length(test),0,sd(test))
  #lines(density(theo),col="red",lwd=2)
  #length(test[which(test <= 1 & test >= -1)])/length(test)
  ref[length(ref)+1]<-median(testh2,na.rm=T)
  #print(wilcox.test(testh1,testh2,paired = F))
  cuth1<-names(head(sort(abs(testh1),decreasing  =T),50))
  cuth2<-names(head(sort(abs(testh2),decreasing  =T),50))
  cut<-unique(c(cuth1,cuth2))
  h1<-testh1[-which(names(testh1)%in%cut)]
  h2<-testh2[-which(names(testh2)%in%cut)]
  h1<-h1[order(names(h1))]
  h2<-h2[order(names(h2))]
  cor<-cor.test(h1,h2,method="spearman")
  cor<-c(cor$estimate,cor$parameter,cor$statistic,cor$p.value)
  cor2<-rbind.data.frame(cor2,cor)
  new<-c(median(testh1),median(testh2),(abs(median(testh2))-abs(median(testh1))))
  new2<-rbind(new2,new)
  plot(h2~h1,xlab=expression(italic(S[Cg])),asp=1,ylab=expression(italic(S[Co])),main=paste(titre," (",length(h1),")",sep=""))
  abline(lm(h2~h1),lty=2,col="red")
  #r_test<-cor[4]
  #if(r_test>=0.05){out<-paste(round(cor[1],2),"ns",sep=" ")}else{if(r_test<0.001){out<-paste(round(cor[1],2),"***",sep=" ")}else{if(r_test<0.01){out<-paste(round(cor[1],2),"**",sep=" ")}else{out<-paste(round(cor[1],2),"*",sep=" ")}}}
  #legend("topleft",paste("r = ",out,sep=""),bty = "n",y.intersp = 0)
  #text(-3.5,1.7,paste("r = ",out,sep=""),pos=4)
  #length(which(test<0 & test >= -1))/length(test[which(test <= 1 & test >= -1)])
}
mtext(side = c(3),text= "Differentially expressed genes",font=2,outer = T)

cor2
new2
#################### ratios #########

tissue<-c("F","L","R")
par(mfrow=c(2,3),oma=c(1,0,0,0,0))
for( i in tissue){print(i)
  if(i=="F"){
    F_d_phased<-dataf[,which(sub(".*_","",colnames(dataf))==i |
                               sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                               sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]} else {
                                 if(i=="L"){F_d_phased<-dataf[,which(sub(".*_","",colnames(dataf))==i |
                                                                       sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                       sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]}else{
                                                                         F_d_phased<-dataf[,which(sub(".*_","",colnames(dataf))==i |
                                                                                                    sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                                                    sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]
                                                                       }
                               }
  
  if(i=="F"){titre<-"Flowers"}else{if(i=="L"){titre<-"Leaves"}else{titre<-"Roots"}}
  
  subgenomes<-F_d_phased[,4]/(F_d_phased[,4]+F_d_phased[,3])
  genomes<-F_d_phased[,1]/(F_d_phased[,1]+F_d_phased[,2])
  plot(subgenomes~genomes,ylab=expression(Cbp[Cg]/(Cbp[Cg]+Cbp[Co])),ylim=c(0,1.2),xlab=expression(CG/(CG+CO)),main=paste(titre," (",length(na.omit(genomes)),")",sep=""))
  abline(lm(subgenomes~genomes),col="red",lty=2)
  abline(v=0.5,col="gray40",lty=3)
  abline(h=0.5,col="gray40",lty=3)
  beta<-round(summary(lm(subgenomes~genomes))$coefficients[2,1],2)
  text(0.0,1.2,expression(italic(beta)),col="red",pos=4,lwd=1.3)
  text(0.02,1.2,paste(" = ",beta,sep=""),col="red",pos=4,lwd=1.3)
}

for( i in tissue){print(i) 
  if(i=="F"){titre<-"Flowers"}else{if(i=="L"){titre<-"Leaves"}else{titre<-"Roots"}}
  if(i=="F"){
    F_d_phased<-dataf[-which(rownames(dataf)%in%liste_F),which(sub(".*_","",colnames(dataf))==i |
                                                                 sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                 sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]} else {
                                                                   if(i=="L"){F_d_phased<-dataf[-which(rownames(dataf)%in%liste_L),which(sub(".*_","",colnames(dataf))==i |
                                                                                                                                           sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                                                                                           sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]}else{
                                                                                                                                             F_d_phased<-dataf[-which(rownames(dataf)%in%liste_R),which(sub(".*_","",colnames(dataf))==i |
                                                                                                                                                                                                          sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Co",sep="") |
                                                                                                                                                                                                          sub(".*[0-9]_","",colnames(dataf))==paste(i,"_Cg",sep=""))]
                                                                                                                                           }
                                                                 }
  subgenomes<-F_d_phased[,4]/(F_d_phased[,4]+F_d_phased[,3])
  genomes<-F_d_phased[,1]/(F_d_phased[,1]+F_d_phased[,2])
  plot(subgenomes~genomes,ylab=expression(Cbp[Cg]/(Cbp[Cg]+Cbp[Co])),ylim=c(0,1.2),xlab=expression(CG/(CG+CO)),main=paste(titre," (",length(na.omit(genomes)),")",sep=""))
  abline(lm(subgenomes~genomes),col="red",lty=2)
  abline(v=0.5,col="gray40",lty=3)
  abline(h=0.5,col="gray40",lty=3)
  beta<-round(summary(lm(subgenomes~genomes))$coefficients[2,1],2)
  text(0.0,1.2,expression(italic(beta)),col="red",pos=4,lwd=1.3)
  text(0.02,1.2,paste(" = ",beta,sep=""),col="red",pos=4,lwd=1.3)
  
}




subgenomes<-F_d_phased$Cbp1_F_Cg/(F_d_phased$Cbp1_F_Cg+F_d_phased$Cbp1_F_Co)
genomes<-F_d_phased$CG_F/(F_d_phased$CG_F+F_d_phased$CO_F)



subgenomes_fc<-log2(F_d_phased$Cbp1_F_Cg+1)-log2(F_d_phased$Cbp1_F_Co+1)
genomes_fc<-log2(F_d_phased$CG_F+1)-log2(F_d_phased$CO_F+1)                 

par(mfrow=c(1,2))
plot(subgenomes~genomes,ylab=expression(Cbp[Cg]/(Cbp[Cg]+Cbp[Co])),xlab=expression(CG/(CG+CO)),main=titre)
abline(lm(subgenomes~genomes),col="red",lty=2)
abline(v=0.5,col="gray40",lty=3)
abline(h=0.5,col="gray40",lty=3)
summary(lm(subgenomes~genomes))

plot(subgenomes_fc~genomes_fc,ylab=expression(Cbp[Cg]/(Cbp[Cg]+Cbp[Co])),xlab=expression(CG/(CG+CO)),main=titre)
abline(lm(subgenomes_fc~genomes_fc),col="red",lty=2)
abline(v=0.5,col="green")
abline(h=0.5,col="green")




subgdelta<-F_d_phased[,4]-F_d_phased[,3]
gdelta<-F_d_phased[,1]-F_d_phased[,2]


abline(h=c(-1000,1000),col="green")
abline(v=c(-1000,1000),col="green")
plot(subgdelta~gdelta,xlim=c(-20000,20000),ylim=c(-20000,20000))
points(seq(-20000,20000,100),seq(-20000,20000,100),type="l",lty=3,col="gray40")
abline(lm(subgdelta~gdelta),col="red",lty=2)

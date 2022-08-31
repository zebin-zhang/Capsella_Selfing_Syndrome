setwd("~/Documents/Post-doc/papier5/")

dataCGCbp<-read.table("~/Documents/Post-doc/Marion_distib/DEgenes_CGvsCbp_F.csv",header=T)
dataCGCO<-read.table("~/Documents/Post-doc/Marion_distib/DEgenes_CGvsCO_F.csv",header=T)
dataCOCbp<-read.table("~/Documents/Post-doc/Marion_distib/DEgenes_COvsCbp_F.csv",header=T)
head(dataCGCbp)
CGCbp<-rownames(dataCGCbp)[which(dataCGCbp$FDR>=0.01)]
CGCO<-rownames(dataCGCO)[which(dataCGCO$FDR>=0.05)]
COCbp<-rownames(dataCOCbp)[which(dataCOCbp$FDR>=0.01)]
liste_F<-CGCO

dataCGCbp<-read.table("~/Documents/Post-doc/Marion_distib/DEgenes_CGvsCbp_L.csv",header=T)
dataCGCO<-read.table("~/Documents/Post-doc/Marion_distib/DEgenes_CGvsCO_L.csv",header=T)
dataCOCbp<-read.table("~/Documents/Post-doc/Marion_distib/DEgenes_COvsCbp_L.csv",header=T)
head(dataCGCbp)
CGCbp<-rownames(dataCGCbp)[which(dataCGCbp$FDR>=0.01)]
CGCO<-rownames(dataCGCO)[which(dataCGCO$FDR>=0.05)]
COCbp<-rownames(dataCOCbp)[which(dataCOCbp$FDR>=0.01)]
liste_L<-CGCO

dataCGCbp<-read.table("~/Documents/Post-doc/Marion_distib/DEgenes_CGvsCbp_R.csv",header=T)
dataCGCO<-read.table("~/Documents/Post-doc/Marion_distib/DEgenes_CGvsCO_R.csv",header=T)
dataCOCbp<-read.table("~/Documents/Post-doc/Marion_distib/DEgenes_COvsCbp_R.csv",header=T)
head(dataCGCbp)
CGCbp<-rownames(dataCGCbp)[which(dataCGCbp$FDR>=0.01)]
CGCO<-rownames(dataCGCO)[which(dataCGCO$FDR>=0.05)]
COCbp<-rownames(dataCOCbp)[which(dataCOCbp$FDR>=0.01)]
liste_R<-CGCO


dataf<-read.table("Cbp_diploids_scaled_FRL_masked_RNA_ASE.csv",header = T)
colnames(dataf)
nrow(na.omit(dataf))
#trans<-dataf$gene#rownames(dataf)
d <- read.table("Cbp_diploids_scaled_FRL_masked_RNA_ASE.csv",header=TRUE, row.names="gene")
annot <- read.csv("28genomes_annot_R.csv", header = T,sep=";")

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
dataf<-dfpop
label<-paste(annot$accession,annot$populations,sep=":")
colnames(dataf)<-label
rownames(dataf)<-trans

CR_F<-rowMeans(dataf[,which(sub(".*_","",colnames(dataf))=="F:CR")],na.rm=T)
CG_F<-rowMeans(dataf[,which(sub(".*_","",colnames(dataf))=="F:CG")],na.rm=T)
CO_F<-rowMeans(dataf[,which(sub(".*_","",colnames(dataf))=="F:CO")],na.rm=T)
Cbp1_F_Co<-rowMeans(dataf[,which(sub("[^_]*","",sub(":.*","",colnames(dataf)))=="_F_Co")],na.rm=T)
Cbp1_F_Cg<-rowMeans(dataf[,which(sub("[^_]*","",sub(":.*","",colnames(dataf)))=="_F_Cg")],na.rm=T)

CR_L<-rowMeans(dataf[,which(sub(".*_","",colnames(dataf))=="L:CR")],na.rm=T)
CG_L<-rowMeans(dataf[,which(sub(".*_","",colnames(dataf))=="L:CG")],na.rm=T)
CO_L<-rowMeans(dataf[,which(sub(".*_","",colnames(dataf))=="L:CO")],na.rm=T)
Cbp1_L_Co<-rowMeans(dataf[,which(sub("[^_]*","",sub(":.*","",colnames(dataf)))=="_L_Co")],na.rm=T)
Cbp1_L_Cg<-rowMeans(dataf[,which(sub("[^_]*","",sub(":.*","",colnames(dataf)))=="_L_Cg")],na.rm=T)

CR_R<-rowMeans(dataf[,which(sub(".*_","",colnames(dataf))=="R:CR")],na.rm=T)
CG_R<-rowMeans(dataf[,which(sub(".*_","",colnames(dataf))=="R:CG")],na.rm=T)
CO_R<-rowMeans(dataf[,which(sub(".*_","",colnames(dataf))=="R:CO")],na.rm=T)
Cbp1_R_Co<-rowMeans(dataf[,which(sub("[^_]*","",sub(":.*","",colnames(dataf)))=="_R_Co")],na.rm=T)
Cbp1_R_Cg<-rowMeans(dataf[,which(sub("[^_]*","",sub(":.*","",colnames(dataf)))=="_R_Cg")],na.rm=T)


dataf<-cbind.data.frame(CR_F,CR_L,CR_R,CG_F,CG_L,CG_R,CO_F,CO_L,CO_R,Cbp1_F_Co,Cbp1_L_Co,Cbp1_R_Co,Cbp1_F_Cg,Cbp1_L_Cg,Cbp1_R_Cg)
rownames(dataf)<-trans

################################# Ratio ###########################
color<- c("#D55E00","#CC79A7","#E69F00","#F0E442", "#56B4E9","#0072B2","#009E73","#999999")
tissue<-c("F","L","R")
par(mfrow=c(1,3),xpd=F)
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
  genome<-c(rep(0,1),rep(2,1),rep(3,1),rep(1,2))
  F_d_phased<-na.omit(F_d_phased)
  #F_d_phased<-as.data.frame(t(F_d_phased))
  d<-DGEList(counts=F_d_phased,group=genome)
  tmmy<-calcNormFactors(d,method = "TMM")
  d<- estimateCommonDisp(tmmy)
  F_d_phased<-d$pseudo.counts
  if(length(which(F_d_phased==0))>0){
    F_d_phased<-F_d_phased[-which(F_d_phased[,1]==0 |F_d_phased[,2]==0 |F_d_phased[,3]==0 |F_d_phased[,4]==0 |F_d_phased[,5]==0),]
  }
  genome2<-genome<-c(0,2,3,1)
  F_d_phased_t<-t.default(F_d_phased)
  #which(F_d_phased==0)
  
  if(i=="F"){titre<-"Flowers"}else{if(i=="L"){titre<-"Leaves"}else{titre<-"Roots"}}
  
  subgenomes<-F_d_phased[,1]/(F_d_phased[,1]+F_d_phased[,3])
  var(subgenomes)
  genomes<-F_d_phased[,2]/(F_d_phased[,2]+F_d_phased[,3])
  var(genomes)
  print(var.test(subgenomes,genomes))
  #write.table(cbind(subgenomes),paste("~/Documents/Post-doc/Marion_distib/ratio_sub_",i,".txt",sep=""),col.names = F,row.names = T,sep = "\t",quote = F)
  #write.table(cbind(genomes),paste("~/Documents/Post-doc/Marion_distib/ratio_par_",i,".txt",sep=""),col.names = F,row.names = T,sep = "\t",quote = F)
  
  plot(subgenomes~genomes,ylab=expression(CR/(CR+CO)),pch=21,bg=adjustcolor("gray70",alpha.f = 0),col="gray70",
       ylim=c(0,1),xlim=c(0,1),xlab=expression(CG/(CG+CO)),main=titre,las=1,cex.axis=1.3,cex.lab=1.3,cex.main=1.2)
  p<-dataEllipse(genomes,subgenomes,center.pch = F,levels=0.95,lty=4,col="red",plot.points =F,draw = F)
  points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="red",lty=2)
  #points(c(1:0),c(0:1),type="l",col="red",lty=2)
  #abline(lm(subgenomes~genomes-1),col="red",lty=2)
  xpr<-genomes*cos((pi/4))-subgenomes*sin((pi/4))
  ypr<-genomes*sin((pi/4))+subgenomes*cos((pi/4))
  pxpr<-p[,1]*cos((pi/4))-p[,2]*sin((pi/4))
  pypr<-p[,1]*sin((pi/4))+p[,2]*cos((pi/4))
  list<-names(xpr[which(xpr < min(pxpr) | xpr > max(pxpr))])
  if(i=="F"){listF<-names(xpr[which(xpr < min(pxpr) | xpr > max(pxpr))])}else{if(i=="L"){listL<-names(xpr[which(xpr < min(pxpr) | xpr > max(pxpr))])}else{listR<-names(xpr[which(xpr < min(pxpr) | xpr > max(pxpr))])}}

  print(length(list))
  #write.table(list,paste("~/Desktop/Marion_distib/cand_",i,".txt",sep=""),col.names = F,row.names = F,quote = F)
  points(subgenomes[which(names(subgenomes)%in%list)]~genomes[which(names(genomes)%in%list)],bg="red2",pch=21)
  points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="black",lty=2)
  points(seq(0,1,0.1),rep(0.5,11),type="l",col="black",lty=2)
  p<-dataEllipse(genomes,subgenomes,center.pch = F,levels=0.95,lty=4,col="red2",lwd = 2,plot.points =F)

  #coef<-round(cor.test(subgenomes,genomes,method="spearman")$estimate,2)
  coef<-round(summary(lm(subgenomes~genomes))$coefficients[2,1],2)
  #coef2<-round(summary(lm(subgenomes~genomes))$coefficients[1,1],2)
  #reg<-seq(0,1,0.1)*coef+coef2
  #points(reg~seq(0,1,1/(length(reg)-1)),type="l",col="gray30")
  abline(lm(subgenomes~genomes),xlim=c(0,1),col="gray40")
  text(0.0,1.2,expression(italic(beta)),col="red",pos=4,lwd=1.3)
  #text(0.02,1.2,paste(" = ",beta,sep=""),col="red",pos=4,lwd=1.3)
  mtext(bquote(italic(beta)*" = "*.(coef)),col="black",cex = 0.8,side=3,adj=0)
  letter<-sub("[a-z].*","",titre)
  # titre2<-c("No difference #1", "Legacy #2","Reverse #3","Intermediate #4","Compensatory-drift #5","Dominance CG #6a","Dominance CO #6b","Transgressive #7")
  # cat<-read.table(paste("~/Documents/Post-doc/Marion_distib/pattern/category_",letter,".txt",sep=""),header = T)
  # 
  # par(mfrow=c(3,3),mar=c(5,5,3,3))
  # plot(subgenomes~genomes,ylab=expression(Cbp[Cg]/(Cbp[Cg]+Cbp[Co])),pch=21,bg=adjustcolor("gray80",alpha.f = 0),col="gray80",
  #      ylim=c(0,1),xlim=c(0,1),xlab=expression(CG/(CG+CO)),main=titre,las=1,cex.axis=1.3,cex.lab=1.3,cex.main=1.2)
  # points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="black",lty=2)
  # points(seq(0,1,0.1),rep(0.5,11),type="l",col="black",lty=2)
  # coef<-round(summary(lm(subgenomes~genomes))$coefficients[2,1],2)
  # coef2<-round(summary(lm(subgenomes~genomes))$coefficients[1,1],2)
  # reg<-seq(0,1,0.1)*coef+coef2
  # points(reg~seq(0,1,1/(length(reg)-1)),type="l",col="gray20",lwd=2)
  # mtext(bquote(italic(beta)*" = "*.(coef)),col="black",cex = 0.8,side=3,adj=0)
  # l<-0
  # for(i in c(7,1,2,3,4,5,6,8)){
  #   l<-l+1
  #   plot(subgenomes~genomes,ylab=expression(Cbp[Cg]/(Cbp[Cg]+Cbp[Co])),pch=21,bg=adjustcolor("gray80",alpha.f = 0),col="gray80",
  #        ylim=c(0,1),xlim=c(0,1),xlab=expression(CG/(CG+CO)),main=titre2[l],las=1,cex.axis=1.3,cex.lab=1.3,cex.main=1.2)
  #   points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="black",lty=2)
  #   points(seq(0,1,0.1),rep(0.5,11),type="l",col="black",lty=2)
  #   coef<-round(summary(lm(subgenomes~genomes))$coefficients[2,1],2)
  #   coef2<-round(summary(lm(subgenomes~genomes))$coefficients[1,1],2)
  #   reg<-seq(0,1,0.1)*coef+coef2
  #   points(reg~seq(0,1,1/(length(reg)-1)),type="l",col="gray20",lwd=2)
  #   mtext(bquote(italic(beta)*" = "*.(coef)),col="black",cex = 0.8,side=3,adj=0)
  #   if(i!=8){points(subgenomes[which(names(subgenomes)%in%cat$transcript[which(cat$Category==i)])]~genomes[which(names(genomes)%in%cat$transcript[which(cat$Category==i)])],bg=color[i],pch=21)
  #     coef<-round(summary(lm(subgenomes[which(names(subgenomes)%in%cat$transcript[which(cat$Category==i)])]~genomes[which(names(genomes)%in%cat$transcript[which(cat$Category==i)])]))$coefficients[2,1],2)
  #     coef2<-round(summary(lm(subgenomes[which(names(subgenomes)%in%cat$transcript[which(cat$Category==i)])]~genomes[which(names(genomes)%in%cat$transcript[which(cat$Category==i)])]))$coefficients[1,1],2)
  #     reg<-seq(0,1,0.01)*coef+coef2
  #     if(i==4){points(reg[-which(reg < 0 | reg > 1)]~seq(0,1,1/(length(reg)-1))[-which(reg < 0 | reg > 1)],type="l",col=color[i],lwd=2)}else{
  #       points(reg~seq(0,1,1/(length(reg)-1)),type="l",col=color[i],lwd=2)}}
  #   if(i==8){
  #     list<-dir("~/Documents/Post-doc/Marion_distib/pattern/",pattern = paste("Transgressive.*_",letter,".txt",sep=""))
  #     up<-read.table(paste("~/Documents/Post-doc/Marion_distib/pattern/",list[6],sep=""))
  #     down<-read.table(paste("~/Documents/Post-doc/Marion_distib/pattern/",list[3],sep=""))
  #     CGup<-read.table(paste("~/Documents/Post-doc/Marion_distib/pattern/",list[4],sep=""))
  #     CGdown<-read.table(paste("~/Documents/Post-doc/Marion_distib/pattern/",list[1],sep=""))
  #     COup<-read.table(paste("~/Documents/Post-doc/Marion_distib/pattern/",list[5],sep=""))
  #     COdown<-read.table(paste("~/Documents/Post-doc/Marion_distib/pattern/",list[2],sep=""))
  #     points(subgenomes[which(names(subgenomes)%in%up$V1 | names(subgenomes)%in%down$V1)]~genomes[which(names(genomes)%in%up$V1 | names(subgenomes)%in%down$V1)],bg="gray50",pch=21)
  #     points(subgenomes[which(names(subgenomes)%in%CGup$V1 | names(subgenomes)%in%COup$V1)]~genomes[which(names(genomes)%in%CGup$V1 | names(subgenomes)%in%COup$V1)],bg="gray50",pch=21)
  #     points(subgenomes[which(names(subgenomes)%in%CGdown$V1 | names(subgenomes)%in%COdown$V1)]~genomes[which(names(genomes)%in%CGdown$V1 | names(subgenomes)%in%COdown$V1)],bg="white",pch=21)
  #     coef<-round(summary(lm(subgenomes[which(names(subgenomes)%in%cat$transcript[which(cat$Category==i)])]~genomes[which(names(genomes)%in%cat$transcript[which(cat$Category==i)])]))$coefficients[2,1],2)
  #     coef2<-round(summary(lm(subgenomes[which(names(subgenomes)%in%cat$transcript[which(cat$Category==i)])]~genomes[which(names(genomes)%in%cat$transcript[which(cat$Category==i)])]))$coefficients[1,1],2)
  #     reg<-seq(0,1,0.01)*coef+coef2
  #     points(reg~seq(0,1,1/(length(reg)-1)),type="l",col=color[i],lwd=2)
  #   }
    
    #mtext(bquote(italic(beta)[pro]*" = "*.(coef)),cex = 0.8,side=3,adj=1,col=color[i])
}


################################# conv_div ###########################

domh1<-function(x){
  mean<-c(mean(x[which(F_d_phased_t$genome==0)],na.rm = T),mean(x[which(F_d_phased_t$genome==4)],na.rm = T),mean(x[which(F_d_phased_t$genome==2)],na.rm = T))
  mod2<-lm(mean[c(1,3)]~c(0,2))
  pred2<-mean(mod2$fitted.values)
  obs2<-mean[2]
  
  sign<-mod2$coefficients[2]
  if(sign<0){
    coeff2<-(pred2-obs2)/pred2} else {
      coeff2<-(obs2-pred2)/pred2
    }
}

d1CG<-function(x){
  mean<-c(mean(x[which(F_d_phased_t$genome==0)],na.rm=T),mean(x[which(F_d_phased_t$genome==2)],na.rm = T),mean(x[which(F_d_phased_t$genome==3)],na.rm = T),mean(x[which(F_d_phased_t$genome==1)],na.rm = T))
  d1CgCo<-mean[4]-mean[3]
  d1CgCO<-mean[4]-mean[2]
  d1CGCO<-mean[1]-mean[2]
  return(c(d1CgCo,d1CgCO,d1CGCO))
}

d2CG<-function(x){
  mean<-c(mean(x[which(F_d_phased_t$genome==0)],na.rm=T),mean(x[which(F_d_phased_t$genome==2)],na.rm = T),mean(x[which(F_d_phased_t$genome==3)],na.rm = T),mean(x[which(F_d_phased_t$genome==1)],na.rm = T))
  d2CoCg<-mean[3]-mean[4]
  d2CoCG<-mean[3]-mean[1]
  d2CGCO<-mean[1]-mean[2]
  return(c(d2CoCg,d2CoCG,d2CGCO))
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
out_mod<-c()
test_prop_d<-c()
test_prop_sub<-c()
ref<-c()
cor2<-c()
par(mfrow=c(2,3),oma=c(0,0,2,0),mar=c(4,4,4,1))
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
  #F_d_phased<-F_d_phased[,-c(4,5)]
  genome<-c(rep(0,1),rep(2,1),rep(1,1),rep(3,2))
  F_d_phased<-na.omit(F_d_phased)
  #F_d_phased<-as.data.frame(t(F_d_phased))
  d<-DGEList(counts=F_d_phased,group=genome)
  tmmy<-calcNormFactors(d,method = "TMM")
  d<- estimateCommonDisp(tmmy)
  relation<-as.factor(rownames(d$samples))
  design<-model.matrix(~relation,datat=d)
  
  dGLM<-estimateGLMCommonDisp(d,design)
  
  
  F_d_phased<-d$pseudo.counts
  genome<-c(4,0,2,3,1) ## CR, CG, CO, Co, Cg
  F_d_phased_t<-t.default(F_d_phased)
  F_d_phased_t<-cbind.data.frame(genome,F_d_phased_t)
  if(i=="F"){titre<-"Flowers"}else{if(i=="L"){titre<-"Leaves"}else{titre<-"Roots"}}
  
  testh1<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh1(x))
  # d1<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)d1CG(x))
  # d1<-d1[,which(abs(d1[1,])<=10000 & abs(d1[2,])<=10000)]
  # sub1<-d1[,which(abs(d1[1,])<=2000 & abs(d1[2,])<=2000)]
  # 
  # p1_full<-lm(d1[1,]~0+d1[2,])$coefficients
  # p1_full_abs<-lm(abs(d1[1,])~0+abs(d1[2,]))$coefficients
  # size1full<-round(length(d1[1,])/length(testh1),2)
  # p1_2000<-lm(sub1[1,]~0+sub1[2,])$coefficients
  # p1_2000abs<-lm(abs(sub1[1,])~0+abs(sub1[2,]))$coefficients
  # size1_2000<-round(length(sub1[1,])/length(testh1),2)
  # out_mod_temp<-cbind(p1_full,p1_full_abs,p1_2000,p1_2000abs,size1full,size1_2000)
  # out_mod<-rbind(out_mod,out_mod_temp)
  # test_d1<-prop.test(length(which(abs(d1[1,])<abs(d1[2,]))),ncol(d1),p=0.5)
  # test_sub1<-prop.test(length(which(abs(sub1[1,])<abs(sub1[2,]))),ncol(sub1),p=0.5)
  
  #test2<-t(test)
  #colnames(test2)<-c("coeff","coeff2")
  #hist(test,breaks=100,freq = F)
  hist(testh1[which(testh1 < 1 & testh1 > -1)],ylim=c(0,2.5),xlim=c(-1,1),breaks=50,freq = F, xlab = expression(CR),main = paste(titre," (",length(which(testh1!=0)),")",sep=""))
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
  
  # testh2<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh2(x))
  # 
  # d2<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)d2CG(x))
  # d2<-d2[,which(abs(d2[1,])<=10000 & abs(d2[2,])<=10000)]
  # sub2<-d2[,which(abs(d2[1,])<=2000 & abs(d2[2,])<=2000)]
  # size2full<-round(length(d2[1,])/length(testh2),2)
  # p2_full<-lm(d2[1,]~0+d2[2,])$coefficients
  # p2_full_abs<-lm(abs(d2[1,])~0+abs(d2[2,]))$coefficients
  # p2_2000<-lm(sub2[1,]~0+sub2[2,])$coefficients
  # p2_2000abs<-lm(abs(sub2[1,])~0+abs(sub2[2,]))$coefficients
  # size2_2000<-round(length(sub2[1,])/length(testh2),2)
  # out_mod_temp<-cbind(p2_full,p2_full_abs,p2_2000,p2_2000abs,size2full,size2_2000)
  # out_mod<-rbind(out_mod,out_mod_temp)
  # test_d2<-prop.test(length(which(abs(d2[1,])<abs(d2[2,]))),ncol(d2),p=0.5)
  # test_sub2<-prop.test(length(which(abs(sub2[1,])<abs(sub2[2,]))),ncol(sub2),p=0.5)
  # test_temp_d<-cbind(ncol(d2),test_d2$estimate,test_d2$conf.int[1],test_d2$conf.int[2],test_d2$statistic,test_d2$p.value,
  #                    ncol(d1),test_d1$estimate,test_d1$conf.int[1],test_d1$conf.int[2],test_d1$statistic,test_d1$p.value)
  # test_prop_d<-rbind(test_prop_d,test_temp_d)
  # test_temp_sub<-cbind(ncol(sub2),test_sub2$estimate,test_sub2$conf.int[1],test_sub2$conf.int[2],test_sub2$statistic,test_sub2$p.value,
  #                      ncol(sub1),test_sub1$estimate,test_sub1$conf.int[1],test_sub1$conf.int[2],test_sub1$statistic,test_sub1$p.value)
  # test_prop_sub<-rbind(test_prop_sub,test_temp_sub)
  # 
  # 
  # #test2<-t(test)
  # #colnames(test2)<-c("coeff","coeff2")
  # #hist(test,breaks=100,freq = F)
  # 
  # hist(testh2[which(testh2 < 1 & testh2 > -1)],ylim = c(0,2.5),xlim=c(-1,1),breaks=50,freq = F, xlab = expression(italic(S[Co])),main = paste(titre," (",length(which(testh2!=0)),")",sep=""))
  # abline(v=0, col="red2",lwd=2,lty=2)
  # abline(v=median(testh2,na.rm=T),col="green",lwd=2,lty=3)
  # coeff_plus<-mean(testh2[which(testh2 <= 1 & testh2 > 0)])
  # coeff_moins<-mean(testh2[which(testh2 <0 & testh2 >= -1)])
  # ref2<-length(which(testh2<0))/length(testh2)
  # b_test<-binom.test(c(length(which(testh2<0)),length(which(testh2>0))),p=0.5)
  # if(b_test$p.value>=0.05){out<-paste(round(ref2,2),"ns",sep=" ")}else{if(b_test$p.value<0.001){out<-paste(round(ref2,2),"***",sep=" ")}else{if(b_test$p.value<0.01){out<-paste(round(ref2,2),"**",sep=" ")}else{out<-paste(round(ref2,2),"*",sep=" ")}}}
  # text(-1,1.5,lab=out,pos=4)#," / ",round(coeff_moins,2),sep=""),pos=4)
  # text(1,1.5,lab=paste(round(1-ref2,2)),pos=2)#," / ",round(coeff_plus,2),sep=""),pos=2)
  # #theo<-rnorm(length(test),0,sd(test))
  # #lines(density(theo),col="red",lwd=2)
  # #length(test[which(test <= 1 & test >= -1)])/length(test)
  # ref[length(ref)+1]<-median(testh2,na.rm=T)
  # cuth1<-names(head(sort(abs(testh1),decreasing  =T),50))
  # cuth2<-names(head(sort(abs(testh2),decreasing  =T),50))
  # cut<-unique(c(cuth1,cuth2))
  # h1<-testh1[-which(names(testh1)%in%cut)]
  # h2<-testh2[-which(names(testh2)%in%cut)]
  # h1<-h1[order(names(h1))]
  # h2<-h2[order(names(h2))]
  # cor<-cor.test(h1,h2,method="spearman")
  # cor<-c(cor$estimate,cor$parameter,cor$statistic,cor$p.value)
  # cor2<-rbind.data.frame(cor2,cor)
  # 
  # plot(h2~h1,xlab=expression(italic(S[Cg])),asp=1,ylab=expression(italic(S[Co])),main=paste(titre," (",length(h1),")",sep=""))
  # abline(lm(h2~h1),lty=2,col="red")
  #length(which(test<0 & test >= -1))/length(test[which(test <= 1 & test >= -1)])
  #print(wilcox.test(testh1,testh2,paired = F))
  
}
#mtext(side = c(3),text= "All genes",font=2,outer = T)
cor2
out_mod<-round(out_mod,2)
rownames(out_mod)<-c("F_Cg","F_Co","L_Cg","L_Co","R_Cg","R_Co")
colnames(out_mod)<-c("10K","10K_abs","2K","2K_abs","Size_10K","Size_2K")

test_prop_d
test_prop_sub
###################### DE gene #####################
out_s<-c()
out_d2<-c()
out_d1<-c()
out_d31<-c()
out_d32<-c()
test_prop_sub<-c()
out_mod<-c()
tissue<-c("F","L","R")
ref<-c()
cor2<-c()
#par(mfrow=c(3,1),oma=c(0,0,2,0),mar=c(4,4,4,1))
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
  
  genome<-c(rep(0,1),rep(2,1),rep(1,1),rep(3,2))
  F_d_phased<-na.omit(F_d_phased)
  #F_d_phased<-as.data.frame(t(F_d_phased))
  d<-DGEList(counts=F_d_phased,group=genome)
  tmmy<-calcNormFactors(d,method = "TMM")
  d<- estimateCommonDisp(tmmy)
  relation<-as.factor(rownames(d$samples))
  design<-model.matrix(~relation,datat=d)
  
  dGLM<-estimateGLMCommonDisp(d,design)
  
  
  F_d_phased<-d$pseudo.counts
  genome<-c(4,0,2,3,1) ## CR, CG, CO, Co, Cg
  F_d_phased_t<-t.default(F_d_phased)
  F_d_phased_t<-cbind.data.frame(genome,F_d_phased_t)
  if(i=="F"){titre<-"Flowers"}else{if(i=="L"){titre<-"Leaves"}else{titre<-"Roots"}}
  
  testh1<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh1(x))
  #d1<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)d1CG(x))
  #d1<-d1[,which(abs(d1[1,])<=10000 & abs(d1[2,])<=10000)]
  #d1<-d1[,which(abs(d1[1,])<=2000 & abs(d1[2,])<=2000)]
  
  
  
  #p1_full<-lm(d1[1,]~0+d1[2,])$coefficients
  #p1_full_abs<-lm(abs(d1[1,])~0+abs(d1[2,]))$coefficients
  #size1full<-round(length(d1[1,])/length(testh1),2)
  #p1_2000<-lm(sub1[1,]~0+sub1[2,])$coefficients
  #p1_2000abs<-lm(abs(sub1[1,])~0+abs(sub1[2,]))$coefficients
  #size1_2000<-round(length(sub1[1,])/length(testh1),2)
  
  #out_mod_temp<-cbind(p1_full,p1_full_abs,p1_2000,p1_2000abs,size1full,size1_2000)
  #out_mod<-rbind(out_mod,out_mod_temp)
  
  # prop1<-(abs(d1[1,])-abs(d1[2,]))/max(abs(d1[1,]),abs(d1[2,]))
  # propbis<-(abs(d1[1,])-abs(d1[3,]))/max(abs(d1[1,]),abs(d1[3,]))
  # dist<-c(t(d1[3,]))
  # distr<-rank(abs(dist))/length(dist)
  # print(cor.test(propbis*-1,abs(distr),method = "spearman"))
  # for(z in c(0,0.000001,0.00001,0.0001,0.001,0.005,0.01,0.05,0.1,0.15,0.2)){
  #   test_d1<-prop.test(length(which(prop1 < 0& abs(prop1)> z )),length(which(abs(prop1)>z)),p=0.5)
  #   test_temp_d1<-cbind(i,length(which(abs(prop1)>z)),length(which(abs(prop1)>z))/ncol(d1),z,test_d1$estimate,test_d1$conf.int[1],test_d1$conf.int[2],test_d1$statistic,test_d1$p.value)
  #   out_d1<-rbind.data.frame(out_d1,test_temp_d1)
  # }
  
  
  #test2<-t(test)
  #colnames(test2)<-c("coeff","coeff2")
  #hist(test,breaks=100,freq = F)
  hist(testh1[which(testh1 < 1 & testh1 > -1)],ylim=c(0,2.5),xlim=c(-1,1),breaks=50,freq = F, xlab = expression(CR),main = paste(titre," (",length(which(testh1!=0)),")",sep=""))
  abline(v=0, col="red2",lwd=2,lty=2)
  abline(v=median(testh1,na.rm=T),col="green",lwd=2,lty=3)
  coeff_plus<-mean(testh1[which(testh1 <= 1 & testh1 > 0)])
  coeff_moins<-mean(testh1[which(testh1 <0 & testh1 >= -1)])
  ref2<-length(which(testh1<0))/length(testh1)
  b_test1<-binom.test(c(length(which(testh1<0)),length(which(testh1>0))),p=0.5)
  if(b_test1$p.value>=0.05){out<-paste(round(ref2,2),"ns",sep=" ")}else{if(b_test1$p.value<0.001){out<-paste(round(ref2,2),"***",sep=" ")}else{if(b_test1$p.value<0.01){out<-paste(round(ref2,2),"**",sep=" ")}else{out<-paste(round(ref2,2),"*",sep=" ")}}}
  text(-1,1.5,lab=out,pos=4)#," / ",round(coeff_moins,2),sep=""),pos=4)
  text(1,1.5,lab=paste(round(1-ref2,2)),pos=2)#," / ",round(coeff_plus,2),sep=""),pos=2)
  out_s1<-cbind(i,"Cg",length(testh1),b_test1$estimate,b_test1$conf.int[1],b_test1$conf.int[2],b_test1$p.value,median(testh1,na.rm=T))
  #theo<-rnorm(length(test),0,sd(test))
  #lines(density(theo),col="red",lwd=2)
  ref[length(ref)+1]<-median(testh1,na.rm=T)
  
  ################################ H2_flower ##################
  #abs(mean[3]-mean[1])
  
  # testh2<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)domh2(x))
  # #test2<-t(test)
  # #colnames(test2)<-c("coeff","coeff2")
  # #hist(test,breaks=100,freq = F)
  # 
  # hist(testh2[which(testh2 < 1 & testh2 > -1)],ylim = c(0,2.5),xlim=c(-1,1),breaks=50,freq = F, xlab = expression(italic(S[Co])),main = paste(titre," (",length(which(testh1!=0)),")",sep=""))
  # 
  # d2<-apply(F_d_phased_t[,c(2:ncol(F_d_phased_t))],2,function(x)d2CG(x))
  # p2_full<-lm(d2[1,]~0+d2[2,])$coefficients
  # p2_full_abs<-lm(abs(d2[1,])~0+abs(d1[2,]))$coefficients
  # size1full<-round(length(d1[1,])/length(testh1),2)
  # p1_2000<-lm(sub1[1,]~0+sub1[2,])$coefficients
  # p1_2000abs<-lm(abs(sub1[1,])~0+abs(sub1[2,]))$coefficients
  # size1_2000<-round(length(sub1[1,])/length(testh1),2)
  # out_mod_temp<-cbind(p1_full,p1_full_abs,p1_2000,p1_2000abs,size1full,size1_2000)
  # out_mod<-rbind(out_mod,out_mod_temp)
  # 
  # prop2<-(abs(d2[1,])-abs(d2[2,]))/max(abs(d2[1,]),abs(d2[2,]))
  # for(z in c(0,0.000001,0.00001,0.0001,0.001,0.005,0.01,0.05,0.1,0.15,0.2)){
  #   test_d2<-prop.test(length(which(prop2 < 0 & abs(prop2)> z )),length(which(abs(prop2)>z)),p=0.5)
  #   test_temp_d2<-cbind(i,length(which(abs(prop2)>z)),length(which(abs(prop2)>z))/ncol(d2),z,test_d2$estimate,test_d2$conf.int[1],test_d2$conf.int[2],test_d2$statistic,test_d2$p.value)
  #   out_d2<-rbind.data.frame(out_d2,test_temp_d2)
  # }
  # 
  # prop31<-(abs(d1[1,])-abs(d1[3,]))/max(abs(d1[1,]),abs(d1[3,]))
  # for(z in c(0,0.000001,0.00001,0.0001,0.001,0.005,0.01,0.05,0.1,0.15,0.2)){
  #   test_d31<-prop.test(length(which(prop31 < 0& abs(prop31)> z )),length(which(abs(prop31)>z)),p=0.5)
  #   test_temp_d31<-cbind(i,length(which(abs(prop31)>z)),length(which(abs(prop31)>z))/ncol(d1),z,test_d31$estimate,test_d31$conf.int[1],test_d31$conf.int[2],test_d31$statistic,test_d31$p.value)
  #   out_d31<-rbind.data.frame(out_d31,test_temp_d31)
  # } 
  # 
  # prop32<-(abs(d2[1,])-abs(d2[3,]))/max(abs(d2[1,]),abs(d2[3,]))
  # for(z in c(0,0.000001,0.00001,0.0001,0.001,0.005,0.01,0.05,0.1,0.15,0.2)){
  #   test_d32<-prop.test(length(which(prop32 < 0& abs(prop32)> z )),length(which(abs(prop32)>z)),p=0.5)
  #   test_temp_d32<-cbind(i,length(which(abs(prop32)>z)),length(which(abs(prop32)>z))/ncol(d2),z,test_d32$estimate,test_d32$conf.int[1],test_d32$conf.int[2],test_d32$statistic,test_d32$p.value)
  #   out_d32<-rbind.data.frame(out_d32,test_temp_d32)
  # } 
  # 
  # test_d2<-prop.test(length(which(abs(d2[1,])<abs(d2[2,]))),ncol(d2),p=0.5)
  # test_sub2<-prop.test(length(which(abs(sub2[1,])<abs(sub2[2,]))),ncol(sub2),p=0.5)
  # test_temp_d<-cbind(ncol(d2),test_d2$estimate,test_d2$conf.int[1],test_d2$conf.int[2],test_d2$statistic,test_d2$p.value,
  #                    ncol(d1),test_d1$estimate,test_d1$conf.int[1],test_d1$conf.int[2],test_d1$statistic,test_d1$p.value)
  # test_prop_d<-rbind(test_prop_d,test_temp_d)
  # test_temp_sub<-cbind(ncol(sub2),test_sub2$estimate,test_sub2$conf.int[1],test_sub2$conf.int[2],test_sub2$statistic,test_sub2$p.value,
  #                      ncol(sub1),test_sub1$estimate,test_sub1$conf.int[1],test_sub1$conf.int[2],test_sub1$statistic,test_sub1$p.value)
  # test_prop_sub<-rbind(test_prop_sub,test_temp_sub)
  # 
  # abline(v=0, col="red2",lwd=2,lty=2)
  # abline(v=median(testh2,na.rm=T),col="green",lwd=2,lty=3)
  # coeff_plus<-mean(testh2[which(testh2 <= 1 & testh2 > 0)])
  # coeff_moins<-mean(testh2[which(testh2 <0 & testh2 >= -1)])
  # ref2<-length(which(testh2<0))/length(testh2)
  # b_test2<-binom.test(c(length(which(testh2<0)),length(which(testh2>0))),p=0.5)
  # if(b_test2$p.value>=0.05){out<-paste(round(ref2,2),"ns",sep=" ")}else{if(b_test2$p.value<0.001){out<-paste(round(ref2,2),"***",sep=" ")}else{if(b_test2$p.value<0.01){out<-paste(round(ref2,2),"**",sep=" ")}else{out<-paste(round(ref2,2),"*",sep=" ")}}}
  # text(-1,1.5,lab=out,pos=4)#," / ",round(coeff_moins,2),sep=""),pos=4)
  # text(1,1.5,lab=paste(round(1-ref2,2)),pos=2)#," / ",round(coeff_plus,2),sep=""),pos=2)
  # #theo<-rnorm(length(test),0,sd(test))
  # #lines(density(theo),col="red",lwd=2)
  # #length(test[which(test <= 1 & test >= -1)])/length(test)
  # ref[length(ref)+1]<-median(testh2,na.rm=T)
  # #print(wilcox.test(testh1,testh2,paired = F))
  # cuth1<-names(head(sort(abs(testh1),decreasing  =T),50))
  # cuth2<-names(head(sort(abs(testh2),decreasing  =T),50))
  # cut<-unique(c(cuth1,cuth2))
  # h1<-testh1[-which(names(testh1)%in%cut)]
  # h2<-testh2[-which(names(testh2)%in%cut)]
  # h1<-h1[order(names(h1))]
  # h2<-h2[order(names(h2))]
  # cor<-cor.test(h1,h2,method="spearman")
  # cor<-c(cor$estimate,cor$parameter,cor$statistic,cor$p.value)
  # cor2<-rbind.data.frame(cor2,cor)
  # out_s2<-cbind(i,"Co",length(testh2),b_test2$estimate,b_test2$conf.int[1],b_test2$conf.int[2],b_test2$p.value,median(testh2,na.rm=T))
  # out_s<-rbind.data.frame(out_s,out_s1,out_s2)
  # 
  # plot(h2~h1,xlab=expression(italic(S[Cg])),asp=1,ylab=expression(italic(S[Co])),main=paste(titre," (",length(h1),")",sep=""))
  # abline(lm(h2~h1),lty=2,col="red")
  #r_test<-cor[4]
  #if(r_test>=0.05){out<-paste(round(cor[1],2),"ns",sep=" ")}else{if(r_test<0.001){out<-paste(round(cor[1],2),"***",sep=" ")}else{if(r_test<0.01){out<-paste(round(cor[1],2),"**",sep=" ")}else{out<-paste(round(cor[1],2),"*",sep=" ")}}}
  #legend("topleft",paste("r = ",out,sep=""),bty = "n",y.intersp = 0)
  #text(-3.5,1.7,paste("r = ",out,sep=""),pos=4)
  #length(which(test<0 & test >= -1))/length(test[which(test <= 1 & test >= -1)])
  
}
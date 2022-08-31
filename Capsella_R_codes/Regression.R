########################         Regression         ########################
## 1. CbpCg / (CbpCg + CbpCo) ~ CG / (CG + CO)
## 2. CR / (CR + CO) ~ CG / (CG + CO)

head(S_Flower)
head(S_Leaf)
head(S_Root)

require(edgeR)
require(car)
require(heatmap3)

Cbp_subgenome_F <- S_Flower$mean_Cbp_Cg_F / (S_Flower$mean_Cbp_Cg_F + S_Flower$mean_Cbp_Co_F)
var(Cbp_subgenome_F) # variance of data
Parents_genome_F <- S_Flower$mean_CG_F / (S_Flower$mean_CG_F + S_Flower$mean_CO_F)
var(Parents_genome_F)
var.test(Cbp_subgenome_F, Parents_genome_F) # F Test To Compare Two Variances

opar <- par
par(mfrow=c(2,3))
plot(Cbp_subgenome_F~Parents_genome_F,ylab=expression(CbpCg/(CbpCg+CbpCo)),pch=21,bg=adjustcolor("gray70",alpha.f = 0),col="gray70",
     ylim=c(0,1),xlim=c(0,1),xlab=expression(CG/(CG+CO)),main=paste0("Flower", " (", dim(S_Flower)[1], ")") )
p<- dataEllipse(Cbp_subgenome_F,Parents_genome_F,center.pch = F,levels=0.95,lty=4,col="red",plot.points =F,draw = F)
points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="red",lty=2)
xpr<-Parents_genome_F*cos((pi/4))-Cbp_subgenome_F*sin((pi/4))
ypr<-Parents_genome_F*sin((pi/4))+Cbp_subgenome_F*cos((pi/4))
pxpr<-p[,1]*cos((pi/4))-p[,2]*sin((pi/4))
pypr<-p[,1]*sin((pi/4))+p[,2]*cos((pi/4))
list<- which(xpr < min(pxpr) | xpr > max(pxpr))
length(list)
points(Cbp_subgenome_F[list]~Parents_genome_F[list],bg="red2",pch=21)
points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="black",lty=2)
points(seq(0,1,0.1),rep(0.5,11),type="l",col="black",lty=2)
p<-dataEllipse(Cbp_subgenome_F,Parents_genome_F,center.pch = F,levels=0.95,lty=4,col="red2",lwd = 2,plot.points =F)

#coef<-round(cor.test(subgenomes,genomes,method="spearman")$estimate,2)
coef<-round(summary(lm(Cbp_subgenome_F~Parents_genome_F))$coefficients[2,1],2)
#coef2<-round(summary(lm(subgenomes~genomes))$coefficients[1,1],2)
#reg<-seq(0,1,0.1)*coef+coef2
#points(reg~seq(0,1,1/(length(reg)-1)),type="l",col="gray30")
abline(lm(Cbp_subgenome_F~Parents_genome_F),xlim=c(0,1),col="gray40")
text(0.0,1.2,expression(italic(beta)),col="red",pos=4,lwd=1.3)
#text(0.02,1.2,paste(" = ",beta,sep=""),col="red",pos=4,lwd=1.3)
mtext(bquote(italic(beta)*" = "*.(coef)),col="black",cex = 0.8,side=3,adj=0)
letter<-sub("[a-z].*","",titre)
# titre2<-c("No difference #1", "Legacy #2","Reverse #3

## Leaf ##
Cbp_subgenome_L <- S_Leaf$mean_Cbp_Cg_L / (S_Leaf$mean_Cbp_Cg_L + S_Leaf$mean_Cbp_Co_L)
var(Cbp_subgenome_L) # variance of data
Parents_genome_L <- S_Leaf$mean_CG_L / (S_Leaf$mean_CG_L + S_Leaf$mean_CO_L)
var(Parents_genome_L)
var.test(Cbp_subgenome_L, Parents_genome_L) # F Test To Compare Two Variances
plot(Cbp_subgenome_L~Parents_genome_L,ylab=expression(CbpCg/(CbpCg+CbpCo)),pch=21,bg=adjustcolor("gray70",alpha.f = 0),col="gray70",
     ylim=c(0,1),xlim=c(0,1),xlab=expression(CG/(CG+CO)),main=paste0("Leaf", " (", dim(S_Leaf)[1], ")") )
p<- dataEllipse(Cbp_subgenome_L,Parents_genome_L,center.pch = F,levels=0.95,lty=4,col="red",plot.points =F,draw = F)
points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="red",lty=2)
xpr<-Parents_genome_L*cos((pi/4))-Cbp_subgenome_L*sin((pi/4))
ypr<-Parents_genome_L*sin((pi/4))+Cbp_subgenome_L*cos((pi/4))
pxpr<-p[,1]*cos((pi/4))-p[,2]*sin((pi/4))
pypr<-p[,1]*sin((pi/4))+p[,2]*cos((pi/4))
list<- which(xpr < min(pxpr) | xpr > max(pxpr))
length(list)
points(Cbp_subgenome_L[list]~Parents_genome_L[list],bg="red2",pch=21)
points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="black",lty=2)
points(seq(0,1,0.1),rep(0.5,11),type="l",col="black",lty=2)
p<-dataEllipse(Cbp_subgenome_L,Parents_genome_L,center.pch = F,levels=0.95,lty=4,col="red2",lwd = 2,plot.points =F)

#coef<-round(cor.test(subgenomes,genomes,method="spearman")$estimate,2)
coef<-round(summary(lm(Cbp_subgenome_L~Parents_genome_L))$coefficients[2,1],2)
#coef2<-round(summary(lm(subgenomes~genomes))$coefficients[1,1],2)
#reg<-seq(0,1,0.1)*coef+coef2
#points(reg~seq(0,1,1/(length(reg)-1)),type="l",col="gray30")
abline(lm(Cbp_subgenome_L~Parents_genome_L),xlim=c(0,1),col="gray40")
text(0.0,1.2,expression(italic(beta)),col="red",pos=4,lwd=1.3)
#text(0.02,1.2,paste(" = ",beta,sep=""),col="red",pos=4,lwd=1.3)
mtext(bquote(italic(beta)*" = "*.(coef)),col="black",cex = 0.8,side=3,adj=0)
# titre2<-c("No difference #1", "Legacy #2","Reverse #3

## Root ##
Cbp_subgenome_R <- S_Root$mean_Cbp_Cg_R / (S_Root$mean_Cbp_Cg_R + S_Root$mean_Cbp_Co_R)
var(Cbp_subgenome_R) # variance of data
Parents_genome_R <- S_Root$mean_CG_R / (S_Root$mean_CG_R + S_Root$mean_CO_R)
var(Parents_genome_R)
var.test(Cbp_subgenome_R, Parents_genome_R) # F Test To Compare Two Variances
plot(Cbp_subgenome_R~Parents_genome_R,ylab=expression(CbpCg/(CbpCg+CbpCo)),pch=21,bg=adjustcolor("gray70",alpha.f = 0),col="gray70",
     ylim=c(0,1),xlim=c(0,1),xlab=expression(CG/(CG+CO)),main=paste0("Root", " (", dim(S_Root)[1], ")") )
p<- dataEllipse(Cbp_subgenome_R,Parents_genome_R,center.pch = F,levels=0.95,lty=4,col="red",plot.points =F,draw = F)
points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="red",lty=2)
xpr<-Parents_genome_R*cos((pi/4))-Cbp_subgenome_R*sin((pi/4))
ypr<-Parents_genome_R*sin((pi/4))+Cbp_subgenome_R*cos((pi/4))
pxpr<-p[,1]*cos((pi/4))-p[,2]*sin((pi/4))
pypr<-p[,1]*sin((pi/4))+p[,2]*cos((pi/4))
list<- which(xpr < min(pxpr) | xpr > max(pxpr))
length(list)
points(Cbp_subgenome_R[list]~Parents_genome_R[list],bg="red2",pch=21)
points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="black",lty=2)
points(seq(0,1,0.1),rep(0.5,11),type="l",col="black",lty=2)
p<-dataEllipse(Cbp_subgenome_R,Parents_genome_R,center.pch = F,levels=0.95,lty=4,col="red2",lwd = 2,plot.points =F)

#coef<-round(cor.test(subgenomes,genomes,method="spearman")$estimate,2)
coef<-round(summary(lm(Cbp_subgenome_R~Parents_genome_R))$coefficients[2,1],2)
#coef2<-round(summary(lm(subgenomes~genomes))$coefficients[1,1],2)
#reg<-seq(0,1,0.1)*coef+coef2
#points(reg~seq(0,1,1/(length(reg)-1)),type="l",col="gray30")
abline(lm(Cbp_subgenome_R~Parents_genome_R),xlim=c(0,1),col="gray40")
text(0.0,1.2,expression(italic(beta)),col="red",pos=4,lwd=1.3)
#text(0.02,1.2,paste(" = ",beta,sep=""),col="red",pos=4,lwd=1.3)
mtext(bquote(italic(beta)*" = "*.(coef)),col="black",cex = 0.8,side=3,adj=0)
# titre2<-c("No difference #1", "Legacy #2","Reverse #3

##########     2. CR / (CR + CO) ~ CG / (CG + CO)    ##########
tissue<-c("F","L","R") # set tissue

# Loop start
for (i in tissue){
        print(i)
        # set a sub loop for viable name
        if (i == "F"){
                Tissue <- "Flower"
        }else{
                if (i == "L"){
                        Tissue <- "Leaf"
                }else{
                        Tissue <- "Root"
                }
        }
        
        tmm <- read.table(paste0("OutputData/TMM_CGCRCO_", i, ".txt"), header = T, sep = "\t")
        tmm$meanCG <- rowMeans(tmm[, which(sub("[0-9]_.*", "", colnames(tmm)) == "CG")])
        tmm$meanCR <- rowMeans(tmm[, which(sub("[0-9]_.*", "", colnames(tmm)) == "CR")])
        tmm$meanCO <- rowMeans(tmm[, which(sub("[0-9]_.*", "", colnames(tmm)) == "CO")])
        
        CR2CRCO <- tmm$meanCR / (tmm$meanCR + tmm$meanCO)
        CG2CGCO <- tmm$meanCG / (tmm$meanCG + tmm$meanCO)
        
        var.test(CR2CRCO, CG2CGCO) # F Test To Compare Two Variances
        plot(CR2CRCO ~ CG2CGCO, ylab=expression(CR/(CR+CO)),pch=21,bg=adjustcolor("gray70",alpha.f = 0),col="gray70",
             ylim=c(0,1),xlim=c(0,1),xlab=expression(CG/(CG+CO)),main=paste0(Tissue, " (", dim(tmm)[1], ")"))
        p<- dataEllipse(CR2CRCO,CG2CGCO,center.pch = F,levels=0.95,lty=4,col="red",plot.points =F,draw = F)
        points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="red",lty=2)
        xpr<-CG2CGCO*cos((pi/4))-CR2CRCO*sin((pi/4))
        ypr<-CG2CGCO*sin((pi/4))+CR2CRCO*cos((pi/4))
        pxpr<-p[,1]*cos((pi/4))-p[,2]*sin((pi/4))
        pypr<-p[,1]*sin((pi/4))+p[,2]*cos((pi/4))
        list<- which(xpr < min(pxpr) | xpr > max(pxpr))
        length(list)
        points(CR2CRCO[list]~CG2CGCO[list],bg="red2",pch=21)
        points(seq(0,1,0.1),seq(0,1,0.1),type="l",col="black",lty=2)
        points(seq(0,1,0.1),rep(0.5,11),type="l",col="black",lty=2)
        p<-dataEllipse(CR2CRCO,CG2CGCO,center.pch = F,levels=0.95,lty=4,col="red2",lwd = 2,plot.points =F)
        
        #coef<-round(cor.test(subgenomes,genomes,method="spearman")$estimate,2)
        coef<-round(summary(lm(CR2CRCO~CG2CGCO))$coefficients[2,1],2)
        #coef2<-round(summary(lm(subgenomes~genomes))$coefficients[1,1],2)
        #reg<-seq(0,1,0.1)*coef+coef2
        #points(reg~seq(0,1,1/(length(reg)-1)),type="l",col="gray30")
        abline(lm(CR2CRCO~CG2CGCO),xlim=c(0,1),col="gray40")
        text(0.0,1.2,expression(italic(beta)),col="red",pos=4,lwd=1.3)
        #text(0.02,1.2,paste(" = ",beta,sep=""),col="red",pos=4,lwd=1.3)
        mtext(bquote(italic(beta)*" = "*.(coef)),col="black",cex = 0.8,side=3,adj=0)
        
}

par(opar)
dev.off()


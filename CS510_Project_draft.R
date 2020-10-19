########################################################################################
# THE DATA USED FOR THESE ANALYSES IS CONFIDENTIAL & CANNOT BE SHARED
# GWAS Analysis in R
# CIDR ER / PR Subtype Analysis
# Subsets both Covariate & Dosage Files
# Creates MH & QQ plots
# Yesha Patel
# CS510 Project
########################################################################################

rm(list=ls())

# Attach R package statmod
library("statmod")

# Set working directories
Dosage_Directory = "F:/Project/R_Results/Imputed/"
Fam_Directory = "F:/Project/Files_yp/imputed/"
Cov_Directory = "F:/Project/Files_yp/imputed/"
Output_Directory= "F:/Project/R_Results/Imputed/"

filename ="chr10_p_noage_score_zero.gprobs"
fam_filename = "chr10.sample"

# Read Covariate files 
Covar = read.table(paste(Cov_Directory,"CIDR_ER_PR_R_cov.txt", sep=""),header=TRUE)
CovIDs=Covar[,2]

# Read Fam file (sample file) and make sure subject order is the same
Fam=read.table(paste(Fam_Directory,fam_filename,sep=""),skip=2,header=FALSE)
FamIDs=Fam[,2]

IDs=match(FamIDs,CovIDs)
Cov_keep <-Covar[IDs,]

temp = 0 
for (i in 1:dim(Cov_keep)[1]){
  if (as.character(Cov_keep[i,2])==as.character(Fam[i,2]))
    temp = temp + 1
}
if (temp == dim(Cov_keep)[1]){
  print ("Covariates file and FAM file share subjects order!")
}
if (temp < dim(Cov_keep)[1]){
  print ("Covariates file and FAM file have different subjects order!")
}  


# Construct an exclusion list and contain observations to be excluded.
Status = Cov_keep$er

# Accepted list of case and control group numberings
Accepted_Notation = c(0,1)

# Construct an exclusion list and contain observations to be excluded.
Exclusion_List = c()
for (i in 1:length(Status)){
  if (!Status[i] %in% Accepted_Notation){
    Exclusion_List = c(Exclusion_List,i)
  }
}

# Subset Cov file
Cov = Cov_keep[-Exclusion_List,]

# Read in all of the covariates
Y = Cov$er
Region <- as.factor(Cov$region_num)
PC1 = Cov$EV1
PC2 = Cov$EV2
PC3 = Cov$EV3
PC4 = Cov$EV4
PC5 = Cov$EV5
PC6 = Cov$EV6
PC7 = Cov$EV7
PC8 = Cov$EV8
PC9 = Cov$EV9
PC10 = Cov$EV10


# Number of lines to process simultaneously
N_lines = 1000
#Dimension = as.numeric(read.table(paste(Dimension_Directory,filename,".Dimension.txt", sep=""))[1])
Dimension=38
D = ceiling(Dimension/N_lines)


# fit a single null model for all of these SNPs
reg_noage.null = glm (Y ~ factor(Region) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, family=binomial)

for (d in 1:D){
  Dosage_Raw=read.table(paste(Dosage_Directory,filename,sep=""),header=FALSE,sep="",nrows=N_lines,skip=N_lines*(d-1))
 
  #Which subjects are to be kept -- look at the overlap of cov and fam files
  include <- which(Fam[,2] %in% Cov$IID)
  keep_columns <- matrix(rbind(6+3*(include-1), 6+3*(include-1)+1, 6+3*(include-1)+2),ncol=1,byrow=F)
  
  #Keep only the columns associated with samples to keep + 5 intro columns of dosage file
  dat.2 <- Dosage_Raw[,c(1:5,keep_columns)] 
  
  rm(Dosage_Raw) #Remove the data to make room and make processing faster
  
  # Genearate dosage variable
  Dosage = matrix(nrow=dim(Cov)[1],ncol=1)
  apply (dat.2,1,FUN=function(L){
    Chr = L[1]
    MarkerName = L[2]
    Position = L[3]
    A1 = L[4]
    A2 = L[5]
     
    for (i in 1:dim(Cov)[1]){
      Dosage[i,1]= 2*as.numeric(L[6+3*(i-1)])+as.numeric(L[6+3*(i-1)+1])
    }
    
    # Determine if this SNP has MAF<0.01
    Frequency = mean(Dosage)/2
    if (Frequency<0.01 | Frequency>0.99){
      Rare = 1
    }
    if (Frequency>=0.01 && Frequency<=0.99){
      Rare = 0
    }
    
    
    # Logistic regression 
    # do one step of Fisher's scoring  using maxit=1  
    reg_noage = glm(Y ~ Dosage + Region + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=Cov,
                    maxit=1,start=c(reg_noage.null$coefficients,0), family=binomial)
    
 
    if (!"Dosage" %in% rownames(summary(reg_noage)$coefficients)){
      Estimate_noage = NA
      SE_noage = NA
      P_noage.Wald = NA
    }
    
    if ("Dosage" %in% rownames(summary(reg_noage)$coefficients)){
      Estimate_noage = coef(summary(reg_noage))["Dosage","Estimate"]
      SE_noage = coef(summary(reg_noage))["Dosage","Std. Error"]
      P_noage.Wald = coef(summary(reg_noage))["Dosage","Pr(>|z|)"]
    }
    
    Score_noage.Z=Estimate_noage/SE_noage
    P_noage.Score=P_noage.Wald
    
    Score_noage.Z.1=glm.scoretest(reg_noage.null,Dosage,dispersion=NULL)
    P_noage.Score.1=2*(1-pnorm(abs(Score_noage.Z.1)))
    

    Sample <- length(Y)
    N.cases <- sum(Y)
    N.controls <- Sample - N.cases
    
    Output = cbind(MarkerName,Chr,Position,A1,A2,N.cases,N.controls,Frequency,Estimate_noage,SE_noage,Score_noage.Z,P_noage.Score,
                   Score_noage.Z.1,P_noage.Score.1,Rare)
    write.table(Output, paste(Output_Directory,filename,".Associationd.txt", sep=""), append=T, quote=F, sep="\t", row.names=F, col.names=F)    
    
  } ) }



#####################################################################################################
# I MADE CHANGES --- SEARCH FOR  "annotate2" & "SNPlist2" and "annotate3" & "SNPlist3"
#####################################################################################################

# Stephen Turner
# http:/StephenTurner.us/
# http:/GettingGeneticsDone.blogspot.com/
# See license at http:/gettinggeneticsdone.blogspot.com/p/copyright.html


# Last updated: Tuesday, April19, 2011
# R code for making manhattan plots and QQ plots from plink output files.
# manhattan() with GWAS data this can take a lot of memory, recommended for use on 64bit machines only, for now.
# Altnernatively, use bmanhattan() , i.e., base manhattan. uses base graphics. way faster.


## This is for testing purposes.
# set.seed(42)
# nchr=23
# nsnps=1000
# d=data.frame(
# SNP=sapply(1:(nchr*nsnps), function(x) paste("rs",x,sep='')),
# CHR=rep(1:nchr,each=nsnps),
# BP=rep(1:nsnps,nchr),
# P=runif(nchr*nsnps)
# )
# annotatesnps <- d$SNP[7550:7750]

# manhattan plot using base graphics

# manhattan plot using base graphics
manhattan = function(dataframe, colors=c("gray10", "gray50"), ymax="max", cex.x.axis=1, limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, annotate2=NULL, annotate3=NULL, ...) {
  
  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
  d$logp = -log10(d$P)
  d$pos=NA
  ticks=NULL
  lastbase=0
  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  if (ymax=="max") ymax<-ceiling(max(d$logp))
  if (ymax<8) ymax<-8
  
  numchroms=length(unique(d$CHR))
  if (numchroms==1) {
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
  } else {
    for (i in unique(d$CHR)) {
      if (i==1) {
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
      }
      ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    }
  }
  
  if (numchroms==1) {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  } else {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(p))), xlab="Chromosome", xaxt="n", type="n", ...))
    axis(1, at=ticks, lab=unique(d$CHR), cex.axis=cex.x.axis)
    icol=1
    for (i in unique(d$CHR)) {
      with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
      icol=icol+1
    }
  }
  
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, logp, col="green", ...))
  }
  
  if (!is.null(annotate2)) {
    d.annotate2=d[which(d$SNP %in% annotate2), ]
    with(d.annotate2, points(pos, logp, col="red", ...))
  }
  
  if (!is.null(annotate3)) {
    d.annotate3=d[which(d$SNP %in% annotate3), ]
    with(d.annotate3, points(pos, logp, col="darkolivegreen", ...))
  }
  
  if (suggestiveline) abline(h=suggestiveline, col="red")
  if (genomewideline) abline(h=genomewideline, col="azure4")
}

# Base graphics qq plot
qq = function(pvector, ...) {
  if (!is.numeric(pvector)) stop("D'oh! P value vector is not numeric.")
  pvector <- pvector[!is.na(pvector) & pvector<1 & pvector>0]
  o = -log10(sort(pvector,decreasing=F))
  #e = -log10( 1:length(o)/length(o) )
  e = -log10( ppoints(length(pvector) ))
  plot(e,o,pch=19,cex=1, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
  abline(0,1,col="red")
}


### OLD GGPLOT2 CODE ###

# manhattan plot using ggplot2
gg.manhattan = function(dataframe, title=NULL, max.y="max", suggestiveline=0, genomewideline=-log10(5e-8), size.x.labels=9, size.y.labels=10, annotate=F, SNPlist=NULL, annotate2=F, SNPlist2=NULL, annotate3=F, SNPlist3=NULL) {
  library(ggplot2)
  if (annotate & is.null(SNPlist)) stop("You requested annotation but provided no SNPlist!")
  d=dataframe
  #limit to only chrs 1-23?
  d=d[d$CHR %in% 1:23, ]
  if ("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d) ) {
    d=na.omit(d)
    d=d[d$P>0 & d$P<=1, ]
    d$logp = -log10(d$P)
    d$pos=NA
    ticks=NULL
    lastbase=0
    #new 2010-05-10
    numchroms=length(unique(d$CHR))
    if (numchroms==1) {
      d$pos=d$BP
    } else {
      
      for (i in unique(d$CHR)) {
        if (i==1) {
          d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
        } else {
          lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
          d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
        }
        ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
      }
      ticklim=c(min(d$pos),max(d$pos))
      
    }
    mycols=rep(c("gray10","gray60"),max(d$CHR))
    if (max.y=="max") maxy=ceiling(max(d$logp)) else maxy=max.y
    if (maxy<10) maxy=10
    if (annotate) d.annotate=d[as.numeric(substr(d$SNP,3,100)) %in% SNPlist, ]
    if (annotate2) d.annotate2=d[as.numeric(substr(d$SNP,3,100)) %in% SNPlist2, ]
    if (annotate3) d.annotate3=d[as.numeric(substr(d$SNP,3,100)) %in% SNPlist3, ]
    if (numchroms==1) {
      plot=qplot(pos,logp,data=d,ylab=expression(-log[10](italic(p))), xlab=paste("Chromosome",unique(d$CHR),"position"))
    } else {
      plot=qplot(pos,logp,data=d, ylab=expression(-log[10](italic(p))) , colour=factor(CHR))
      plot=plot+scale_x_continuous(name="Chromosome", breaks=ticks, labels=(unique(d$CHR)))
      plot=plot+scale_y_continuous(limits=c(0,maxy), breaks=1:maxy, labels=1:maxy)
      plot=plot+scale_colour_manual(value=mycols)
    }
    if (annotate) plot=plot + geom_point(data=d.annotate, colour=I("blue"))
    plot=plot + opts(legend.position = "none")
    plot=plot + opts(title=title)
    plot=plot+opts(
      panel.background=theme_blank(),
      panel.grid.minor=theme_blank(),
      axis.text.x=theme_text(size=size.x.labels, colour="grey50"),
      axis.text.y=theme_text(size=size.y.labels, colour="grey50"),
      axis.ticks=theme_segment(colour=NA)
    )
    
    if (annotate2) plot=plot + geom_point(data=d.annotate2, colour=I("red"))
    plot=plot + opts(legend.position = "none")
    plot=plot + opts(title=title)
    plot=plot+opts(
      panel.background=theme_blank(),
      panel.grid.minor=theme_blank(),
      axis.text.x=theme_text(size=size.x.labels, colour="grey50"),
      axis.text.y=theme_text(size=size.y.labels, colour="grey50"),
      axis.ticks=theme_segment(colour=NA)
    )
    
    if (annotate3) plot=plot + geom_point(data=d.annotate3, colour=I("darkolivegreen"))
    plot=plot + opts(legend.position = "none")
    plot=plot + opts(title=title)
    plot=plot+opts(
      panel.background=theme_blank(),
      panel.grid.minor=theme_blank(),
      axis.text.x=theme_text(size=size.x.labels, colour="grey50"),
      axis.text.y=theme_text(size=size.y.labels, colour="grey50"),
      axis.ticks=theme_segment(colour=NA)
    )
    
    
    if (suggestiveline) plot=plot+geom_hline(yintercept=suggestiveline,colour="red", alpha=I(1/3))
    if (genomewideline) plot=plot+geom_hline(yintercept=genomewideline,colour="azure4")
    plot
  } else {
    stop("Make sure your data frame contains columns CHR, BP, and P")
  }
}

gg.qq = function(pvector, title=NULL, spartan=F) {
  library(ggplot2)
  o = -log10(sort(pvector,decreasing=F))
  #e = -log10( 1:length(o)/length(o) )
  e = -log10( ppoints(length(pvector) ))
  plot=qplot(e,o, xlim=c(0,max(e)), ylim=c(0,max(o))) + stat_abline(intercept=0,slope=1, col="black")
  plot=plot+opts(title=title)
  plot=plot+scale_x_continuous(name=expression(Expected~~-log[10](italic(p))))
  plot=plot+scale_y_continuous(name=expression(Observed~~-log[10](italic(p))))
  if (spartan) plot=plot+opts(panel.background=theme_rect(col="grey50"), panel.grid.minor=theme_blank())
  plot
}

gg.qqman = function(data="plinkresults") {
  myqqplot = ggqq(data$P)
  mymanplot = ggmanhattan(data)
  ggsave(file="qqplot.png",myqqplot,w=5,h=5,dpi=100)
  ggsave(file="manhattan.png",mymanplot,width=12,height=9,dpi=100)
}

gg.qqmanall= function(command="ls *assoc") {
  filelist=system(command,intern=T)
  datalist=NULL
  for (i in filelist) {datalist[[i]]=read.table(i,T)}
  highestneglogp=ceiling(max(sapply(datalist, function(df) max(na.omit(-log10(df$P))))))
  print(paste("Highest -log10(P) = ",highestneglogp),quote=F)
  start=Sys.time()
  for (i in names(datalist)) {
    myqqplot=ggqq(datalist[[i]]$P, title=i)
    ggsave(file=paste("qqplot-", i, ".png", sep=""),myqqplot, width=5, height=5,dpi=100)
    mymanplot=ggmanhattan(datalist[[i]], title=i, max.y=highestneglogp)
    ggsave(file=paste("manhattan-", i, ".png", sep=""),mymanplot,width=12,height=9,dpi=100)
  }
  end=Sys.time()
  print(elapsed<-end-start)
}


####################################################################
####################################################################
####################################################################

## SET WORKING DIRECTORY
setwd("F:/Project/R_Results/Imputed/ER")
getwd()



####################################################################
## READ IN THE DATA


red<-read.table("Overall_ERpos_ERneg_ALL.metrics.txt", fill=TRUE, header=TRUE)
red[1:10, ]
dim(red)

names(red)<-c("SNP","A1","A2", "Chr","BP", "N_Cases", "N_Controls", "Frequency", "Estimate", "SE", "P_Score", "Rare",
              "chrom", "SNP2", "allele1", "allele2", "Est_erneg", "SE_erneg", "score_Erneg", "P_erneg", "position", "a0",
              "freq_a1", "info", "type", "Est_erpos", "SE_erpos", "score_Erpos", "P_erpos" )

one <- na.omit(red)
one[1:10, ]
dim(one)

P = one$P_erneg
BP = one$BP
CHR = one$Chr

dat2= cbind.data.frame(CHR,BP,P)

tiff(filename = "F:/Project/R_Results/P_Score_ERneg_ALL_MH.tiff", type = "cairo")
manhattan(dat2, colors=c("midnightblue","gray59","red3"), pch=20, genomewideline=-log10(1e-9), suggestiveline=F)
title(main="Manhattan-plot Pvalue_Score ER Neg - ALL")
dev.off()

tiff(filename = "F:/Project/R_Results/P_Score_ERneg_ALL_QQ.tiff", type = "cairo")
qq(one$P_erneg)
title(main="QQ-plot Pvalue_Score ER Neg- ALL")
dev.off()

######## calculate the lambda
chisq_erneg <- qchisq (1-one$P_erneg, 1)

#lambda erneg= 1.042257
med_erneg <- median (chisq_erneg / qchisq(0.5,1))
med_erneg

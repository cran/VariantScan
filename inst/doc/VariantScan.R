## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- install-----------------------------------------------------------------
#Install from CRAN
#install.packages("VariantScan")
## or you can get the latest version of HierDpart from github
#library(devtools)

#install_github("xinghuq/VariantScan")
library("ggplot2")
library("VariantScan")


## ---- library-----------------------------------------------------------------
# example genepop file
f <- system.file('extdata',package='VariantScan')
infile <- file.path(f, "sim1.csv")
## read genotype file
geno=read.csv(infile)

# traits
traitq=geno[,14]
genotype=geno[,-c(1:14)]

# get PCs as covariates

PCs=prcomp(genotype)
summary(PCs$x[,1:2])


## -----------------------------------------------------------------------------
loessW=VScan(x=genotype,y=(traitq),methods ="loess")


## -----------------------------------------------------------------------------

loessWcv=VScan(x=genotype,y=(traitq),U=PCs$x[,1:2],methods ="loess")


## -----------------------------------------------------------------------------
lmW=VScan(x=genotype,y=(traitq),methods ="lm")

lmWcv=VScan(x=genotype,y=(traitq),U=PCs$x[,1:2],methods ="lm")


## ----ggplot2, fig1, fig.height = 5, fig.width = 10.5, fig.align = "center"----
## 
Loci<-rep("Neutral", 1000)
Loci[c(201,211,221,231,241,251,261,271,281,291)]<-"QTL"
Selected_Loci<-Loci[-which(Loci=="Neutral")]


g1=ggplot() +
  geom_point(aes(x=which(Loci=="Neutral"), y=-log10(lmWcv$p_norm$p.value[-which(Loci!="Neutral")])), col = "gray83") +
  geom_point(aes(x=which(Loci!="Neutral"), y=-log10(lmWcv$p_norm$p.value[-which(Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNPs") + ylab("-log10(p-value)") +ylim(c(0,35))+theme_bw()

g1




## ----ggpot2, fig2, fig.height = 5, fig.width = 10.5, fig.align = "center"-----

## Manhattan plot


g2=ggplot() +
  geom_point(aes(x=which(Loci=="Neutral"), y=-log10(loessWcv$p_norm$p.value[-which(Loci!="Neutral")])), col = "gray83") +
  geom_point(aes(x=which(Loci!="Neutral"), y=-log10(loessWcv$p_norm$p.value[-which(Loci=="Neutral")]), colour = Selected_Loci)) +
  xlab("SNPs") + ylab("-log10(p-value)") +ylim(c(0,35))+theme_bw()

g2



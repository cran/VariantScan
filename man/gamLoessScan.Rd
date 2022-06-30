\name{gamLoessScan}
\alias{gamLoessScan}
\title{Variants (Biomarkers, e.g., genomic (genetic loci), transcriptomic (gene expression), epigenomic (methylations), proteomic(protein), metabolomic (metabolites) variants) Scanning and Association Tests Using Generalized Additive Model with  Local Polynomial Regression (LOESS).

}
\description{Fitting a Generalized Additive Mixed Models (GAMM) with Local Polynomial Regression in association testing.

}
\usage{
gamLoessScan(genotype, traits, U, cv_method = "adaptive_cv", 
model_metric = "RMSE", n_hyperparameter_search = 10,verbose=TRUE, ...)
}

\arguments{
   \item{genotype}{Varants/genotypes matrix coding in reference allele (0,1,2) or variant count

}
  \item{traits}{Traits

}
  \item{U}{ Covariates or confounding factors

}
  \item{cv_method}{Cross-validation

}
  \item{model_metric}{Model performance metrics, based on which the optimal model is determined.

}
  \item{n_hyperparameter_search}{Number of hyperparameters for tuning

}
 \item{verbose}{whether print training messages.

}
  \item{\dots}{other arguments passing to generalized additive mixed models (gam)

}
 
}
\details{Fits the specified generalized additive mixed model (GAMM) with LOESS smoothness.

}
\value{The weights of variants as well as their p-values

}
\references{
Wood S.N. (2006b) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press.

Wang, Y. (1998) Mixed effects smoothing spline analysis of variance. J.R. Statist. Soc. B 60, 159-174.

Lin, X and Zhang, D. (1999) Inference in generalized additive mixed models by using smoothing splines. JRSSB. 55(2):381-400.

}

\examples{
# not run
\donttest{f <- system.file('extdata',package='VariantScan')
infile <- file.path(f, "sim1.csv")
geno=read.csv(infile)

traitq=geno[,14]
genotype=geno[,-c(1:14)]
PCs=prcomp(genotype)
test=gamLoessScan(genotype =genotype,traits =(traitq),U=PCs$x[,1:2],n_hyperparameter_search=5)
}
}

\keyword{gamLoessScan}
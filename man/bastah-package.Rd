\name{bastah}
\alias{bastah}
\docType{package}
\title{
\packageTitle{bastah}
}
\description{
\packageDescription{bastah}
}

\usage{
bastah (X, y, family = "gaussian", mcorr = "holm", N = 10000, ncores = 4, verbose = FALSE)
}
\arguments{
  \item{X}{
    An n by p numeric design matrix with p columns for p predictor variables and n rows corresponding to n observations.
  }
  \item{y}{
    A numeric response variable of length n.
  }
  \item{family}{
    Family of the response variable. It should be either "gaussian" or "binomial". (default = "gaussian")
  }
  \item{mcorr}{
    Multiple correction method. It can be either "WY" or any of p.adjust.methods. (default = "holm")
  }
  \item{N}{
    It is the number of samples to take for the empirical distribution which is used to correct the pvalues if multiple correction method is "WY" (Westfall-Young). (default = 10000)
  }
  \item{ncores}{
    Maximum number of cores to be used for parallel execution. (default = 4)
  }
  \item{verbose}{
    Prints more information if this is set to TRUE. (default = FALSE)
  }
}
\details{
In this package lasso.proj function of hdi package is updated for application on big data. The original lasso.proj is updated by replacing node-wise regression with scaled lasso. BigQUIC is used for sparse precision matrix calculation. Data is always normalized before processing. Normalization technique used by Vlaming and Groenen (2014) is used.
The method has been successfully used on large SNP (Single Nucleotide Polymorphism) datasets for GWAS (Genomewide Association Study).

The package can use scikit-learn (http://scikit-learn.org) for a better performance. It is advised to install doMC, rPython, python, numpy and scikit-learn. The package uses scikit-learn at runtime, therefore, python, numpy and scikit-learn are not required for package installation and can be installed after installation of the package.

NOTE: We have noticed that lars package in R crashes, so it is recommended to use scikit-learn.

NOTE: In preprocessing step, variables having a constant value are not considered. The list of varaibles used is returned in selection varaible of the result.
}

\value{
An object with Class \code{"bastah"}
  \item{pval}{Calculated p-values}
  \item{pval.corr}{Corrected p-values}
  \item{sigmahat}{Estimated standard deviation}
  \item{bhat}{Estimated coefficients}
  \item{selection}{Indicies of variables selected for analysis}
}

\author{
\packageAuthor{bastah}

Maintainer: \packageMaintainer{bastah}
}
\references{
C. Hsieh, M. Sustik, I. Dhillon, P. Ravikumar, R. Poldrack. In Neural Information Processing Systems (NIPS), December 2013. (Oral)

S. van de Geer, P. Buhlmann, Y. Ritov and R. Dezeure (2014) On asymptotically optimal confidence regions and tests for high-dimensional models. Annals of Statistics
42, 1166-1202.

C. Zhang, S. Zhang(2014) Confidence intervals for low dimensional parameters in high dimensional linear models. Journal of the Royal Statistical Society: Series B76, 217-242.

P. Buhlmann and S. van de Geer(2015) High-dimensional inference in misspecified linear models. Electronic Journal of Statistics 9, 1449-1473.

R. de Vlaming and P. J. F. Groenen (2014) The current and future use of ridge regression for prediction in quantitative genetics. BioMed Research International, 2015, 143712.
}
\keyword{BigQUIC}
\keyword{package}
\keyword{lasso}
\keyword{projection}
\keyword{p-values}
\keyword{lasso.proj}
\keyword{GWAS}

\examples{
# The package is accompanied with a simulated genome-wide association
# study dataset "snps" containing n=100 observations of p=500 predictors
   data(snps)
# The association of SNPs to the phenotype can be identified using bastah
# NOTE: We have noticed that lars package in R crashes,
# so it is recommended to use scikit-learn (see package details).
\dontrun{
      result = bastah(X = snps$X, y = snps$y, family = "binomial", verbose = TRUE)
   }
}

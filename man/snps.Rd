\name{snps}
\alias{snps}
\docType{data}
\title{SNP data set}
\description{Simulated Single Polymorphism Nucleotide (SNP) dataset
  containing \eqn{n=100} observations of \eqn{p=500} predictors (SNPs, 1=Homozygote1, 0=Heterozygote, -1=Homozygote2) and a one-dimensional response (1=case, 0=control). The dataset is generated using GWAsimulator.}
\usage{data(snps)}
\format{
  \describe{
    \item{y}{Phenotype (1=case, 0=control) of 100 individuals.}
    \item{x}{SNP genotype data (SNPs, 1=Homozygote1, 0=Heterozygote, -1=Homozygote2) of 500 simulated SNPs.}
  }
}
\references{
C. Li and M. Li (2008) GWAsimulator: a rapid whole-genome simulation program. Bioinformatics 24 (1): 140-142.
}
\examples{data(snps)}
\keyword{datasets}

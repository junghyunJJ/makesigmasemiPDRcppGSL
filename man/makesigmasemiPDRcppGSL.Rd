\name{makesigmasemiPDRcppGSL}
\alias{makesigmasemiPDRcppGSL}
\title{positive semidefinite}
\description{
Compute the positive semidefinite
}
\usage{
makesigmasemiPDRcppGSL(geno)
}
\arguments{
  \item{geno}{a m x n numeric matrix}
}
\details{
Compute the positive semidefinite using RcppGSL
}
\references{
  Leveraging allelic heterogeneity to increase power of association testing
(https://www.biorxiv.org/content/10.1101/498360v1)
}
\examples{
devtools::install_github("junghyunJJ/marsR")
library(marsR)
data(testdata)
makesigmasemiPDRcppGSL(geno)
}
\author{Junghyun Jung}



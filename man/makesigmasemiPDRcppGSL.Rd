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
data("testdata")
makesigmasemiPDRcppGSL(testdata$geno)[1:5,1:5]
}
\author{Junghyun Jung}



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimateBACsizeFromVDV.R
\name{estimateBACsizeFromVDV}
\alias{estimateBACsizeFromVDV}
\title{Estimate the size of a BAC using the length of all VDV reads}
\usage{
estimateBACsizeFromVDV(vdvLength, nclust = 2:10, method = c("jenks", "pam"))
}
\arguments{
\item{vdvLength}{integer vector of lengths of VDV reads}

\item{nclust}{integer vector of numbers of clusters to use (default to 2:10)}

\item{method}{character string indicating which method to use ("jenks" or "pam"). Default is "jenks" (faster)}
}
\value{
tibble
}
\description{
Cluster the sizes of VDV read for different number of clusters
(use \code{LongestDNA} rather than \code{ReadLength}, see details).
For each clustering, obtain the median size for each cluster and take the largest one.
The estimated BAC size is the median of these larger medians
}
\details{
In order to use the result of this function with the \code{\link{FilterBACreads}} function
          you should input the "LongestDNA" values instead of the ReadLength values for VDV reads.
}
\examples{
## Generate some random insert sizes
set.seed(12345)
InsertLengths <- c(sample(102e3:106e3, 30),
                   sample(43e3:45e3, 10),
                   sample(23e3:25e3, 10))
## Estimate the BAC size (smaller lengths are considered likely recombinants)
estimateBACsizeFromVDV(InsertLengths, method = "jenks")$BACsize

}
\author{
Pascal GP Martin
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{pm_fluidity_all}
\alias{pm_fluidity_all}
\title{Genome Fluidity}
\usage{
pm_fluidity_all(Panmatrix)
}
\arguments{
\item{Panmatrix}{Panmatrix produced by make_panmatrix functions}
}
\description{
This function computes fluidity without sampling.
}
\details{
Fluidity takes values in [0,1], with 1 denoting no common genes.This metric was introduced is called the Sorensen distance (Deza & Deza, 2009) and was first introduced in the
  context of a pangenome analysis in (Kislyuk et al ,2011).
}
\examples{
\dontrun{pm_fluidity_all(Panmatrix)}

}
\references{
A. O. Kislyuk, B. Haegeman, N. H. Bergman, and J. S. Weitz, "Genomic fluidity:
 an integrative view of gene diversity within microbial populations," BMC genomics, pp. 12-32, 2011.
M. M. Deza and E. Deza, Encyclopedia of Distances. Springer, 2009.
}

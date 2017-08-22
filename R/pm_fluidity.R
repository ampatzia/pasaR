#'  Genome Fluidity with sampling
#'
#' This function computes fluidity with sampling, as implemented on package micropan (optimized for speed).
#'   Fluidity takes values in [0,1], with 1 denoting no common genes.
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param n.sim Number of simulations,defaults to 10
#' @details Fluidity takes values in [0,1], with 1 denoting no common genes.This metric was introduced is called the Sorensen distance (Deza & Deza, 2009) and was first introduced in the
#' context of a pangenome analysis in (Kislyuk et al ,2011).
#' @export
#' @examples \dontrun{pm_fluidity(Panmatrix, n.sim=100)}
#'
#' @references
#' A. O. Kislyuk, B. Haegeman, N. H. Bergman, and J. S. Weitz, "Genomic fluidity???: an integrative view of gene diversity within microbial populations," BMC genomics, pp. 12-32, 2011.
#' M. M. Deza and E. Deza, Encyclopedia of Distances. Springer, 2009.

pm_fluidity <- function (Panmatrix, n.sim = 10)
{
  panm <- sapply(as.data.frame(Panmatrix),function(x) as.logical(x))
  ng <- dim(panm)[1]
  flu <- rep(0, n.sim)
  for (i in 1:n.sim) {
    ii <- sample(ng, 2)
    g1 <- panm[ii[1], ]
    g2 <- panm[ii[2], ]
    flu[i] <- sum(abs(g1-g2))/(sum(g1) + sum(g2))
  }
  flu.list <- list(Mean = mean(flu), Std = stats::sd(flu))
  return(flu.list)
}

#'  Binomial mixture fitting
#'
#' This function determines the closedness of the pangenome using a binomial mixture fit of the data.
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param K.range Range of model components to be tested, defaults to search from 2 to 8
#' @param core.detect.prob Detection probability of core genes. Defaults to 1.0
#'
#' @details The function returns pangenome and core size estimation along with the propabilities
#'  of the components.This is an optimized version of the binomixEstimate() function from package
#'  micropan  (Snipen & Hiland, 2015).
#'
#' @export
#' @examples \dontrun{pm_binom(panm, K.range =2:8)}
#'
#' @references
#' L. Snipen and K. H. Liland, "micropan: an R-package for microbial pan-genomics.," BMC bioinformatics, vol. 16, p. 79, 2015

pm_binom <- function (Panmatrix, K.range = 2:8, core.detect.prob = 1)
{
  pan.matrix <- sapply(Panmatrix, function(x) as.logical(x))
  y <- table(factor(colSums(pan.matrix), levels = 1:dim(pan.matrix)[1]))
  bic.tab <- matrix(NA, nrow = length(K.range), ncol = 3)
  colnames(bic.tab) <- c("Core.size", "Pan.size", "BIC")
  rownames(bic.tab) <- paste(K.range, "components")
  mix.list <- vector("list", length(K.range))
  for (i in 1:length(K.range)) {
    lst <- binomixMachine(y, K.range[i], core.detect.prob)
    bic.tab[i, ] <- lst[[1]]
    mix.list[[i]] <- lst[[2]]
  }
  if (bic.tab[length(K.range), 3] == min(bic.tab[, 3]))
    warning("Minimum BIC at maximum K, increase upper limit of K.range")
  binomix <- list(BIC.table = bic.tab, Mix.list = mix.list)

  return(binomix)
}

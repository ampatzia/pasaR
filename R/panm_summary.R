#'  Panmatrix Summary
#'
#' This function produces a frequency table of Genome participation in clusters
#' @param Panmatrix  Panmatrix produced by make_panmatrix functions
#' @export
#' @examples \dontrun{panm_summary(Panmatrix)}
#'

panm_summary <- function (Panmatrix)
{
  object <- sapply(Panmatrix, function(x) as.logical(x))
  levs <- 1:nrow(object)
  y <- as.data.frame(table(factor(colSums(object), levels = levs)))
  colnames(y) <- c("Genomes", "Clusters")
  return(y)
}

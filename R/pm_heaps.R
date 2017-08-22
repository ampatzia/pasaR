#'  Heap's Law fitting
#'
#' This function determines the closedness of the pangenome using a curve fit of the data
#'    based on Heaps Law.
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param n_perm Number of permutations
#'
#' @details The regression fit returns two estimated parameters, intercept and decay parameter a.  If a<1 then the pangenome is considered to be open. This is an optimized version
#'    of the heaps() function from package micropan. The theoretical aspects are discussed in the canonical work of Tettelin et  al. (2008)
#'
#' @export
#' @examples  \dontrun{pm_heaps(Panmatrix,n_perm=100)}
#'
#' @references Tettelin, H., Riley, D., Cattuto, C., Medini, D. (2008). Comparative genomics: the bacterial pan-genome. Current Opinions in Microbiology, 12:472-477.

pm_heaps <- function (Panmatrix, n_perm){
  if (missing(n_perm)) {n_perm = 100}
  genomes<-NULL
  genes<-NULL
  pan.matrix <- sapply(Panmatrix, function(x) as.logical(x))
  ng <- nrow(Panmatrix)
  nmat <- matrix(0, nrow = (ng - 1), ncol = n_perm)

  nmat<-replicate(n_perm,{

    cm <- apply(pan.matrix[sample(ng), ], 2, cumsum)
    rowSums((cm == 1)[2:ng, ] & (cm == 0)[1:(ng -  1), ])
  })

  nmat<-t(nmat)
  colnames(nmat) <- c(2:(ncol(nmat)+1))
  nmat <- tidyr::gather(as.data.frame(nmat), genomes, genes)%>%transform(., genomes=as.numeric(genomes))


  p0 <- c(mean(nmat$genes[nmat$genomes == 2]), 1)


  objectFun<-function (p, x, y)
  {
    y.hat <- p[1] * x^(-p[2])
    J <- sqrt(sum((y - y.hat)^2))/length(x)
    return(J)
  }

  fit <- stats::optim(p0, objectFun, gr = NULL, nmat$genomes, nmat$genes, method = "L-BFGS-B",
                      lower = c(0, 0), upper = c(10000, 2))
  p.hat <- fit$par
  names(p.hat) <- c("Intercept", "alpha")
  return(p.hat)
}

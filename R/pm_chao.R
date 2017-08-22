#'  Chao lower bound estimator
#'
#' This function computes two versions of the lower bound Chao Estimator for the pangenome
#'    size.
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param biased produces a biased estimation, same as the chao() function in micropan
#'
#' @details The function prodused a biased & unbiased version of the metric along with variance
#'    and a 95% CI . This is a conservative estimator, ie. it is likely to produce a small
#'    estimation of the pangenome compared to other estimators.
#'    This function is an optimized and expanded version of chao() from package micropan.
#'
#'
#' @export
#' @examples \dontrun{pm_chao(Panmatrix)}
#'
#' @references
#' Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43, 783-791.
#' Gotelli, N.J. and Colwell, R.K., 2011. Estimating species richness. Biological diversity: frontiers in measurement and assessment, 12, pp.39-54.

pm_chao <- function (Panmatrix,biased=FALSE) {
  panm <- sapply(as.data.frame(Panmatrix),function(x) as.logical(x))
  y <- table(factor(colSums(panm), levels = 1:dim(panm)[1]))

  if (y[2] != 0){
    pan.size.biased <- round(sum(y) + y[1]^2/(2 * y[2]))
  }


  pan.size <- round(sum(y) + (y[1]*(y[1]-1))/(2 * (y[2]+1))) #is bias corrected
  fratio <- y[1]/y[2]
  chao.variance <- y[2]*(0.5*(fratio)^2 + fratio^3 + 0.25*fratio^4 ) #chao 1987

  C_var <- exp(1.96*sqrt(log(1+chao.variance/(pan.size-sum(y)))))

  l_bound <- sum(y) +(pan.size-sum(y))/C_var
  u_bound <- sum(y) +(pan.size-sum(y))*C_var

  if(is.null(pan.size.biased)) {pan.size.biased<-NA}
  if( biased==TRUE){
    results <- c(pan.size.biased, pan.size, chao.variance, l_bound, u_bound)

    names(results) <- c("Estimated pangenome size - Biased","Estimated pangenome size", " Estimator Variance","CI (95%) -Lower Bound", "CI (95%)- Upper Bound")
  }else{
    results <- c(pan.size,chao.variance,l_bound,u_bound)

    names(results) <- c("Estimated pangenome size", " Estimator Variance","CI (95%) -Lower Bound", "CI (95%)- Upper Bound")
  }

  return(results)

}

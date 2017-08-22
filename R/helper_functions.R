#Helper Functions
# 1) gtools_comb
# 2) binomixMachine
# 3) negTruncLogLike



#'  Helper function: gtools::combination
#'
#' This helper function is produces all possible combinations for a set of numbers.
#'   This is an direct copy of the combination() function from package gtools, and is
#'   used in function pm_fludity_all().
#'
#' @keywords internal




gtools_comb <- function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE)
{
  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) !=
      0)
    stop("bad value of n")
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) !=
      0)
    stop("bad value of r")
  if (!is.atomic(v) || length(v) < n)
    stop("v is either non-atomic or too short")
  if ((r > n) & repeats.allowed == FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n)
      stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  if (repeats.allowed)
    sub <- function(n, r, v) {
      if (r == 0)
        v0
      else if (r == 1)
        matrix(v, n, 1)
      else if (n == 1)
        matrix(v, 1, r)
      else rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n -
                                                            1, r, v[-1]))
    }
  else sub <- function(n, r, v) {
    if (r == 0)
      v0
    else if (r == 1)
      matrix(v, n, 1)
    else if (r == n)
      matrix(v, 1, n)
    else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])),
               Recall(n - 1, r, v[-1]))
  }
  sub(n, r, v[1:n])
}



#'  Binomix machine
#'
#' This function is a helper borrowed from package micropan to be used to compute
#'  binomial mixture pangenome models
#' @keywords internal


binomixMachine <- function (y, K, core.detect.prob = 1)
{
  n <- sum(y)
  G <- length(y)
  ctr <- list(maxit = 200*K, reltol = 1e-8)
  np <- K - 1
  pmix0 <- rep(1, np)/K
  pdet0 <- (1:np)/(np + 1)
  p.initial <- c(pmix0, pdet0)
  A <- rbind(c(rep(1, np), rep(0, np)), c(rep(-1, np), rep(0,
                                                           np)), diag(np + np), -1 * diag(np + np))
  b <- c(0, -1, rep(0, np + np), rep(-1, np + np))
  est <- stats::constrOptim(theta = p.initial, f = negTruncLogLike,
                            grad = NULL, method = "Nelder-Mead", control = ctr, ui = A,
                            ci = b, y = y, core.p = core.detect.prob)
  estimates <- numeric(3)
  names(estimates) <- c("Core.size", "Pan.size", "BIC")
  estimates[3] <- 2 * est$value + log(n) * (np + K -1)
  p.mix <- c(1 - sum(est$par[1:np]), est$par[1:np])
  p.det <- c(core.detect.prob, est$par[(np + 1):length(est$par)])
  ixx <- order(p.det)
  p.det <- p.det[ixx]
  p.mix <- p.mix[ixx]
  theta_0 <- choose(G, 0) * sum(p.mix * (1 - p.det)^G)
  y_0 <- n * theta_0/(1 - theta_0)
  estimates[2] <- n + round(y_0)
  ixx <- which(p.det >= core.detect.prob)
  estimates[1] <- round(estimates[2] * sum(p.mix[ixx]))
  mixmod <- matrix(c(p.det, p.mix), nrow = 2, byrow = T)
  rownames(mixmod) <- c("Detection.prob", "Mixing.prop")
  colnames(mixmod) <- paste("Comp_", 1:K, sep = "")
  return(list(estimates, mixmod))
}

#'   Truncated log likelihood
#'
#' This function is a helper borrowed from package micropan to be used to compute
#'  truncated log likelihood
#' @keywords internal

negTruncLogLike <- function (p, y, core.p)
{
  np <- length(p)/2
  p.det <- c(core.p, p[(np + 1):length(p)])
  p.mix <- c(1 - sum(p[1:np]), p[1:np])
  G <- length(y)
  K <- length(p.mix)
  n <- sum(y)
  theta_0 <- choose(G, 0) * sum(p.mix * (1 - p.det)^G)
  L <- -n * log(1 - theta_0)
  for (g in 1:G) {
    theta_g <- choose(G, g) * sum(p.mix * p.det^g * (1 -
                                                       p.det)^(G - g))
    L <- L + y[g] * log(theta_g)
  }
  return(-L)
}

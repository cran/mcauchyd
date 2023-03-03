kldcauchy <- function(sigma1, sigma2, eps = 1e-06) {
  #' Kullback-Leibler Divergence between Centered Multivariate Cauchy Distributions
  #'
  #' Computes the Kullback-Leibler divergence between two random vectors distributed
  #' according to multivariate Cauchy distributions (MCD) with zero location vector.
  #'
  #' @aliases kldcauchy
  #'
  #' @usage kldcauchy(sigma1, sigma2, eps = 1e-06)
  #' @param sigma1 symmetric, positive-definite matrix. The scatter matrix of the first distribution.
  #' @param sigma2 symmetric, positive-definite matrix. The scatter matrix of the second distribution.
  #' @param eps numeric. Precision for the computation. Default: 1e-06.
  #' @return A  numeric value: the Kullback-Leibler divergence between the two distributions,
  #' with two attributes \code{attr(, "epsilon")} (precision of the result) 
  #' and \code{attr(, "k")} (number of iterations).
  #'
  #' @details Given \eqn{X_1}, a random vector of \eqn{R^p} distributed according to the MCD
  #' with parameters \eqn{(0, \Sigma_1)}
  #' and \eqn{X_2}, a random vector of \eqn{R^p} distributed according to the MCD
  #' with parameters \eqn{(0, \Sigma_2)}.
  #' 
  #' Let \eqn{\lambda_1, \dots, \lambda_p} the eigenvalues of the square matrix \eqn{\Sigma_1 \Sigma_2^{-1}}
  #' sorted in increasing order: \deqn{\lambda_1 < \dots < \lambda_{p-1} < \lambda_p}
  #' Depending on the values of these eigenvalues,
  #' the computation of the Kullback-Leibler divergence of \eqn{X_1} from \eqn{X_2}
  #' is given by:
  #' \itemize{
  #' \item if \eqn{\lambda_1 < 1} and \eqn{\lambda_p > 1}:\cr
  #' \eqn{ \displaystyle{ KL(X_1||X_2) = -\frac{1}{2} \ln{ \prod_{i=1}^p{\lambda_i}} + \frac{1+p}{2} \bigg( \ln{\lambda_p} } } \cr
  #' \eqn{ \displaystyle{ - \frac{\partial}{\partial a} \bigg\{ F_D^{(p)} \bigg( a, \underbrace{\frac{1}{2}, \dots, \frac{1}{2}, a + \frac{1}{2}}_p ; a + \frac{1+p}{2} ; 1 - \frac{\lambda_1}{\lambda_p}, \dots, 1 - \frac{\lambda_{p-1}}{\lambda_p}, 1 - \frac{1}{\lambda_p} \bigg) \bigg\}\bigg|_{a=0} \bigg) } }
  #' 
  #' \item if \eqn{\lambda_p < 1}:\cr
  #' \eqn{ \displaystyle{ KL(X_1||X_2) = -\frac{1}{2} \ln{ \prod_{i=1}^p{\lambda_i}} } }
  #' \eqn{ \displaystyle{ - \frac{1+p}{2} \frac{\partial}{\partial a} \bigg\{ F_D^{(p)} \bigg( a, \underbrace{\frac{1}{2}, \dots, \frac{1}{2}}_p ; a + \frac{1+p}{2} ; 1 - \lambda_1, \dots, 1 - \lambda_p \bigg) \bigg\}\bigg|_{a=0} } }
  #' 
  #' \item if \eqn{\lambda_1 > 1}:\cr
  #' \eqn{ \displaystyle{ KL(X_1||X_2) = -\frac{1}{2} \ln{ \prod_{i=1}^p{\lambda_i}} - \frac{1+p}{2} \prod_{i=1}^p\frac{1}{\sqrt{\lambda_i}} } } \cr
  #' \eqn{ \displaystyle{ \times \frac{\partial}{\partial a} \bigg\{ F_D^{(p)} \bigg( \frac{1+p}{2}, \underbrace{\frac{1}{2}, \dots, \frac{1}{2}}_p ; a + \frac{1+p}{2} ; 1 - \frac{1}{\lambda_1}, \dots, 1 - \frac{1}{\lambda_p} \bigg) \bigg\}\bigg|_{a=0} } }
  #' }
  #' 
  #' where \eqn{F_D^{(p)}} is the Lauricella \eqn{D}-hypergeometric function defined for \eqn{p} variables:
  #' \deqn{ \displaystyle{ F_D^{(p)}\left(a; b_1, ..., b_p; g; x_1, ..., x_p\right) = \sum\limits_{m_1 \geq 0} ... \sum\limits_{m_p \geq 0}{ \frac{ (a)_{m_1+...+m_p}(b_1)_{m_1} ... (b_p)_{m_p} }{ (g)_{m_1+...+m_p} } \frac{x_1^{m_1}}{m_1!} ... \frac{x_p^{m_p}}{m_p!} } } }
  #' 
  #' The computation of the partial derivative uses the \code{\link{pochhammer}} function.
  #' 
  #' @author Pierre Santagostini, Nizar Bouhlel
  #' @references N. Bouhlel, D. Rousseau, A Generic Formula and Some Special Cases for the Kullback–Leibler Divergence between Central Multivariate Cauchy Distributions.
  #' Entropy, 24, 838, July 2022.
  #' \doi{10.3390/e24060838}
  #'
  #' @examples
  #' \donttest{
  #' sigma1 <- matrix(c(1, 0.6, 0.2, 0.6, 1, 0.3, 0.2, 0.3, 1), nrow = 3)
  #' sigma2 <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.4, 0.1, 0.4, 1), nrow = 3)
  #' kldcauchy(sigma1, sigma2)
  #' kldcauchy(sigma2, sigma1)
  #' 
  #' sigma1 <- matrix(c(0.5, 0, 0, 0, 0.4, 0, 0, 0, 0.3), nrow = 3)
  #' sigma2 <- diag(1, 3)
  #' # Case when all eigenvalues of sigma1 %*% solve(sigma2) are < 1
  #' kldcauchy(sigma1, sigma2)
  #' # Case when all eigenvalues of sigma1 %*% solve(sigma2) are > 1
  #' kldcauchy(sigma2, sigma1)
  #' }
  #' 
  #' @export
  
  # sigma1 and Sigma2 must be matrices
  if (is.numeric(sigma1) & !is.matrix(sigma1))
    sigma1 <- matrix(sigma1)
  if (is.numeric(sigma2) & !is.matrix(sigma2))
    sigma2 <- matrix(sigma2)
  
  # Number of variables
  p <- nrow(sigma1)
  
  # sigma1 and sigma2 must be square matrices with the same size
  if (ncol(sigma1) != p | nrow(sigma2) != p | ncol(sigma2) != p)
    stop("sigma1 et sigma2 doivent must be square matrices with rank p.")
  
  # IS sigma1 symmetric, positive-definite?
  if (!isSymmetric(sigma1))
    stop("sigma1 must be a symmetric, positive-definite matrix.")
  lambda1 <- eigen(sigma1, only.values = TRUE)$values
  if (any(lambda1 < .Machine$double.eps))
    stop("sigma1 must be a symmetric, positive-definite matrix.")
  
  # IS sigma2 symmetric, positive-definite?
  if (!isSymmetric(sigma2))
    stop("sigma2 must be a symmetric, positive-definite matrix.")
  lambda2 <- eigen(sigma2, only.values = TRUE)$values
  if (any(lambda2 < .Machine$double.eps))
    stop("sigma2 must be a symmetric, positive-definite matrix.")
  
  # Eigenvalues of sigma1 %*% inv(sigma2)
  lambda <- sort(eigen(sigma1 %*% solve(sigma2), only.values = TRUE)$values, decreasing = FALSE)
  
  prodlambda <- prod(lambda)
  
  k <- 5
  
  # M: data.frame of the indices for the nested sums
  # (i.e. all arrangements of n elements from {0:k})
  M <- expand.grid(rep(list(0:k), p))
  M <- M[-1, , drop = FALSE]
  Msum <- apply(M, 1, sum)
  kstep <- 5
  
  if (lambda[p] < 1) { # equations 109 & 102
    
    
    # The first 5 elements of the sum
    d <- 0
    for (i in 1:length(Msum)) {
      commun <- prod(
        sapply(1:p, function(j) {
          pochhammer(0.5, M[i, j])*(1 - lambda[j])^M[i, j]/factorial(M[i, j])
        })
      )
      d <- d + commun * pochhammer(1, Msum[i]) / ( pochhammer((1 + p)/2, Msum[i]) * Msum[i] )
    }
    
    # Next elements of the sum, until the expected precision
    k1 <- 1:k
    derive <- 0
    while (abs(d) > eps/10 & !is.nan(d)) {
      epsret <- signif(abs(d), 1)*10
      k <- k1[length(k1)]
      k1 <- k + (1:kstep)
      derive <- derive + d
      
      # M: data.frame of the indices for the nested sums
      M <- as.data.frame(matrix(nrow = 0, ncol = p))
      if (p > 1) {
        for (i in 1:(p-1)) {
          Mlist <- c( rep(list(0:k), p-i), rep(list(k1), i) )
          M <- rbind( M, expand.grid(Mlist) )
          for (j in 1:(p-1)) {
            Mlist <- Mlist[c(p, 1:(p-1))]
            M <- rbind(M, expand.grid(Mlist))
          }
        }
      }
      M <- rbind( M, expand.grid(rep(list(k1), p)) )
      Msum <- apply(M, 1, sum)
      
      d <- 0
      for (i in 1:length(Msum)) {
        commun <- prod(
          sapply(1:p, function(j) {
            pochhammer(0.5, M[i, j])*(1 - lambda[j])^M[i, j]/factorial(M[i, j])
          })
        )
        d <- d + commun * pochhammer(1, Msum[i]) / ( pochhammer((1 + p)/2, Msum[i]) * Msum[i] )
      }
      
    }
    
    # Next elements of the sum, with step=1, while not NaN
    if (is.nan(d)) {
      k1 <- k
      d <- 0
      while (!is.nan(d)) {
        if (d > 0)
          epsret <- signif(abs(d), 1)*10
        k <- k1
        k1 <- k + 1
        derive <- derive + d
        
        # M: data.frame of the indices for the nested sums
        M <- as.data.frame(matrix(nrow = 0, ncol = p))
        if (p > 1) {
          for (i in 1:(p-1)) {
            Mlist <- c( rep(list(0:k), p-i), rep(list(k1), i) )
            M <- rbind( M, expand.grid(Mlist) )
            for (j in 1:(p-1)) {
              Mlist <- Mlist[c(p, 1:(p-1))]
              M <- rbind(M, expand.grid(Mlist))
            }
          }
        }
        M <- rbind( M, rep(k1, p) )
        Msum <- apply(M, 1, sum)
        
        d <- 0
        for (i in 1:length(Msum)) {
          commun <- prod(
            sapply(1:p, function(j) {
              pochhammer(0.5, M[i, j])*(1 - lambda[j])^M[i, j]/factorial(M[i, j])
            })
          )
          d <- d + commun * pochhammer(1, Msum[i]) / ( pochhammer((1 + p)/2, Msum[i]) * Msum[i] )
        }
        
      }
    }
    
    result <- -0.5 * log(prodlambda) - (1 + p)/2 * derive
    
  } else if (lambda[1] > 1) { # equations 110
    
    # The first 5 elements of the sum
    d <- 0
    for (i in 1:length(Msum)) {
      commun <- prod(
        sapply(1:p, function(j) {
          pochhammer(0.5, M[i, j])*(1 - 1/lambda[j])^M[i, j]/factorial(M[i, j])
        })
      )
      A <- sum(1/(0:(Msum[i]-1) + (1+p)/2))
      d <- d - commun * A # / pochhammer((1 + p)/2, Msum[i])
    }
    
    # Next elements of the sum, until the expected precision
    k1 <- 1:k
    derive <- 0
    while (abs(d) > eps/10 & !is.nan(d)) {
      epsret <- signif(abs(d), 1)*10
      k <- k1[length(k1)]
      k1 <- k + (1:kstep)
      derive <- derive + d
      
      # M: data.frame of the indices for the nested sums
      M <- as.data.frame(matrix(nrow = 0, ncol = p))
      if (p > 1) {
        for (i in 1:(p-1)) {
          Mlist <- c( rep(list(0:k), p-i), rep(list(k1), i) )
          M <- rbind( M, expand.grid(Mlist) )
          for (j in 1:(p-1)) {
            Mlist <- Mlist[c(p, 1:(p-1))]
            M <- rbind(M, expand.grid(Mlist))
          }
        }
      }
      M <- rbind( M, expand.grid(rep(list(k1), p)) )
      Msum <- apply(M, 1, sum)
      
      d <- 0
      for (i in 1:length(Msum)) {
        commun <- prod(
          sapply(1:p, function(j) {
            pochhammer(0.5, M[i, j])*(1 - 1/lambda[j])^M[i, j]/factorial(M[i, j])
          })
        )
        A <- sum(1/(0:(Msum[i]-1) + (1+p)/2))
        d <- d - commun * A # / pochhammer((1 + p)/2, Msum[i])
      }
      
    }
    
    # Next elements of the sum, with step=1, while not NaN
    if (is.nan(d)) {
      k1 <- k
      d <- 0
      while (!is.nan(d)) {
        if (d > 0)
          epsret <- signif(abs(d), 1)*10
        k <- k1
        k1 <- k + 1
        derive <- derive + d
        
        # M: data.frame of the indices for the nested sums
        M <- as.data.frame(matrix(nrow = 0, ncol = p))
        if (p > 1) {
          for (i in 1:(p-1)) {
            Mlist <- c( rep(list(0:k), p-i), rep(list(k1), i) )
            M <- rbind( M, expand.grid(Mlist) )
            for (j in 1:(p-1)) {
              Mlist <- Mlist[c(p, 1:(p-1))]
              M <- rbind(M, expand.grid(Mlist))
            }
          }
        }
        M <- rbind( M, rep(k1, p) )
        Msum <- apply(M, 1, sum)
        
        d <- 0
        for (i in 1:length(Msum)) {
          commun <- prod(
            sapply(1:p, function(j) {
              pochhammer(0.5, M[i, j])*(1 - 1/lambda[j])^M[i, j]/factorial(M[i, j])
            })
          )
          A <- sum(1/(0:(Msum[i]-1) + (1+p)/2))
          d <- d - commun * A # / pochhammer((1 + p)/2, Msum[i])
        }
        
      }
    }
    
    result <- -0.5 * log(prodlambda) - (1 + p)/2 * prod(lambda)^-0.5 * derive
    
  } else { # equations 106 & 101
    
    # The first 5 elements of the sum
    d <- 0
    for (i in 1:length(Msum)) {
      commun <- prod(
        sapply(1:(p-1), function(j) {
          pochhammer(0.5, M[i, j])*(1 - lambda[j]/lambda[p])^M[i, j]/factorial(M[i, j])
        })
      )
      commun <- commun*(1 - 1/lambda[p])^M[i, p]/factorial(M[i, p])
      d <- d + commun * pochhammer(0.5, M[i, p])*pochhammer(1, Msum[i]) / ( pochhammer((1 + p)/2, Msum[i]) * Msum[i] )
    }
    
    # Next elements of the sum, until the expected precision
    k1 <- 1:k
    derive <- 0
    while (abs(d) > eps/10 & !is.nan(d)) {
      epsret <- signif(abs(d), 1)*10
      k <- k1[length(k1)]
      k1 <- k + (1:kstep)
      derive <- derive + d
      
      # M: data.frame of the indices for the nested sums
      M <- as.data.frame(matrix(nrow = 0, ncol = p))
      if (p > 1) {
        for (i in 1:(p-1)) {
          Mlist <- c( rep(list(0:k), p-i), rep(list(k1), i) )
          M <- rbind( M, expand.grid(Mlist) )
          for (j in 1:(p-1)) {
            Mlist <- Mlist[c(p, 1:(p-1))]
            M <- rbind(M, expand.grid(Mlist))
          }
        }
      }
      M <- rbind( M, expand.grid(rep(list(k1), p)) )
      Msum <- apply(M, 1, sum)
      
      d <- 0
      for (i in 1:length(Msum)) {
        commun <- prod(
          sapply(1:(p-1), function(j) {
            pochhammer(0.5, M[i, j])*(1 - lambda[j]/lambda[p])^M[i, j]/factorial(M[i, j])
          })
        )
        commun <- commun*(1 - 1/lambda[p])^M[i, p]/factorial(M[i, p])
        d <- d + commun * pochhammer(0.5, M[i, p])*pochhammer(1, Msum[i]) / ( pochhammer((1 + p)/2, Msum[i]) * Msum[i] )
      }
      
    }
    
    # Next elements of the sum, with step=1, while not NaN
    if (is.nan(d)) {
      k1 <- k
      d <- 0
      while (!is.nan(d)) {
        if (d > 0)
          epsret <- signif(abs(d), 1)*10
        k <- k1
        k1 <- k + 1
        derive <- derive + d
        
        # M: data.frame of the indices for the nested sums
        M <- as.data.frame(matrix(nrow = 0, ncol = p))
        if (p > 1) {
          for (i in 1:(p-1)) {
            Mlist <- c( rep(list(0:k), p-i), rep(list(k1), i) )
            M <- rbind( M, expand.grid(Mlist) )
            for (j in 1:(p-1)) {
              Mlist <- Mlist[c(p, 1:(p-1))]
              M <- rbind(M, expand.grid(Mlist))
            }
          }
        }
        M <- rbind( M, rep(k1, p) )
        Msum <- apply(M, 1, sum)
        
        d <- 0
        for (i in 1:length(Msum)) {
          commun <- prod(
            sapply(1:(p-1), function(j) {
              pochhammer(0.5, M[i, j])*(1 - lambda[j]/lambda[p])^M[i, j]/factorial(M[i, j])
            })
          )
          commun <- commun*(1 - 1/lambda[p])^M[i, p]/factorial(M[i, p])
          d <- d + commun * pochhammer(0.5, M[i, p])*pochhammer(1, Msum[i]) / ( pochhammer((1 + p)/2, Msum[i]) * Msum[i] )
        }
        
      }
    }
    
    result <- -0.5 * log(prodlambda) + (1 + p)/2 * (log(lambda[p]) - derive)
    
  }
  
  if (is.nan(d)) {
    epsret <- signif(epsret, 1)
    warning("Cannot reach the precision ", eps, " due to NaN\n",
            "Number of iteration: ", k, "\n",
            "Precision reached:", epsret)
    attr(result, "epsilon") <- epsret
  } else {
    attr(result, "epsilon") <- eps
  }
  attr(result, "k") <- k
  
  return(result)
}

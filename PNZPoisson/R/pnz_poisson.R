#' Cumulative Distribution Function of the PNZ-Poisson Distribution
#'
#' Computes the cumulative distribution function (CDF) of the PNZ-Poisson distribution for a vector of integer values.
#'
#' @param k_vec A vector of non-negative integers at which to evaluate the CDF.
#' @param lambda A positive parameter representing the mean of the distribution.
#' @param theta A positive parameter representing the dispersion of the distribution.
#'
#' @return A numeric vector of the same length as \code{k_vec}, containing the CDF evaluated at each element of \code{k_vec}.
#' @details The PNZ-Poisson distribution is a generalization of the Poisson distribution that introduces arbitrary underdispersion.
#' @examples
#' # Evaluate the CDF at k = 0:10 with lambda = 2 and theta = 1
#' ppnzpois(0:10, lambda = 2, theta = 1)
#' @export
ppnzpois <- function(k_vec, lambda, theta) {
  sapply(k_vec, ppnzpois_scalar, lambda = lambda, theta = theta)
}

# Internal function: Scalar CDF of the PNZ-Poisson distribution
ppnzpois_scalar <- function(k, lambda, theta) {
  if (k < 0 ) {
    return(0)
  }
  if ( k != floor(k)){
    k <- floor(k)
  }
  if (lambda <= 0 || theta <= 0) {
    stop("Parameters 'lambda' and 'theta' must be positive.")
  }
  
  kp1 <- k + 1
  if (k == 0) {
    return((kp1 - lambda) * (pgamma(lambda, shape = kp1 / theta, scale = theta, lower.tail = FALSE)) +
             kp1 * (dgamma(lambda / theta, shape = 1 + kp1 / theta, scale = 1)))
  }
  return(
    (kp1 - lambda) * (pgamma(lambda, shape = kp1 / theta, scale = theta, lower.tail = FALSE)) +
      kp1 * (dgamma(lambda / theta, shape = 1 + kp1 / theta, scale = 1)) -
      (k - lambda) * (pgamma(lambda, shape = k / theta, scale = theta, lower.tail = FALSE)) -
      k * (dgamma(lambda / theta, shape = 1 + k / theta, scale = 1))
  )
}

#' Probability Mass Function of the PNZ-Poisson Distribution
#'
#' Computes the probability mass function (PMF) of the PNZ-Poisson distribution for a vector of integer values.
#'
#' @param k_vec A vector of non-negative integers at which to evaluate the PMF.
#' @param lambda A positive parameter representing the mean of the distribution.
#' @param theta A positive parameter representing the dispersion of the distribution.
#' 
#' @return A numeric vector of the same length as \code{k_vec}, containing the PMF evaluated at each element of \code{k_vec}.
#' @details The PNZ-Poisson distribution allows for greater flexibility in modeling count data with over-dispersion compared to the standard Poisson distribution.
#' @examples
#' # Evaluate the PMF at k = 0:10 with lambda = 2 and theta = 1
#' dpnzpois(0:10, lambda = 2, theta = 1)
#' @export
dpnzpois <- function(k_vec, lambda, theta) {
  sapply(k_vec, dpnzpois_scalar, lambda = lambda, theta = theta)
}

# Internal function: Scalar PMF of the PNZ-Poisson distribution
dpnzpois_scalar <- function(k, lambda, theta) {
  if (k != floor(k)){
    stop("k must be an integer")
  }
  
  if (k < 0) {
    return(0)
  }
  
  return(ppnzpois_scalar(k, lambda, theta) - ppnzpois_scalar(k - 1, lambda, theta))
}

#' Random Generation from the PNZ-Poisson Distribution
#'
#' Generates random samples from the PNZ-Poisson distribution.
#'
#' @param n Integer specifying the number of random samples to generate.
#' @param lambda A positive parameter representing the mean of the distribution.
#' @param theta A positive parameter representing the dispersion of the distribution.
#'
#' @return An integer vector of length \code{n}, containing random samples from the PNZ-Poisson distribution.
#' @details The generated random sample follow the PNZ-Poisson distribution.
#' @examples
#' # Generate 10 random samples with lambda = 2 and theta = 1
#' rpnzpois(10, lambda = 2, theta = 1)
#' @export
rpnzpois <- function(n, lambda, theta) {
  # n: number of samples to generate
  # lambda, theta: parameters of the distribution
  sample <- integer(n)
  
  for (i in seq_len(n)) {
    u <- runif(1)  # Uniform random number between 0 and 1
    k <- -1        # Start from k = 0 after increment
    repeat {
      k <- k + 1
      if (ppnzpois_scalar(k, lambda, theta) >= u) {
        sample[i] <- k
        break
      }
    }
  }
  return(sample)
}








# R/qpnzpois.R

#' Quantile Function of the PNZ-Poisson Distribution
#'
#' Computes the quantile function of the PNZ-Poisson distribution for a vector of probabilities.
#'
#' @param p_vec A vector of probabilities (each between 0 and 1) at which to evaluate the quantiles.
#' @param lambda A positive parameter representing the mean of the distribution.
#' @param theta A positive parameter representing the dispersion of the distribution.
#' @param lower.tail Logical; if \code{TRUE} (default), probabilities are \eqn{P(X \leq x)}, otherwise, \eqn{P(X > x)}.
#' @param log.p Logical; if \code{TRUE}, probabilities \code{p} are given as \code{log(p)}.
#'
#' @return An integer vector of the same length as \code{p_vec}, containing the quantiles corresponding to each probability in \code{p_vec}.
#' @details The quantile function finds the smallest integer \eqn{k} such that the CDF at \eqn{k} is greater than or equal to the specified probability \eqn{p}.
#' @examples
#' # Find the 25th, 50th, and 75th percentiles with lambda = 2 and theta = 1
#' qpnzpois(c(0.25, 0.5, 0.75), lambda = 2, theta = 1)
#' @export
qpnzpois <- function(p_vec, lambda, theta, lower.tail = TRUE, log.p = FALSE) {
  # Validate input probabilities
  if (log.p) {
    p_vec <- exp(p_vec)
  }
  if (!lower.tail) {
    p_vec <- 1 - p_vec
  }
  
  if (any(p_vec < 0 | p_vec > 1)) {
    stop("All probabilities 'p_vec' must be between 0 and 1.")
  }
  
  sapply(p_vec, function(p) {
    if (p == 0) {
      return(0)
    }
    if (p == 1) {
      return(Inf)
    }
    k <- 0
    while (ppnzpois(k, lambda, theta) < p) {
      k <- k + 1
      # To prevent infinite loops in edge cases
      if (k > 1e6) {
        stop("Unable to find quantile. Please check the input parameters.")
      }
    }
    return(k)
  })
}


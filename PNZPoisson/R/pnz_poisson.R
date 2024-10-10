ppnzpois_scalar <- function(k, lambda,theta) {
  if (k < 0 ) {
    return(0)
  }
  if ( k != floor(k)){
    k <- floor(k)
  }
  if (lambda <= 0 || theta <= 0) {
    stop("Parameters 'lambda' and 'theta' must be positive.")
  }

  kp1=k+1
  if (k==0) {
    return((kp1-lambda)*(pgamma(lambda, shape = kp1 / theta, scale = theta,lower.tail = FALSE))+
             kp1*(dgamma(lambda/theta, shape=1+kp1/theta, scale=1)) )
  }
  return(
    (kp1-lambda)*(pgamma(lambda, shape = kp1 / theta, scale = theta,lower.tail = FALSE))+
      kp1*(dgamma(lambda/theta, shape=1+kp1/theta, scale=1))-
      (k-lambda)*(pgamma(lambda, shape = k / theta, scale = theta,lower.tail = FALSE))-
      k*(dgamma(lambda/theta, shape=1+k/theta, scale=1))

  )
}

# CDF of the PNZ-Poisson
ppnzpois <- function(k_vec, lambda, theta) {


  sapply(k_vec, ppnzpois_scalar, lambda = lambda, theta = theta)
}

dpnzpois_scalar <- function(k, lambda, theta) {
  if (k != floor(k)){
    stop("k must be an integer")
  }

  if (k < 0) {
    return(0)
  }

  return(ppnzpois_scalar(k, lambda, theta)-ppnzpois_scalar(k - 1, lambda, theta))
}

# PMF of the PNZ-Poisson
dpnzpois <- function(k_vec, lambda, theta) {
  sapply(k_vec, dpnzpois_scalar, lambda = lambda, theta = theta)
}

# Function to generate random samples from the PNZ-Poisson
rpnzpois <- function(n,lambda, theta) {
  # n: number of samples to generate
  # lambda, theta: parameters of the distribution
  sample <- integer(n)

  for (i in seq(n)) {
    u <- runif(1)  # Uniform random number between 0 and 1
    k <- -1        # Start from k = 0 after increment
    repeat {
      k <- k + 1
      if (ppnzpois_scalar(k, lambda, theta) >= u) {
        sample[i] <- k
        break}
    }
  }
  return(sample)
}


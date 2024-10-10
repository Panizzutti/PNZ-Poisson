





neg_log_likelihood_loglink <- function(params, X, y) {
  # Extract beta coefficients and theta
  p <- length(params) - 1
  betas <- params[1:p]
  theta <- params[length(params)]

  # Ensure theta is positive
  if (theta <= 0) {
    return(Inf)  # Return infinity to penalize non-positive theta
  }
  #print("gaherqgreger")
  #print(c(betas,theta))

  # Compute linear predictor and lambda
  eta <- X %*% betas
  lambda <- exp(eta)


  # Use sapply to compute dpnzpois for each observation
  Li <- sapply(seq_along(y), function(i) {
    if (lambda[i] == 0) {
      1
    } else {
      dpnzpois_scalar(y[i], lambda[i], theta)
    }
  })

  # Handle zero or negative PMF values to avoid log(0) or log of negative numbers
  # Add a small epsilon to Li to prevent log(0)
  epsilon <- 1e-64
  Li[Li <= 0 | is.na(Li)] <- epsilon

  # Compute the negative log-likelihood
  negLL <- -sum(log(Li))
  return(negLL)
}

# Define the glm.pnz function
glm.pnz <- function(formula, data, ...) {
  # Extract response and design matrix
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  X <- model.matrix(formula, data)

  # Fit the model using mle_pnzpois_loglink
  fit_result <- mle_pnzpois_loglink(X, y, ...)

  # Store additional information
  fit_result$call <- match.call()
  fit_result$formula <- formula
  fit_result$data <- data

  # Set the class of the object
  class(fit_result) <- "pnz_glm"

  return(fit_result)
}

# Summary method for pnz_glm objects
summary.pnz_glm <- function(object, ...) {
  # Extract estimates and standard errors
  estimates <- object$betas
  std_errors <- object$std_errors[1:length(estimates)]

  # Calculate z-values and p-values
  z_values <- estimates / std_errors
  p_values <- 2 * (1 - pnorm(abs(z_values)))

  # Create a coefficients table
  coefficients_table <- data.frame(
    Estimate = estimates,
    'Std. Error' = std_errors,
    'z value' = z_values,
    'Pr(>|z|)' = p_values,
    check.names = FALSE  # Prevents R from altering column names
  )

  # Add row names as parameter names
  param_names <- colnames(object$X)
  rownames(coefficients_table) <- param_names

  # Extract additional measures
  logLik <- object$logLik
  AIC <- object$AIC
  BIC <- object$BIC
  Res_deviance<-object$Res_deviance
  Null_deviance<-object$Null_deviance
  df <- object$df

  # Output the summary
  summary_output <- list(
    call = object$call,
    coefficients = coefficients_table,
    logLik = logLik,
    AIC = AIC,
    BIC = BIC,
    Res_deviance = Res_deviance,
    Null_deviance = Null_deviance,
    df = df,
    theta = object$theta,
    convergence = object$convergence
  )

  class(summary_output) <- "summary.pnz_glm"
  return(summary_output)
}

# Print method for summary.pnz_glm objects
print.summary.pnz_glm <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nCoefficients:\n")
  printCoefmat(x$coefficients, digits = 4, signif.stars = TRUE)

  cat("\nDispersion parameter (theta):")
  print(x$theta[1])

  cat("Log-likelihood:")
  print(x$logLik)

  cat("\nNull deviance:")
  print(x$Null_deviance)

  cat("\nResidual deviance:")
  print(x$Res_deviance)

  cat("\nAIC:")
  print(x$AIC)

  cat("BIC:")
  print(x$BIC)

  if (x$convergence == 0) {
    cat("\nOptimization converged successfully.\n")
  } else {
    cat("\nOptimization did not converge.\n")
  }
}


MLE_find <- function(X, y, start_params = NULL, method = "L-BFGS-B", max_retries = 5) {
  n <- nrow(X)
  p <- ncol(X)

  # Generate initial starting values if not provided
  generate_starting_values <- function() {

    glm_fit <- glm(y ~ X-1, family = poisson(link = "log"))


    start_betas <- coef(glm_fit)
    pearson_residuals <- residuals(glm_fit, type = "pearson")
    start_theta <- sum(pearson_residuals^2) / (n - p)
    if (start_theta <= 0 || is.na(start_theta) || is.infinite(start_theta)) {
      start_theta <- 1
      warning("Calculated start_theta is non-positive or invalid. Using start_theta = 1.")
    }
    c(start_betas, start_theta)
  }

  # Get initial estimates from the Poisson model
  initial_estimates <- if (is.null(start_params)) {
    generate_starting_values()
  } else {
    start_params
  }

  if (length(initial_estimates) != p + 1) {
    stop("Length of start_params must be equal to number of betas plus one (theta).")
  }

  # Set lower and upper bounds
  lower_bounds <- c(rep(-Inf, p), 1e-6)  # Theta > 0
  upper_bounds <- c(rep(Inf, p), Inf)



  # Optimization with retry mechanism
  attempt <- 0
  while (attempt <= max_retries) {
    # Initialize var_cov_matrix and std_errors at the beginning of each attempt
    var_cov_matrix <- NULL
    std_errors <- rep(NA, p + 1)

    # Generate new starting values close to initial estimates
    if (attempt == 1) {
      # Use initial estimates on the first attempt
      start_params <- initial_estimates
    } else {
      # Perturb initial estimates slightly for subsequent attempts
      # Decrease perturbation with each attempt
      start_betas <- initial_estimates[1:p] + 0.1*attempt*rnorm(p)
      start_theta <- initial_estimates[p + 1] + 0.1*attempt*rnorm(1)
      # Ensure theta remains positive
      if (start_theta <= 0) {
        start_theta <- initial_estimates[p + 1] * (1 + 0.1*attempt)
      }
      start_params <- c(start_betas, start_theta)
    }


    fit <- tryCatch({
      optim(
        par = start_params,
        fn = neg_log_likelihood_loglink,
        X = X,
        y = y,
        method = method,
        lower = lower_bounds,
        upper = upper_bounds,
        control = list(
          maxit = 1000,
          parscale = abs(start_params),
          pgtol = 1e-8,
          REPORT = 1  # Reports progress every iteration
        ),
        hessian = TRUE
      )
    }, error = function(e) {
      NULL
    })

    # Check if optimization was successful
    if (!is.null(fit) && fit$convergence == 0 && is.finite(fit$value)) {
      # Attempt to compute the variance-covariance matrix
      if (!is.null(fit$hessian) && all(dim(fit$hessian) == c(p + 1, p + 1))) {
        var_cov_matrix <- tryCatch(solve(fit$hessian), error = function(e) NULL)
        if (!is.null(var_cov_matrix)) {
          # Hessian inversion successful
          std_errors <- sqrt(diag(var_cov_matrix))
          # Optimization and Hessian inversion successful, exit loop
          break
        } else {
          # Hessian inversion failed, retry optimization
          warning(sprintf("Optimization attempt %d: Hessian is singular. Trying new starting values.", attempt))
        }
      } else {
        # Hessian not available or incorrect dimensions, retry optimization
        warning(sprintf("Optimization attempt %d: Hessian not available or has incorrect dimensions. Trying new starting values.", attempt))
      }
    } else {
      # Optimization failed, retry with new starting values
      warning(sprintf("Optimization attempt %d failed. Trying new starting values.", attempt))
    }

    attempt <- attempt + 1
  }

  # If optimization still failed after max_retries
  if (attempt > max_retries || is.null(fit) || fit$convergence != 0 || !is.finite(fit$value) || is.null(var_cov_matrix)) {
    stop("Optimization failed after multiple attempts. Consider checking the data or model specification.")
  }

  fit$std_errors= std_errors
  fit$var_cov_matrix= var_cov_matrix
  return(fit)


}




mle_pnzpois_loglink <- function(X, y, start_params = NULL, method = "L-BFGS-B", max_retries = 5) {
  n <- nrow(X)
  p <- ncol(X)

  fit= MLE_find(X, y,start_params=start_params, method = method, max_retries = max_retries)
  NULLfit= MLE_find( as.matrix(rep(1,n),ncol=1), y,start_params=start_params, method = method, max_retries = max_retries)
  # Extract estimated parameters
  est_params <- fit$par
  betas_est <- est_params[1:p]
  theta_est <- est_params[p + 1]

  logLikelihood = -fit$value

  satLi <- sapply(seq_along(y), function(i) {
    if (y[i] == 0) {
      1
    } else {
      dpnzpois_scalar(y[i], y[i], theta_est)
    }
  })


  epsilon <- 1e-64
  satLi[satLi <= 0 | is.na(satLi)] <- epsilon
  saturatedLL= sum(log(satLi))

  resdeviance= 2*(saturatedLL - logLikelihood)


  nullLL=-NULLfit$value
  nulldeviance=2*(saturatedLL-nullLL)




  # Compute AIC and BIC
  AIC = 2 * (p+1) - 2 * logLikelihood
  BIC = log(n) * (p+1) - 2 * logLikelihood

  # Compile results
  result <- list(
    betas = betas_est,
    theta = theta_est,
    std_errors = fit$std_errors,
    logLik = logLikelihood,
    AIC = 2 * (p+1) - 2 * logLikelihood,
    BIC = log(n) * (p+1) - 2 * logLikelihood,
    Res_deviance = resdeviance,
    Null_deviance = nulldeviance,
    convergence = fit$convergence,
    var_cov_matrix = fit$var_cov_matrix,
    X = X,
    y = y
  )

  return(result)
}





# Jacob Townson

## With how often we need to simulate data for our model to test with, I thought
## it would be useful to have the simulation code in functions that we can call upon
## whenever it's needed. This way we don't have to keep re-coding the same stuff
## repeatedly. 

require(plyr)
require(dplyr)
require(ggplot2)

fcurve <- function(MW = 1:20, C, b){
  indiv_flow = C*(as.numeric(MW)^b)
  return(indiv_flow)
}

RSSfunc <- function(dat,x,C,b) sum((dat-(C*(x^b)))^2)

# Simulating data (no changepoint)

## n = Number of x values, N = number of observations at each x,
## C = The constant value in our flow function, b = exponent value in flow function, 
## sd = standard deviation of the simulated data

no_cp_sim <- function(N = 20, n = 20, C = 10, b = .3, sd = 1){
  
  x = sort(rep(1:n, N))
  
  # Simulate Data
  
  f = C*(x^b)
  sim_flows=f+rnorm(N*n,sd=sd)
  sim_dat = data.frame(x, sim_flows)
  nums = unique(x) # Note because of this line, cp cannot be 2 or N-2 for these sims
  
  return(sim_dat)
}

# Simulating data (no changepoint), with x values between 0 and 1

## n = Number of x values, N = number of observations at each x,
## C = The constant value in our flow function, b = exponent value in flow function, 
## cp = changepoint, sd = standard deviation of the simulated data

no_cp_sim_zero2one <- function(N = 20, n = 20, C = 10, b = .3, sd = 1){
  
  x = sort(rep((1:n)/n, N))
  
  # Simulate Data
  
  f = C*(x^b)
  sim_flows=f+rnorm(N*n,sd=sd)
  sim_dat = data.frame(x, sim_flows)
  nums = unique(x) # Note because of this line, cp cannot be 2 or N-2 for these sims
  
  return(sim_dat)
}

# Simulating data (with changepoint)

## n = Number of x values, N = number of observations at each x,
## C1-2 = The constant value in our flow function before and after the cp,
## b1-2 = exponent value in flow function before and after the cp, 
## cp = changepoint, sd = standard deviation of the simulated data

cp_sim <- function(N = 20, n = 20, C1 = 10, b1 = .3, C2 = 15, b2 = .5, 
                   cp=5, sd = 1){
  
  x = sort(rep(1:n, N))
  
  # Simulate Data
  
  f = ifelse(x>cp, C2, C1)*(x^(ifelse(x>cp, b2, b1)))
  sim_flows=f+rnorm(N*n,sd=sd)
  sim_dat = data.frame(x, sim_flows)
  nums = unique(x) # Note because of this line, cp cannot be 2 or N-2 for these sims
  
  return(sim_dat)
}

# nonlin.NRg

## This function is our version of the NR algorithm optimization. It minimizes the 
## RSS for the b value, and then calculates the correlating C value using our mathematical
## properties that we found for the model. 
## x = observed x values
## y = observed output values
## b.start = the starting value for estimating b
## maxiter = maximum number of iterations to go through when minimizing the model
## print.output = choose whether or not to have the output of each iteration printed to the console
## b_bounds = if you know enough information about the model, this allows you to play with the bounds of b

nonlin.NRg <- function(x, y, b.start=0.5, tol=1e-6, maxiter=1000, print.output=FALSE, b_bounds = c(.001, .99)){
  b=b.start
  if(print.output == TRUE){
    cat("Starting value: b=",b.start,"\n")
  }
  b.last=b+1
  i=0
  X=matrix(x,length(x),length(x))
  Y=matrix(y,length(y),length(y))
  boundtest1 = FALSE
  boundtest2 = FALSE
  doublebound = FALSE
  while ((abs(b-b.last)>tol) & (i<maxiter) & !doublebound){
    b.last=b
    g=sum((X*t(X)^2)^b*Y*log(X/t(X)))
    dg=sum((X*(t(X)^2))^b*log(X*(t(X)^2))*Y*log(X/t(X)))
    b=b-g/dg
    if(b < b_bounds[1] & !boundtest1){
      b = b_bounds[2]
      boundtest1 = TRUE
    } else if(b > b_bounds[2] & !boundtest2){
      b = b_bounds[1]
      boundtest2 = TRUE
    } else if(b < b_bounds[1] & boundtest1){
      b = b_bounds[1]
      doublebound = TRUE
    } else if(b > b_bounds[2] & boundtest2){
      b = b_bounds[2]
      doublebound = TRUE
    }
    i=i+1
    if(print.output == TRUE){
      cat("Iteration",i,": b=",b,"\n")
    }
  }
  if (i==maxiter){
    cat("Warning: reached the maximum number of iterations\n")
  }
  C=sum(x^b*y*log(x))/sum(x^(2*b)*log(x))
  RSS = RSSfunc(y, x, C, b)
  return(list('C' = C, 'b' = b, 'RSS' = RSS))
}

# nonlin.NRg with CP

## This function is our version of the NR algorithm optimization WITH a CP. It minimizes the 
## RSS for the b value on either side, 
### and then calculates the correlating C value using our mathematical
## properties that we found for the model. 
## x = observed x values
## y = observed output values
## b.start = the starting value for estimating b on both sides of the CP
## maxiter = maximum number of iterations to go through when minimizing the model
## print.output = choose whether or not to have the output of each iteration printed to the console
## b_bounds = if you know enough information about the model, this allows you to play with the bounds of b

nonlin.NRg_CP <- function(x, y, b.start=0.5, tol=1e-6, maxiter=1000, print.output=FALSE, b_bounds = c(.001, .99)){
  # Function to find optimal b and C using Newton-Raphson for a given subset of data
  fit_segment <- function(x_seg, y_seg, b.start) {
    b <- b.start
    b.last <- b + 1
    i <- 0
    X <- matrix(x_seg, length(x_seg), length(x_seg))
    Y <- matrix(y_seg, length(y_seg), length(y_seg))
    boundtest1 <- FALSE
    boundtest2 <- FALSE
    doublebound <- FALSE
    
    while ((abs(b-b.last) > tol) & (i < maxiter) & !doublebound) {
      b.last <- b
      g <- sum((X*t(X)^2)^b*Y*log(X/t(X)))
      dg <- sum((X*(t(X)^2))^b*log(X*(t(X)^2))*Y*log(X/t(X)))
      b <- b - g/dg
      
      if(b < b_bounds[1] & !boundtest1) {
        b <- b_bounds[2]
        boundtest1 <- TRUE
      } else if(b > b_bounds[2] & !boundtest2) {
        b <- b_bounds[1]
        boundtest2 <- TRUE
      } else if(b < b_bounds[1] & boundtest1) {
        b <- b_bounds[1]
        doublebound <- TRUE
      } else if(b > b_bounds[2] & boundtest2) {
        b <- b_bounds[2]
        doublebound <- TRUE
      }
      
      i <- i + 1
      if(print.output == TRUE) {
        cat("Iteration", i, ": b=", b, "\n")
      }
    }
    
    if (i == maxiter && print.output) {
      cat("Warning: reached the maximum number of iterations\n")
    }
    
    C <- sum(x_seg^b*y_seg*log(x_seg))/sum(x_seg^(2*b)*log(x_seg))
    return(list('C' = C, 'b' = b))
  }
  
  # Get unique x values to properly identify potential changepoints
  unique_x <- sort(unique(x))
  
  if(length(unique_x) < 5) {
    stop("Not enough unique x values. Need at least 5 different x values to find a changepoint.")
  }
  
  # Find possible changepoint locations based on unique x values
  min_unique_x <- unique_x[3]
  max_unique_x <- unique_x[length(unique_x) - 2]
  
  if(print.output) {
    cat("Testing changepoints between x =", min_unique_x, "and x =", max_unique_x, "\n")
  }
  
  # Initialize variables
  min_total_rss <- Inf
  best_cp <- NULL
  best_left_fit <- NULL
  best_right_fit <- NULL
  
  # Testing all possible changepoints
  for(cp in unique_x[unique_x >= min_unique_x & unique_x <= max_unique_x]) {
    if(print.output) {
      cat("Testing changepoint at x =", cp, "\n")
    }
    
    # Split data
    left_indices <- which(x <= cp)
    right_indices <- which(x > cp)
    
    left_x <- x[left_indices]
    left_y <- y[left_indices]
    right_x <- x[right_indices]
    right_y <- y[right_indices]
    
    # Check our x's
    if(length(unique(left_x)) < 3 || length(unique(right_x)) < 3) {
      if(print.output) {
        cat("  Skipping: Not enough unique x values on both sides\n")
      }
      next
    }
    
    # Model build on each side
    left_fit <- fit_segment(left_x, left_y, b.start)
    right_fit <- fit_segment(right_x, right_y, b.start)
    
    # Calculate RSS
    left_rss <- RSSfunc(left_y, left_x, left_fit$C, left_fit$b)
    right_rss <- RSSfunc(right_y, right_x, right_fit$C, right_fit$b)
    ############# This sum works for the case where SD is the same on both sides of the CP, but if we have different SDs, ############# 
    ############# we will need to make these RSS values weighted based on the properties on each side of the CP.          ############# 
    total_rss <- left_rss + right_rss
    
    if(print.output) {
      cat("  Left RSS:", left_rss, "Right RSS:", right_rss, "Total:", total_rss, "\n")
    }
    
    # Update best changepoint if we found a better split
    if(total_rss < min_total_rss) {
      min_total_rss <- total_rss
      best_cp <- cp
      best_left_fit <- left_fit
      best_right_fit <- right_fit
    }
  }
  
  if(is.null(best_cp)) {
    stop("Could not find a valid changepoint.")
  }
  
  if(print.output) {
    cat("\nBest changepoint found at x =", best_cp, "\n")
    cat("Segment 1 (x ≤", best_cp, "): C =", best_left_fit$C, "b =", best_left_fit$b, "\n")
    cat("Segment 2 (x >", best_cp, "): C =", best_right_fit$C, "b =", best_right_fit$b, "\n")
    cat("Total RSS:", min_total_rss, "\n")
  }
  
  return(list(
    'changepoint' = best_cp,
    'C1' = best_left_fit$C,
    'b1' = best_left_fit$b,
    'C2' = best_right_fit$C,
    'b2' = best_right_fit$b,
    'total_rss' = min_total_rss
  ))
}






# nonlin.NRg with CP - GPU Accelerated

## This function is our version of the NR algorithm optimization WITH a CP.  
## It minimizes the RSS for the b value on either side, 
### and then calculates the correlating C value using our mathematical
## properties that we found for the model. This version uses gpuR to attempt
## speeding up the computation time (and maybe ease the strain on my CPU). 
## x = observed x values
## y = observed output values
## b.start = the starting value for estimating b on both sides of the CP
## maxiter = maximum number of iterations to go through when minimizing the model
## print.output = choose whether or not to have the output of each iteration printed to the console
## b_bounds = if you know enough information about the model, this allows you to play with the bounds of b

nonlin.NRg_CP_GPU <- function(x, y, b.start=0.5, tol=1e-6, maxiter=1000, print.output=FALSE, b_bounds = c(.001, .99)) {
# Check if gpuR is installed
  if (!requireNamespace("gpuR", quietly = TRUE)) {
    stop("Package 'gpuR' is needed for GPU acceleration. Please install it first.")
  }
  
  library(gpuR)
  
  # Function to find optimal b and C using Newton-Raphson on GPU for a segment
  fit_segment_gpu <- function(x_seg, y_seg, b.start) {
    b <- b.start
    b.last <- b + 1
    i <- 0
    
    # Move data to GPU
    X_gpu <- gpuMatrix(matrix(x_seg, length(x_seg), length(x_seg)), type="float")
    Y_gpu <- gpuMatrix(matrix(y_seg, length(y_seg), length(y_seg)), type="float")
    
    # Precompute log(X) on GPU
    logX_gpu <- gpuMatrix(log(matrix(x_seg, length(x_seg), length(x_seg))), type="float")
    
    boundtest1 <- FALSE
    boundtest2 <- FALSE
    doublebound <- FALSE
    
    while ((abs(b-b.last) > tol) & (i < maxiter) & !doublebound) {
      b.last <- b
      
      # Compute X*t(X)^2 on GPU
      XtX2_gpu <- X_gpu * (t(X_gpu)^2)
      
      # Compute (X*t(X)^2)^b on GPU
      XtX2_pow_b_gpu <- XtX2_gpu^b
      
      # Compute log(X/t(X)) on GPU
      logRatio_gpu <- logX_gpu - t(logX_gpu)
      
      # Compute g on GPU
      g_components_gpu <- XtX2_pow_b_gpu * Y_gpu * logRatio_gpu
      g <- sum(as.matrix(g_components_gpu))
      
      # Compute log(X*(t(X)^2)) on GPU
      logXtX2_gpu <- log(XtX2_gpu)
      
      # Compute dg on GPU
      dg_components_gpu <- (XtX2_pow_b_gpu * logXtX2_gpu) * Y_gpu * logRatio_gpu
      dg <- sum(as.matrix(dg_components_gpu))
      
      b <- b - g/dg
      
      # Apply bounds
      if(b < b_bounds[1] & !boundtest1) {
        b <- b_bounds[2]
        boundtest1 <- TRUE
      } else if(b > b_bounds[2] & !boundtest2) {
        b <- b_bounds[1]
        boundtest2 <- TRUE
      } else if(b < b_bounds[1] & boundtest1) {
        b <- b_bounds[1]
        doublebound <- TRUE
      } else if(b > b_bounds[2] & boundtest2) {
        b <- b_bounds[2]
        doublebound <- TRUE
      }
      
      i <- i + 1
      if(print.output == TRUE) {
        cat("Iteration", i, ": b=", b, "\n")
      }
    }
    
    if (i == maxiter && print.output) {
      cat("Warning: reached the maximum number of iterations\n")
    }
    
    # Move x_seg and y_seg to GPU for final C calculation
    x_seg_gpu <- gpuVector(x_seg, type="float")
    y_seg_gpu <- gpuVector(y_seg, type="float")
    
    # Calculate C using GPU computations
    x_pow_b_gpu <- x_seg_gpu^b
    log_x_gpu <- log(x_seg_gpu)
    num <- sum(as.vector(x_pow_b_gpu * y_seg_gpu * log_x_gpu))
    denom <- sum(as.vector(x_pow_b_gpu^2 * log_x_gpu))
    
    C <- num/denom
    
    return(list('C' = C, 'b' = b))
  }
  
  # Get unique x values
  unique_x <- sort(unique(x))
  
  if(length(unique_x) < 5) {
    stop("Not enough unique x values. Need at least 5 different x values to find a changepoint.")
  }
  
  # Find possible changepoint locations
  cp_candidates <- unique_x[3:(length(unique_x) - 2)]
  
  if(print.output) {
    cat("Testing", length(cp_candidates), "potential changepoints\n")
  }
  
  # Initialize variables for tracking results
  min_total_rss <- Inf
  best_cp <- NULL
  best_left_fit <- NULL
  best_right_fit <- NULL
  
  # Calculate RSS function using GPU
  RSSfunc_gpu <- function(y_actual, x_vals, C, b) {
    y_pred <- C * (x_vals^b)
    return(sum((y_actual - y_pred)^2))
  }
  
  # Batch processing of changepoints
  # For larger datasets, process in batches to avoid GPU memory issues
  batch_size <- min(100, length(cp_candidates))
  num_batches <- ceiling(length(cp_candidates) / batch_size)
  
  for(batch_idx in 1:num_batches) {
    start_idx <- (batch_idx - 1) * batch_size + 1
    end_idx <- min(batch_idx * batch_size, length(cp_candidates))
    batch_cp_candidates <- cp_candidates[start_idx:end_idx]
    
    if(print.output) {
      cat("Processing batch", batch_idx, "of", num_batches, "\n")
    }
    
    # Process each changepoint in the batch
    for(cp in batch_cp_candidates) {
      if(print.output) {
        cat("Testing changepoint at x =", cp, "\n")
      }
      
      # Split data
      left_indices <- which(x <= cp)
      right_indices <- which(x > cp)
      
      left_x <- x[left_indices]
      left_y <- y[left_indices]
      right_x <- x[right_indices]
      right_y <- y[right_indices]
      
      # Check for sufficient data points
      if(length(unique(left_x)) < 3 || length(unique(right_x)) < 3) {
        if(print.output) {
          cat("  Skipping: Not enough unique x values on both sides\n")
        }
        next
      }
      
      # Fit segments using GPU acceleration
      left_fit <- fit_segment_gpu(left_x, left_y, b.start)
      right_fit <- fit_segment_gpu(right_x, right_y, b.start)
      
      # Calculate RSS
      left_rss <- RSSfunc_gpu(left_y, left_x, left_fit$C, left_fit$b)
      right_rss <- RSSfunc_gpu(right_y, right_x, right_fit$C, right_fit$b)
      total_rss <- left_rss + right_rss
      
      if(print.output) {
        cat("  Left RSS:", left_rss, "Right RSS:", right_rss, "Total:", total_rss, "\n")
      }
      
      # Update best changepoint
      if(total_rss < min_total_rss) {
        min_total_rss <- total_rss
        best_cp <- cp
        best_left_fit <- left_fit
        best_right_fit <- right_fit
      }
    }
    
    # Clean up GPU memory after each batch
    gpuR::gpuClean()
  }
  
  if(is.null(best_cp)) {
    stop("Could not find a valid changepoint.")
  }
  
  if(print.output) {
    cat("\nBest changepoint found at x =", best_cp, "\n")
    cat("Segment 1 (x ≤", best_cp, "): C =", best_left_fit$C, "b =", best_left_fit$b, "\n")
    cat("Segment 2 (x >", best_cp, "): C =", best_right_fit$C, "b =", best_right_fit$b, "\n")
    cat("Total RSS:", min_total_rss, "\n")
  }
  
  # Final cleanup of GPU resources
  gpuR::gpuClean()
  
  return(list(
    'changepoint' = best_cp,
    'C1' = best_left_fit$C,
    'b1' = best_left_fit$b,
    'C2' = best_right_fit$C,
    'b2' = best_right_fit$b,
    'total_rss' = min_total_rss
  ))
}
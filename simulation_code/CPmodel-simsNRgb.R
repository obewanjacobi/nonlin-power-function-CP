# Jacob Townson
# Nonlinear Changepoint Detection Simulations - Newton-Raphson Algo (Using the g(b) optimization function)
# g(b) = sum(sum(x_i x_j^2)^b y_i ln(x_i / x_j))

# The nonlin.NRg() function was built by Dr. Gill, and optimizes a nonlinear model using the NR method but simplifies
# our model to only have 1 parameter, namely b. This code will run just like the findCP_nonlin() function built in the 
# nonlinear_model-simsLM.R file.

# Predict Nonlinear Output Values
## I didn't like the way the predict function behaves for the nlsLM() function, so I'm making a personalized one

fcurve <- function(MW = 1:20, C, b){
  indiv_flow = C*(as.numeric(MW)^b)
  return(indiv_flow)
}

# Finding the change point and confidence
## R = number of simulations, J = Number of x values, N = number of observations at each x,
## C = The constant value in our flow function, b = exponent value in flow function, delta = change factor at cp
## cp = changepoint

findCP_NRg <- function(R = 500, N = 20, J = 20, C = 10, b = .3, delta = 1.25, cp=4.5, sd = 1, maxiter = 1000,
                       write_results = TRUE){
  
  require(plyr)
  require(dplyr)
  require(ggplot2)
  source('./algorithm_studies/NR_DrGill_algo.R')
  RSSfunc <- function(dat,x,C,b) sum((dat-(C*(x^b)))^2)
  
  n = J # will fix this to be 1 variable later
  j = 1
  x = sort(rep(1:J, N))
  is.CI.correct=rep(0,R)
  intervals = list(c(0,0)) # initialize a list where we can store the interval with the determined cp based on the 
  # minimum RSS.
  p = list(0)              # initializes a list where we can store plots for each simulation
  
  while(j <= R){
    f=C*x^b*delta* ifelse(x>cp, 1, 1/delta)
    sim_flows=f+rnorm(N*J,sd=sd)
    test = list(c(0,0))    # initialize a list where we can store RSS for each test.
    mod_dat = data.frame(x, sim_flows)
    C_left = 5
    b_left = .5
    C_right = 5
    b_right = .5
    nums = 2:(n-2)
    
    for(i in sample(nums)){          # i must be 2 or 1 less than n because you can't build a model off of 1 point of
      lhs = filter(mod_dat, x <= i)  # observations
      rhs = filter(mod_dat, x > i)
      
      left_mod = nonlin.NRg(lhs$x, lhs$sim_flows, b.start=b_left, tol=1e-6, maxiter=maxiter, print.output=FALSE)
      left_RSS = RSSfunc(lhs$sim_flows, lhs$x, left_mod$C, left_mod$b)
      right_mod = nonlin.NRg(rhs$x, rhs$sim_flows, b.start=b_right, tol=1e-6, maxiter=maxiter, print.output=FALSE)
      right_RSS = RSSfunc(rhs$sim_flows, rhs$x, right_mod$C, right_mod$b)
      
      test[[i-1]] = c(left_RSS, right_RSS)
      
      # predicted C and b values to input into start value of g algo next round
      C_left = left_mod$C
      b_left = left_mod$b
      C_right = right_mod$C
      b_right = right_mod$b
      
      if(j == R){ # stores the last simulations plots for study
        mod_dat$predictsL = fcurve(mod_dat$x, left_mod$C, left_mod$b)
        mod_dat$predictsR = fcurve(mod_dat$x, right_mod$C, right_mod$b)
        p[[i-1]] = ggplot(mod_dat, aes(x, sim_flows)) + geom_point() +
          geom_line(data = subset(mod_dat, x <= cp), aes(y = predictsL), size = 1, color = 'red') +
          geom_line(data = subset(mod_dat, x >= cp), aes(y = predictsR), size = 1, color = 'blue') +
          labs(x = '', y = 'Simulated Flows', title = 'Pipeline Changepoint Simulation')
        # in this list p[[cp - 1]] will be the ideal plot
      }
    }
    
    ## Not super clear on how to find minimum
    ## Attempt #1: find the mean of the left and right RSS values and take the lowest one
    k=1
    rs = c(0)
    for(k in 1:length(test)){
      rs[k] = sum(test[[k]])
    }
    intervals[[j]] = c(which.min(rs)+1, which.min(rs)+2)
    if(write_results == TRUE){
      cat(j,"th simulated data set: CI= ",intervals[[j]],"\n")
    }
    is.CI.correct[j] = intervals[[j]][1] <= cp && intervals[[j]][2] >= cp
    
    j = j+1
  }
  
  return(list('intervals' = intervals, 'is.CI.correct' = is.CI.correct, 'plots' = p,
              'confidence' = mean(is.CI.correct)))
}
# Jacob Townson

# CP consistency simulations

# This script will be comprised of a functions that does 2 things. 
# 1 to simulate data with a changepoint. 
# The next part will try and predict the model, assuming that the
# changepoint is plus or minus a certain distance from the changepoint. 
# The goal is to see if Q is minimized when the correct changepoint is chosen.

fcurve <- function(MW = 1:20, C, b){
  indiv_flow = C*(as.numeric(MW)^b)
  return(indiv_flow)
}

## n = Number of x values, N = number of observations at each x,
## C = The constant value in our flow function, b = exponent value in flow function, 
## delta = change factor at cp, cp = changepoint

# NOTE: n must be > 6, and cp must be a number equal to or between 3 and n-3

consistency_tester <- function(N = 20, n = 20, C1 = 10, b1 = .3, C2 = 12, b2 = .25, 
                               cp=5, sd = 1, maxiter = 1000,
                               write_results = FALSE){
  require(plyr)
  require(dplyr)
  require(ggplot2)
  source('./algorithm_studies/NR_DrGill_algo.R')
  RSSfunc <- function(dat,x,C,b) sum((dat-(C*(x^b)))^2)
  
  j = 1
  x = sort(rep(1:n, N))
  p = list(0)              # initializes a list where we can store plots for each simulation
  
  # Simulate Data
  
  f = ifelse(x>cp, C2, C1)*(x^(ifelse(x>cp, b2, b1)))
  sim_flows=f+rnorm(N*n,sd=sd)
  test = list(c(0,0,0))    # initialize a list where we can store RSS for each test.
  mod_dat = data.frame(x, sim_flows)
  C_left = 5
  b_left = .5
  C_right = 5
  b_right = .5
  nums = unique(x) # Note because of this line, cp cannot be 2 or N-2 for these sims
  
  for(i in 3:(length(nums)-3)){
    temp = nums[i]
    lhs = filter(mod_dat, x <= temp)
    rhs = filter(mod_dat, x > temp)
    
    left_mod = nonlin.NRg(lhs$x, lhs$sim_flows, b.start=b_left, tol=1e-6, 
                          maxiter=maxiter, print.output=write_results)
    left_RSS = RSSfunc(lhs$sim_flows, lhs$x, left_mod$C, left_mod$b)
    right_mod = nonlin.NRg(rhs$x, rhs$sim_flows, b.start=b_right, tol=1e-6, 
                           maxiter=maxiter, print.output=write_results)
    right_RSS = RSSfunc(rhs$sim_flows, rhs$x, right_mod$C, right_mod$b)
    total_RSS = left_RSS+right_RSS
    
    test[[i-2]] = c('left_RSS' = left_RSS, 'right_RSS' = right_RSS, 
                    'Q' = total_RSS, 'C_left' = left_mod$C, 'b_left' = left_mod$b,
                    'C_right' = right_mod$C, 'b_right' = right_mod$b)
    if(i == 3){oldQ = c(test[[i-2]][3],1)}
    
    mod_dat$predictsL = fcurve(mod_dat$x, left_mod$C, left_mod$b)
    mod_dat$predictsR = fcurve(mod_dat$x, right_mod$C, right_mod$b)
    p[[i-2]] = ggplot(mod_dat, aes(x, sim_flows)) + geom_point() +
        geom_line(data = subset(mod_dat, x <= cp), aes(y = predictsL), size = 1, color = 'red') +
        geom_line(data = subset(mod_dat, x >= cp), aes(y = predictsR), size = 1, color = 'blue') +
        labs(x = '', y = 'Simulated Flows', title = 'Fluid Flow Changepoint Simulation')
    # in this list p[[cp - 1]] will be the ideal plot
    names(test)[i-2] = paste('misspecified at point x =', i)
    names(p)[i-2] = paste('misspecified at point x =', i)
    
    if(oldQ[1] > test[[i-2]][3]){oldQ = c(test[[i-2]][3], i-2)}
  }
  
  does_CP_minimize = oldQ[2] == cp-2
  names(test)[cp-2] = 'correct'
  names(p)[cp-2] = 'correct'
  
  return(list('vals' = test, 'does_CP_minimize' = does_CP_minimize, 'plots' = p))
}
#Jacob Townson

# Simulating Determinant of Matrix of 2nd Derivatives

## Here we are simulating data to observe the behavior of the determinant of the 
## matrix of 2nd partial derivatives of Q. The determinant of partials is as follows:

## det(M) = 16 \left( \sum_{i = 1}^{j} \sum_{k = 1}^{j}(x_i ^{b_1}x_k^{b_1}\ln(x_k)
##          (y_k-2C_1x_k^{b_1}))(C_1x_i^{b_1}\ln(x_k) + \ln(x_i)(y_i-2C_1x_i^{b_1})) 
##          \right) * \left( \sum_{i = j+1}^{n} \sum_{k = j+1}^{n}(x_i ^{b_2}x_k^{b_2}
##          \ln(x_k)(y_k-2C_2x_k^{b_2}))(C_2x_i^{b_2}\ln(x_k) + \ln(x_i)(y_i-2C_2x_i^{b_2})) \right)
## This is written in latex to make it a little easier to read if copied to a latex editor. 

## So we only need to test simulated data on x_i^{b_1}x_k^{b_1}\ln(x_k)(y_k-2C_1x_k^{b_1}) and
## C_1x_i^{b_1}\ln(x_k) + \ln(x_i)(y_i-2C_1x_i^{b_1})
## If either of these expressions (with their sums) ever equal 0, then the entire determinant will also equal 0.
## So let's start here, sepcifically looking at 
## \left(\sum_{i = 1}^{j}\sum_{k = 1}^{j}(x_i^{b_1}x_k^{b_1}\ln(x_k)(y_k-2C_1x_k^{b_1}))(C_1x_i^{b_1}\ln(x_k)+\ln(x_i)(y_i-2C_1x_i^{b_1}))\right)
## So we won't need to test anything with the cp for this simulation code.

set.seed(82394)
source("./functions/simulate_model_data.R")
## n = Number of x values, N = number of observations at each x,
## C = The constant value in our flow function, b = exponent value in flow function, 
## cp = changepoint, sd = standard deviation of the simulated data



# Standard Example

n = 20
N = 10
C = 10
b = .3
sd = 1

sim_datasets = list()
determ = c()
t=1
for(t in 1:100){
  sim_datasets[[t]] = no_cp_sim(N,n,C,b,sd)
  
  # Couldn't get this to work, but would be faster than for loops
  #double_sum <- function(i, y_i){
  #  sum(sapply(1:n,function(k, y_k)sum((i^b)*(k^b)*log(k)(y_k-2*C*(k^b))*(C*(i^b)*log(k)+log(i)(y_i-2*C*(i^b))))))
  #}
  #determ[t] = sum(sapply(1:n,double_sum))
  i = 1
  k = 1
  determ[t] = 0
  current_sim = sim_datasets[[t]]
  rows = dim(current_sim)[1]
  for(i in 1:rows){
    for(k in 1:rows){
      determ[t] = determ[t] + 
                  -4*((current_sim[i,1]^b)*(current_sim[k,1]^b)*log(current_sim[k,1])*(current_sim[k,2]-2*C*(current_sim[k,1]^b)))*(C*(current_sim[i,1]^b)*log(current_sim[k,1])+log(current_sim[i,1])*(current_sim[i,2]-2*C*(current_sim[i,1]^b)))
      k = k + 1
    }
    i = i + 1
  }
  print(paste("The determinant of simulation", t, "was", determ[t]))
  t = t+1
}

summary(determ)





# x values are all between 0 and 1 Example

n = 20
N = 5
C = 10
b = .9
sd = 1

sim_datasets2 = list()
determ2 = c()
t=1
for(t in 1:100){
  sim_datasets2[[t]] = no_cp_sim_zero2one(N,n,C,b,sd)
  i = 1
  k = 1
  determ2[t] = 0
  current_sim = sim_datasets2[[t]]
  rows = dim(current_sim)[1]
  for(i in 1:rows){
    for(k in 1:rows){
      determ2[t] = determ2[t] + 
        -4*((current_sim[i,1]^b)*(current_sim[k,1]^b)*log(current_sim[k,1])*(current_sim[k,2]-2*C*(current_sim[k,1]^b)))*(C*(current_sim[i,1]^b)*log(current_sim[k,1])+log(current_sim[i,1])*(current_sim[i,2]-2*C*(current_sim[i,1]^b)))
      k = k + 1
    }
    i = i + 1
  }
  print(paste("The determinant of simulation", t, "was", determ2[t]))
  t = t+1
}
#plot(sim_datasets2[[1]]$x, sim_datasets2[[1]]$sim_flows)

summary(determ2)























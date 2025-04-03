#Jacob Townson

# Testing Derivatives of Q

## Here I am testing whether or not the derivatives of Q hold up in our simulations.
## When speaking of the derivatives of Q, I am referencing the notes from 01-14-2022
## contained in the "Notes" folder in the repo. 

source('./consistency_sims.R')
set.seed(82394)
N = 1
n = 20
x = sort(rep(1:n, N))
C = 10
b = .25
delta = 1.25
sd = 0
cp = 10
f=C*x^b*delta* ifelse(x>cp, 1, 1/delta)
sim_flows=f+rnorm(N*n,sd=sd)

# Since sd = 0, we know that there will be no random element to these simulations.
# Thus the above sim_flows values will be the same as below in the tests.

test17 = consistency_tester(N, n, C, b, delta, 
                            cp, sd = 0)
x1 = x[x <= cp]
x2 = x[x > cp]
sim_flows1 = sim_flows[x <= cp]
sim_flows2 = sim_flows[x > cp]

# CORRECT CP
#C1
lhs = unname(test17$vals$correct[4]*sum(x1^(2*test17$vals$correct[5])))
rhs = sum(sim_flows1*x1^(test17$vals$correct[5]))
lhs-rhs  #Notice the amount here is negligible. It's 0, as expected

#C2
lhs = unname(test17$vals$correct[6]*sum(x2^(2*test17$vals$correct[7])))
rhs = sum(sim_flows2*x2^(test17$vals$correct[7]))
lhs-rhs  #Notice the amount here is negligible. It's 0, as expected

#b1
lhs = unname(test17$vals$correct[4]*sum((log(x1))*x1^(2*test17$vals$correct[5])))
rhs = sum((log(x1))*sim_flows1*x1^(test17$vals$correct[5]))
lhs-rhs  #It's 0, as expected

#b2
lhs = unname(test17$vals$correct[6]*sum((log(x2))*x2^(2*test17$vals$correct[7])))
rhs = sum((log(x2))*sim_flows2*x2^(test17$vals$correct[7]))
lhs-rhs  #It's 0, as expected




# INCORRECT CP (LEFT)
#C1
lhs = unname(test17$vals$misspecified_left[4]*sum(x1^(2*test17$vals$misspecified_left[5])))
rhs = sum(sim_flows1*x1^(test17$vals$misspecified_left[5]))
lhs-rhs  #Notice the amount here is negligible. It's 0, as expected

#C2
lhs = unname(test17$vals$misspecified_left[6]*sum(x2^(2*test17$vals$misspecified_left[7])))
rhs = sum(sim_flows2*x2^(test17$vals$misspecified_left[7]))
lhs-rhs  #It's non zero because the right of the cp was not correctly specified

#b1
lhs = unname(test17$vals$misspecified_left[4]*sum((log(x1))*x1^(2*test17$vals$misspecified_left[5])))
rhs = sum((log(x1))*sim_flows1*x1^(test17$vals$misspecified_left[5]))
lhs-rhs  #It's 0, as expected since the left hand side contained no incorrect specifications

#b2
lhs = unname(test17$vals$misspecified_left[6]*sum((log(x2))*x2^(2*test17$vals$misspecified_left[7])))
rhs = sum((log(x2))*sim_flows2*x2^(test17$vals$misspecified_left[7]))
lhs-rhs  #It's non zero because the right of the cp was not correctly specified




# INCORRECT CP (RIGHT)
#C1
lhs = unname(test17$vals$misspecified_right[4]*sum(x1^(2*test17$vals$misspecified_right[5])))
rhs = sum(sim_flows1*x1^(test17$vals$misspecified_right[5]))
lhs-rhs  #It's non zero because the left of the cp was not correctly specified

#C2
lhs = unname(test17$vals$misspecified_right[6]*sum(x2^(2*test17$vals$misspecified_right[7])))
rhs = sum(sim_flows2*x2^(test17$vals$misspecified_right[7]))
lhs-rhs  #Notice the amount here is negligible. It's 0, as expected

#b1
lhs = unname(test17$vals$misspecified_right[4]*sum((log(x1))*x1^(2*test17$vals$misspecified_right[5])))
rhs = sum((log(x1))*sim_flows1*x1^(test17$vals$misspecified_right[5]))
lhs-rhs  #It's non zero because the left of the cp was not correctly specified

#b2
lhs = unname(test17$vals$misspecified_right[6]*sum((log(x2))*x2^(2*test17$vals$misspecified_right[7])))
rhs = sum((log(x2))*sim_flows2*x2^(test17$vals$misspecified_right[7]))
lhs-rhs  #It's 0, as expected since the right hand side contained no incorrect specifications
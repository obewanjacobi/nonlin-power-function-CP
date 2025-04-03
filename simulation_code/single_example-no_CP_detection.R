# create simulated data
# f(C,b) = C(p_1 ^b_1 + ... + p_n^b_n)
require(plyr)
require(dplyr)
require(minpack.lm)
require(ggplot2)
fcurve <- function(MW = 1:20, C, b){
  indiv_flow = C*(as.numeric(MW)^b)
  return(indiv_flow)
}

# Simulate the data

#set.seed(82394)
Cdelta = 1.1
bdelta = 1.1
R = 500
N = 20
J= 20
C = 10
b = .3
cp = 6.5
sd = 2
x = sort(rep(1:J,N))
guess = c(-Inf,Inf)

ne = rep(0,R)
f=C*x^(b * ifelse(x < cp, 1, bdelta)) * Cdelta * ifelse(x>cp, 1, 1/Cdelta)
sim_flows=f+rnorm(N*J,sd=sd)
mod_dat = data.frame(x, sim_flows)
lhs = filter(mod_dat, x < cp)
rhs = filter(mod_dat, x > cp)

# Build a model on the simulated data

left_mod = nlsLM(formula = sim_flows ~ C*x^b, 
                 data = lhs, start = list(C = 5, b=.5),
                 upper = c(20, 1), lower = c(0, 0))

right_mod = nlsLM(formula = sim_flows ~ C*x^b, 
                  data = rhs, start = list(C = 5, b=.5),
                  upper = c(20, 1), lower = c(0, 0))

tot_mod = nlsLM(formula = sim_flows ~ C*x^b, 
                data = mod_dat, start = list(C = 5, b=.5),
                upper = c(20, 1), lower = c(0, 0))
## Note- pline_mod$control$ftol is how this model determined when to terminate its iterations:
## Termination occurs when both the actual and predicted relative reductions in the sum of squares are at most ftol.
## Also note that nlsLM() has a subset parameter, however it is unclear on how to use it, thus it will be avoided.

# Visualize the model on the simulated data

mod_dat$predictsL = fcurve(mod_dat$x, coef(left_mod)[[1]], coef(left_mod)[[2]])
mod_dat$predictsR = fcurve(mod_dat$x, coef(right_mod)[[1]], coef(right_mod)[[2]])
mod_dat$predicts  = fcurve(mod_dat$x, coef(tot_mod)[[1]], coef(tot_mod)[[2]])
p = ggplot(mod_dat, aes(x, sim_flows)) + geom_point() +
    geom_line(data = subset(mod_dat, x <= cp), aes(y = predictsL), size = 1.5, color = 'red') +
    geom_line(data = subset(mod_dat, x >= cp), aes(y = predictsR), size = 1.5, color = 'blue') +
    geom_line(data = mod_dat, aes(y = predicts), size = 1, color = 'darkorchid4') +
    labs(x = '', y = 'Simulated Flows', title = 'Pipeline Changepoint Simulation')
p
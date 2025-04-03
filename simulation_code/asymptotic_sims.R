# Asymptotic Normality Simulations

## This code is made to simulate our model to study the asymptotic properties of our model. 

# Setup Code

set.seed(823)
require(ggplot2)
require(plotly)
require(parallel)
library(pbapply)
source("./functions/simulate_model_data.R")

C.0=c(10,12)
b.0=c(.333,.555)
cp.0=13
sigma.0=5

C=runif(2, min = 0.1, max = C.0 + 20)
b=runif(2)
theta = c(C[1], C[2], b[1], b[2])

################# Simulation code (No CP) ################# 

# Simulation Prep

set.seed(823)
runs = 200000
results_noCP = data.frame(C_hat = numeric(0), b_hat = numeric(0))

n_cores = detectCores() - 2
cl = makeCluster(n_cores)
clusterExport(cl, c("no_cp_sim", "nonlin.NRg", "RSSfunc", "C.0", "b.0", "sigma.0"))

start_time = Sys.time()

## RECALL FOR SIMULATION FUNCTION:
## n = Number of x values, N = number of observations at each x,
## C = The constant value in our flow function, b = exponent value in flow function, 
## sd = standard deviation of the simulated data

# Run the simulations in parallel
results_list <- pblapply(1:runs, function(i) {
  sim_vals = no_cp_sim(N = 20, n = 20, C = C.0[1], b = b.0[1], sd = sigma.0)
  
  opt_sim_params = nonlin.NRg(sim_vals$x, sim_vals$sim_flows, print.output=TRUE)
  
  data.frame(C_hat = opt_sim_params$C, b_hat = opt_sim_params$b, RSS = opt_sim_params$RSS)
}, cl = cl)

stopCluster(cl)
results_noCP = do.call(rbind, results_list)

end_time = Sys.time()
execution_time = end_time - start_time
minutes = as.numeric(execution_time, units = "mins")
print(paste("No CP execution took", round(minutes, 2), "minutes"))

write.csv(results_noCP, 'simParams_noCP.csv', row.names = FALSE)

# Visualize Distribution

p_C <- ggplot(results_noCP, aes(x = C_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.25, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "C Distribution",
       x = "Values",
       y = "Density") + theme_minimal()
ggsave(filename = "./C_histogram_plot_noCP.png", 
       plot = p_C, width = 8, height = 6, dpi = 600)
p_C_plotly = ggplotly(p_C)
p_C_plotly

p_b <- ggplot(results_noCP, aes(x = b_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.005, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "b Distribution",
       x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./b_histogram_plot_noCP.png", 
       plot = p_b, width = 8, height = 6, dpi = 600)
p_b_plotly = ggplotly(p_b)
p_b_plotly
















################# Simulation code (With CP) ################# 

# Prep Simulation

set.seed(823)
runs = 50000
results_CP = data.frame(C1_hat = numeric(0), C2_hat = numeric(0), 
                         b1_hat = numeric(0), b2_hat = numeric(0),
                         tau_hat = numeric(0))

n_cores = detectCores() - 2
cl = makeCluster(n_cores)
clusterExport(cl, c("cp_sim", "nonlin.NRg_CP", "RSSfunc", "C.0", "b.0", "cp.0", "sigma.0"))

start_time = Sys.time()

## RECALL FOR SIMULATION FUNCTION:
## n = Number of x values, N = number of observations at each x,
## C1-2 = The constant value in our flow function before and after the cp,
## b1-2 = exponent value in flow function before and after the cp, 
## cp = changepoint, sd = standard deviation of the simulated data

# Run the simulations in parallel
results_list <- pblapply(1:runs, function(i) {
  sim_vals = cp_sim(N = 20, n = 20, C1 = C.0[1], b1 = b.0[1], C2 = C.0[2], b2 = b.0[2], 
                    cp=cp.0, sd = sigma.0)
  
  opt_sim_params = nonlin.NRg_CP(sim_vals$x, sim_vals$sim_flows, print.output=TRUE)
  
  data.frame(C1_hat = opt_sim_params$C1, C2_hat = opt_sim_params$C2, 
             b1_hat = opt_sim_params$b1, b2_hat = opt_sim_params$b2,
             tau_hat = opt_sim_params$changepoint)
}, cl = cl)

stopCluster(cl)
results_CP = do.call(rbind, results_list)

end_time = Sys.time()
execution_time = end_time - start_time
minutes = as.numeric(execution_time, units = "mins")
print(paste("CP execution took", round(minutes, 2), "minutes"))

write.csv(results_CP, 'simParams_withCP_normalJump.csv', row.names = FALSE)

# Visualize Distribution

p_C1 <- ggplot(results_CP, aes(x = C1_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.25, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "C1 Distribution",
       x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./C1_histogram_plot_CP.png", 
       plot = p_C1, width = 8, height = 6, dpi = 600)
p_C1_plotly = ggplotly(p_C1)
p_C1_plotly

p_C2 <- ggplot(results_CP, aes(x = C2_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.25, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "C2 Distribution",
       x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./C2_histogram_plot_CP.png", 
       plot = p_C2, width = 8, height = 6, dpi = 600)
p_C2_plotly = ggplotly(p_C2)
p_C2_plotly

p_b1 <- ggplot(results_CP, aes(x = b1_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.025, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "b1 Distribution",
       x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./b1_histogram_plot_CP.png", 
       plot = p_b1, width = 8, height = 6, dpi = 600)
p_b1_plotly = ggplotly(p_b1)
p_b1_plotly

p_b2 <- ggplot(results_CP, aes(x = b2_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.025, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "b2 Distribution",
       x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./b2_histogram_plot_CP.png", 
       plot = p_b2, width = 8, height = 6, dpi = 600)
p_b2_plotly = ggplotly(p_b2)
p_b2_plotly

p_tau <- ggplot(results_CP, aes(x = tau_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 1, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "CP Distribution",
       x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./tau_histogram_plot_CP.png", 
       plot = p_tau, width = 8, height = 6, dpi = 600)
p_tau_plotly = ggplotly(p_tau)
p_tau_plotly











################# Simulation code (With CP and very small difference after it) ################# 

# This simulation may give us a better idea of the distribution of the changepoint itself.
# While we know that the changepoint is a discrete value and is NOT continuous, we can still
# study how it behaves in our model.

# Prep Simulation

C.0=c(10,10.2)
b.0=c(.333,.355)
cp.0=9
sigma.0=5

C=runif(2, min = 0.1, max = C.0 + 20)
b=runif(2)
theta = c(C[1], C[2], b[1], b[2])

set.seed(823)
runs = 200000
results_CP = data.frame(C1_hat = numeric(0), C2_hat = numeric(0), 
                        b1_hat = numeric(0), b2_hat = numeric(0),
                        tau_hat = numeric(0))

n_cores = detectCores() - 2
cl = makeCluster(n_cores)
clusterExport(cl, c("cp_sim", "nonlin.NRg_CP", "RSSfunc", "C.0", "b.0", "cp.0", "sigma.0"))

start_time = Sys.time()

## RECALL FOR SIMULATION FUNCTION:
## n = Number of x values, N = number of observations at each x,
## C1-2 = The constant value in our flow function before and after the cp,
## b1-2 = exponent value in flow function before and after the cp, 
## cp = changepoint, sd = standard deviation of the simulated data

# Run the simulations in parallel
results_list <- pblapply(1:runs, function(i) {
  sim_vals = cp_sim(N = 20, n = 20, C1 = C.0[1], b1 = b.0[1], C2 = C.0[2], b2 = b.0[2], 
                    cp=cp.0, sd = sigma.0)
  
  opt_sim_params = nonlin.NRg_CP(sim_vals$x, sim_vals$sim_flows, print.output=TRUE)
  
  data.frame(C1_hat = opt_sim_params$C1, C2_hat = opt_sim_params$C2, 
             b1_hat = opt_sim_params$b1, b2_hat = opt_sim_params$b2,
             tau_hat = opt_sim_params$changepoint)
}, cl = cl)

stopCluster(cl)
results_CP = do.call(rbind, results_list)

end_time = Sys.time()
execution_time = end_time - start_time
minutes = as.numeric(execution_time, units = "mins")
print(paste("CP execution took", round(minutes, 2), "minutes"))

write.csv(results_CP, 'simParams_withCP_smallJump.csv', row.names = FALSE)

# Visualize Distribution

p_C1 <- ggplot(results_CP, aes(x = C1_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.25, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "C1 Distribution",
       x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./C1_histogram_plot_CP.png", 
       plot = p_C1, width = 8, height = 6, dpi = 600)
p_C1_plotly = ggplotly(p_C1)
p_C1_plotly

p_C2 <- ggplot(results_CP, aes(x = C2_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.25, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "C2 Distribution",
       x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./C2_histogram_plot_CP.png", 
       plot = p_C2, width = 8, height = 6, dpi = 600)
p_C2_plotly = ggplotly(p_C2)
p_C2_plotly

p_b1 <- ggplot(results_CP, aes(x = b1_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.025, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "b1 Distribution",
       x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./b1_histogram_plot_CP.png", 
       plot = p_b1, width = 8, height = 6, dpi = 600)
p_b1_plotly = ggplotly(p_b1)
p_b1_plotly

p_b2 <- ggplot(results_CP, aes(x = b2_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.025, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "b2 Distribution",
       x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./b2_histogram_plot_CP.png", 
       plot = p_b2, width = 8, height = 6, dpi = 600)
p_b2_plotly = ggplotly(p_b2)
p_b2_plotly

p_tau <- ggplot(results_CP, aes(x = tau_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 1, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(title = "CP Distribution",
       x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./tau_histogram_plot_CP.png", 
       plot = p_tau, width = 8, height = 6, dpi = 600)
p_tau_plotly = ggplotly(p_tau)
p_tau_plotly



























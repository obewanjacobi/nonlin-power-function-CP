# Asymptotic Normality Simulations

## This code is made to simulate our model to study the asymptotic properties of our model. 

# Setup Code

set.seed(823)
require(ggplot2)
require(plotly)
require(parallel)
library(pbapply)
source("./functions/simulate_model_data.R")

C.0=c(10,10.2)
b.0=c(.333,.355)
cp.0=9
sigma.0=2.5




################# Simulation code (few n's) ################# 

# Prep Simulation

set.seed(823)
runs = 20000
results_CP = data.frame(C1_hat = numeric(0), C2_hat = numeric(0), 
                         b1_hat = numeric(0), b2_hat = numeric(0),
                         tau_hat = numeric(0))

n_cores = detectCores() - 2
cl = makeCluster(n_cores)
clusterExport(cl, c("cp_sim", "nonlin.NRg_CP", "RSSfunc", "C.0", "b.0", "cp.0", "sigma.0"))

start_time = Sys.time()

# Run the simulations in parallel
results_list <- pblapply(1:runs, function(i) {
  sim_vals = cp_sim(N = 5, n = 20, C1 = C.0[1], b1 = b.0[1], C2 = C.0[2], b2 = b.0[2], 
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

write.csv(results_CP, 'simParams_withCP_few_ns.csv', row.names = FALSE)

# Visualize Distribution

p_C1 <- ggplot(results_CP, aes(x = C1_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.25, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./C1_histogram_plot_CP_few_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./C2_histogram_plot_CP_few_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./b1_histogram_plot_CP_few_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./b2_histogram_plot_CP_few_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./tau_histogram_plot_CP_few_ns.png", 
       plot = p_tau, width = 8, height = 6, dpi = 600)
p_tau_plotly = ggplotly(p_tau)
p_tau_plotly











################# Simulation code (medium number of n's) ################# 

# Prep Simulation

set.seed(823)
runs = 20000
results_CP = data.frame(C1_hat = numeric(0), C2_hat = numeric(0), 
                        b1_hat = numeric(0), b2_hat = numeric(0),
                        tau_hat = numeric(0))

n_cores = detectCores() - 2
cl = makeCluster(n_cores)
clusterExport(cl, c("cp_sim", "nonlin.NRg_CP", "RSSfunc", "C.0", "b.0", "cp.0", "sigma.0"))

start_time = Sys.time()

# Run the simulations in parallel
results_list <- pblapply(1:runs, function(i) {
  sim_vals = cp_sim(N = 50, n = 20, C1 = C.0[1], b1 = b.0[1], C2 = C.0[2], b2 = b.0[2], 
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

write.csv(results_CP, 'simParams_withCP_medium_ns.csv', row.names = FALSE)

# Visualize Distribution

p_C1 <- ggplot(results_CP, aes(x = C1_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.25, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./C1_histogram_plot_CP_medium_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./C2_histogram_plot_CP_medium_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./b1_histogram_plot_CP_medium_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./b2_histogram_plot_CP_medium_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./tau_histogram_plot_CP_medium_ns.png", 
       plot = p_tau, width = 8, height = 6, dpi = 600)
p_tau_plotly = ggplotly(p_tau)
p_tau_plotly













################# Simulation code (embiggen n's) ################# 

# Prep Simulation

set.seed(823)
runs = 20000
results_CP = data.frame(C1_hat = numeric(0), C2_hat = numeric(0), 
                        b1_hat = numeric(0), b2_hat = numeric(0),
                        tau_hat = numeric(0))

n_cores = detectCores() - 2
cl = makeCluster(n_cores)
clusterExport(cl, c("cp_sim", "nonlin.NRg_CP", "nonlin.NRg_CP_GPU", "RSSfunc", "C.0", "b.0", "cp.0", "sigma.0"))

start_time = Sys.time()

# Run the simulations in parallel
results_list <- pblapply(1:runs, function(i) {
  sim_vals = cp_sim(N = 150, n = 20, C1 = C.0[1], b1 = b.0[1], C2 = C.0[2], b2 = b.0[2], 
                    cp=cp.0, sd = sigma.0)
  
  opt_sim_params = nonlin.NRg_CP(sim_vals$x, sim_vals$sim_flows, print.output=TRUE)
  #opt_sim_params = nonlin.NRg_CP_GPU(sim_vals$x, sim_vals$sim_flows, print.output=TRUE) # GPU accelerated version
                                                                                         # Complex installation
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

write.csv(results_CP, 'simParams_withCP_LOTSO_ns.csv', row.names = FALSE)

# Visualize Distribution

p_C1 <- ggplot(results_CP, aes(x = C1_hat)) +
  geom_histogram(aes(y = after_stat(density)), 
                 binwidth = 0.25, 
                 fill = "lightblue", 
                 color = "blue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1.2) +
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./C1_histogram_plot_CP_LOTSO_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./C2_histogram_plot_CP_LOTSO_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./b1_histogram_plot_CP_LOTSO_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./b2_histogram_plot_CP_LOTSO_ns.png", 
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
  labs(x = "Values",
       y = "Density") +
  theme_minimal()
ggsave(filename = "./tau_histogram_plot_CP_LOTSO_ns.png", 
       plot = p_tau, width = 8, height = 6, dpi = 600)
p_tau_plotly = ggplotly(p_tau)
p_tau_plotly



























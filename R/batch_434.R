source("lotka_volterra_stan_functions.R")

rda_list <- list()
# for(noise_level in c(0.001,0.002,0.005, 0.007,0.01,0.02,0.05, 0.07)) {
for(noise_level in c(0.01)) {
  print(noise_level)
  data <- build_data(myseed=434, species=3, noise_level = noise_level,folder="stan_lotka")
  # Fitting
  fit <- fit_stan_rk4(data,folder="stan_lotka")
  print(fit)
  
  # Plotting
  make_stan_plots(fit,data,folder="stan_lotka")
  post <- plot_analysis(fit,data,folder="stan_lotka")
  rda_list[[noise]] <- fit
}

# save(rda_list, file="stan_lotka_434_fits.rda")
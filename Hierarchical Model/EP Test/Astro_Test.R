require("rstan")

#  ..... Compile Stan Code .....
# astro_fit <- stan_model(model_code = astro_code, model_name = "Astro Model");

# ..... Determine Priors for the Simulation .....
fitted_a <- matrix(0, nrow = 360, ncol = 8);

for(i in 0:359){
  cat("Degrees: ",i,"\n")
  current_data <- read.csv(paste("../Exponential\ Model/Output/ParamsBy0.1/",i,"-",i+0.1,".csv",sep=""));
  fitted_a[i+1,] <- apply(current_data, 2, median)[-1];
}

rm(current_data)
mu_log_a <- apply(log(fitted_a), 2, median)
log_tau_log_a <- log(apply(log(fitted_a), 2, sd)) #confirmed works well without log at the end!

# ...... Simulate the Data ......
J = 50;
N = 100;
astro_data <- astro_sim(phi = c(mu_log_a, log_tau_log_a, 5.3), J = J, n.samp = N);
plot(astro_data$x, astro_data$y)

#  ..... Priors .....
pri_mu <- c(mu_log_a, log_tau_log_a, 5.3);
pri_sig <- diag(17) * 0.5;
init_data <- list();
for(i in 1:4) { init_data[[i]] <- list(phi = pri_mu, eta = matrix(0, nrow = 8, ncol = J));}

# ...... EP Fit .......
# 1300 seconds with K = 5
# 1350 seconds with K = 10
# 1170 seconds with K = 25
astro_EP_fit_5 <- EP(fit = astro_fit, J = J, K = 5, data = astro_data, S = 2, mc_iter = 100, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
astro_EP_fit_10 <- EP(fit = astro_fit, J = J, K = 10, data = astro_data, S = 2, mc_iter = 100, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
astro_EP_fit_25 <- EP(fit = astro_fit, J = J, K = 25, data = astro_data, S = 2, mc_iter = 100, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
astro_EP_fit_25

# ...... MC Fit ....... 
# 908 seconds
data <- list(N = N * J, B = J, x = astro_data$x, y = astro_data$y, bin = astro_data$bin, Mu_Cav = pri_mu, Sig_Cav = pri_sig);
system.time(astro_MC_fit <- sampling(astro_fit, data = data, iter = 100, chains = 4, init_data = init_data));
colMeans(extract(astro_MC_fit)$phi)
apply(extract(astro_MC_fit)$phi, 2, var)


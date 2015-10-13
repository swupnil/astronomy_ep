require("rstan")

#  ..... Compile Stan Code .....
model4_fit <- stan_model(model_code = model4_code, model_name = "Model 4");

# ..... Determine Priors for the Simulation .....
mu_log_a <- log(c(361, 10, 2, 650))
log_tau_log_a <- log(c(0.34, 0.5, 0.5, .55))

# ...... Simulate the Data ......
J = 50;
N = 500;
model4_data <- model4_sim(phi = c(mu_log_a, log_tau_log_a, 2, .0001), J = J, n.samp = N);
plot(model4_data$x, model4_data$y)

#  ..... Priors .....
pri_mu <- c(mu_log_a, log_tau_log_a, 2, .0001);
pri_sig <- diag(length(pri_mu)) * 0.5;
init_data <- list();
for(i in 1:4) { init_data[[i]] <- list(phi = pri_mu, eta = matrix(0, nrow = 5, ncol = J));}

# ...... EP Fit .......
# 1089 seconds with K = 5
# -- seconds with K = 10
# -- seconds with K = 25
model4_EP_fit_5 <- EP(fit = model4_fit, J = J, K = 5, data = model4_data, S = 2, mc_iter = 100, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
model4_EP_fit_10 <- EP(fit = model4_fit, J = J, K = 10, data = model4_data, S = 2, mc_iter = 100, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
model4_EP_fit_25 <- EP(fit = model4_fit, J = J, K = 25, data = model4_data, S = 2, mc_iter = 100, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
model4_EP_fit_5

# ...... MC Fit ....... 
# 1350 seconds
data <- list(N = N * J, B = J, x = model4_data$x, y = model4_data$y, bin = model4_data$bin, Mu_Cav = pri_mu, Sig_Cav = pri_sig);
system.time(model4_MC_fit <- sampling(model4_fit, data = data, iter = 100, chains = 4, init_data = init_data));
colMeans(extract(model4_MC_fit)$phi)
apply(extract(model4_MC_fit)$phi, 2, sd)


require("rstan")

#  ..... Compile Stan Code .....
model0_fit <- stan_model(model_code = model0_code, model_name = "Model 0");

# ...... Simulate the Data ......
sub_data0 <- model0_sim(phi = 200, n.samp = 5000);

#  ..... Priors .....
pri_mu_cav <- 0;
pri_sig_cav <- 20^2;
init_data <- list();
for(i in 1:4) { init_data[[i]] <- list(phi = pri_mu_cav);}

# ...... EP Fit .......
model0_EP_fit <- Model0_EP(fit = model0_fit, data = sub_data0, S = 3, K = 5, prior.Mu = pri_mu_cav, prior.Sigma.inv = 1/pri_sig_cav);
model0_EP_fit

# ...... MC Fit .......
data <- list(N = length(sub_data0$y), y = sub_data0$y, Mu_Cav = pri_mu_cav, Sig_Cav = pri_sig_cav);
system.time(model0_MC_fit <- sampling(model0_fit, data = data, iter = 100, chains = 4, init_data = init_data));
mean(extract(model0_MC_fit)$phi)
var(extract(model0_MC_fit)$phi)

require("rstan")

#  ..... Compile Stan Code .....
easy_fit <- stan_model(model_code = easy_code, model_name = "Easy Fit");

#  ..... Priors .....
pri_mu_cav <- c(200, 20, 1, 0.1)
init_data <- list();
for(i in 1:4) { init_data[[i]] <- list(eta = matrix(0, nrow=2, ncol=J), phi = pri_mu_cav);}

# ...... Simulate the Data ......
easy_J <- 72;
easy_sub <- easy.sim(phi = c(200, 20, 1, 0.1), J = easy_J, n.samp = 200);

# ...... EP Fit .......
easy_EP_fit <- astro.EP(fit = easy_fit, data = easy_sub, J = easy_J, K = easy_J, S = 2, prior.Mu = pri_mu_cav, SMOOTH = 0.3, parallel = TRUE);
easy_EP_fit$Post.Mu

# ...... MC Fit .......
easy_data <- list(N = length(easy_sub$y), M = 4, B = easy_J, y = easy_sub$y, bin = easy_sub$bin, Mu_Cav = pri_mu_cav, Sig_Cav = diag(4) * 4);
system.time(easy_MC_fit <- sampling(easy_fit, data = easy_data, iter = 100, chains = 4, init = init_data));
colMeans(extract(easy_MC_fit)$phi)

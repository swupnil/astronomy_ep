require("rstan")

#  ..... Compile Stan Code .....
model2_fit <- stan_model(model_code = model2_code, model_name = "Model 2");

# ...... Simulate the Data ......
J = 50;
N = 1000;
sub_data2 <- model2_sim(phi = c(200, log(5)), J = J, n.samp = N);

#  ..... Priors .....
pri_mu <- c(0, 0);
pri_sig <- matrix(c(1000, 0, 0, 1), nrow = 2, ncol = 2);
init_data <- list();
for(i in 1:4) { init_data[[i]] <- list(phi = pri_mu, theta = rep(0, J));}

# ...... EP Fit .......
# 452.5 seconds with K = 5
# 467.3 seconds with K = 10
# 296 seconds with K = 25
model2_EP_fit <- EP(fit = model2_fit, J = J, K = 25, data = sub_data2, S = 1, mc_iter = 1000, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
model2_EP_fit

# ...... MC Fit ....... 
# 1137 seconds
data <- list(N = N*J, B = J, y = sub_data2$y, bin = sub_data2$bin, Mu_Cav = pri_mu, Sig_Cav = pri_sig);
system.time(model2_MC_fit <- sampling(model2_fit, data = data, iter = 1000, chains = 4, init_data = init_data));
model2_MC_fit
mean(extract(model2_MC_fit)$phi[,1])
var(extract(model2_MC_fit)$phi[,1])
mean(extract(model2_MC_fit)$phi[,2])
var(extract(model2_MC_fit)$phi[,2])
cov(extract(model2_MC_fit)$phi[,2], extract(model2_MC_fit)$phi[,1])

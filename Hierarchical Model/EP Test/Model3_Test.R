require("rstan")

#  ..... Compile Stan Code .....
model3_fit <- stan_model(model_code = model3_code, model_name = "Model 3");

# ...... Simulate the Data ......
J = 50;
N = 100;
sub_data3 <- model3_sim(phi = c(-100, 10, log(5), log(4), log(20)), J = J, n.samp = N);

#  ..... Priors .....
pri_mu <- c(0, 0, 0, 0, 0);
pri_sig <- diag(5);
pri_sig[1,1] <- pri_sig[2,2] <- 1000;
pri_sig[3,3] <- pri_sig[4,4] <- pri_sig[5,5] <- 2;
init_data <- list();
for(i in 1:4) { init_data[[i]] <- list(phi = pri_mu, beta0 = rep(0, J), beta1 = rep(0,J));}

# ...... EP Fit .......
# 61 seconds with K = 5
# 48 seconds with K = 10
# 60 seconds with K = 25
model3_EP_fit_5 <- EP(fit = model3_fit, J = J, K = 5, data = sub_data3, S = 3, mc_iter = 1000, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
model3_EP_fit_10 <- EP(fit = model3_fit, J = J, K = 10, data = sub_data3, S = 3, mc_iter = 1000, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
model3_EP_fit_25 <- EP(fit = model3_fit, J = J, K = 25, data = sub_data3, S = 2, mc_iter = 1000, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
model3_EP_fit_10

# ...... MC Fit ....... 
# 81 seconds
data <- list(N = N * J, B = J, x = sub_data3$x, y = sub_data3$y, bin = sub_data3$bin, Mu_Cav = pri_mu, Sig_Cav = pri_sig);
system.time(model3_MC_fit <- sampling(model3_fit, data = data, iter = 1000, chains = 4, init_data = init_data));
model3_MC_fit
mean(extract(model3_MC_fit)$phi[,1])
var(extract(model3_MC_fit)$phi[,1])
mean(extract(model3_MC_fit)$phi[,2])
var(extract(model3_MC_fit)$phi[,2])
mean(extract(model3_MC_fit)$phi[,3])
var(extract(model3_MC_fit)$phi[,3])
mean(extract(model3_MC_fit)$phi[,4])
var(extract(model3_MC_fit)$phi[,4])
mean(extract(model3_MC_fit)$phi[,5])
var(extract(model3_MC_fit)$phi[,5])

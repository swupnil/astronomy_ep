require("rstan")

#  ..... Compile Stan Code .....
model3_fit <- stan_model(model_code = model3_code, model_name = "Model 3");

# ...... Simulate the Data ......
J = 50;
N = 100;
sub_data3 <- model3_sim(phi = c(100, 20, log(20), log(1), log(1.5), log(.5)), J = J, n.samp = N);
plot(sub_data3$x, sub_data3$y)

#  ..... Priors .....
pri_mu <- rep(0, 6);
pri_sig <- diag(c(rep(1000,2), 2, rep(.5,3)));
init_data <- list();
for(i in 1:4) { init_data[[i]] <- list(phi = pri_mu, eta = matrix(0, nrow = 3, ncol = J));}

# ...... EP Fit .......
# 680 seconds with K = 5
# 680 seconds with K = 10
model3_EP_fit_5 <- EP(fit = model3_fit, J = J, K = 5, data = sub_data3, S = 3, mc_iter = 1000, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
model3_EP_fit_10 <- EP(fit = model3_fit, J = J, K = 10, data = sub_data3, S = 2, mc_iter = 1000, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
model3_EP_fit_25 <- EP(fit = model3_fit, J = J, K = 25, data = sub_data3, S = 2, mc_iter = 1000, prior.Mu = pri_mu, prior.Sigma.inv = solve(pri_sig));
model3_EP_fit_25

# ...... MC Fit ....... 
# 535 seconds
data <- list(N = N * J, B = J, x = sub_data3$x, y = sub_data3$y, bin = sub_data3$bin, Mu_Cav = pri_mu, Sig_Cav = pri_sig);
system.time(model3_MC_fit <- sampling(model3_fit, data = data, iter = 1000, chains = 4, init = init_data));
colMeans(extract(model3_MC_fit)$phi)
apply(extract(model3_MC_fit)$phi, 2, sd)

# no log on beta0 and beta1 works with strong prior
# no log on beta0 and beta1 works with wide enough zero-centered prior
# log on beta0 and beta1 works ok with strong prior but not great
# log on beta0 and beta1 doesn't work with flat prior
require("rstan")

#  ..... Compile Stan Code .....
model1_fit <- stan_model(model_code = model1_code, model_name = "Model 1");

# ...... Simulate the Data ......
J = 50;
N = 1000;
sub_data1 <- model1_sim(phi = 200, J = J, n.samp = N);

#  ..... Priors .....
pri_mu <- 0;
pri_sig <- 20^2;
init_data <- list();
for(i in 1:4) { init_data[[i]] <- list(phi = pri_mu, theta = rep(0, J));}

# ...... EP Fit ....... 
#81 seconds with K = 5
#73 seconds with K = 10
#69 seconds with K = 25
model1_EP_fit <- EP_1D(fit = model1_fit, J = J, K = 25, data = sub_data1, S = 3, prior.Mu = pri_mu, prior.Sigma.inv = 1/pri_sig, mc_iter = 1000);
model1_EP_fit

# ...... MC Fit ....... 
#201.0 seconds
#183.0 seconds; strong prior
data <- list(N = N*J, B = J, y = sub_data1$y, bin = sub_data1$bin, Mu_Cav = pri_mu, Sig_Cav = pri_sig);
system.time(model1_MC_fit <- sampling(model1_fit, data = data, iter = 1000, chains = 4, init_data = init_data));
model1_MC_fit
mean(extract(model1_MC_fit)$phi)
var(extract(model1_MC_fit)$phi)

load("Exponential.RData")
require("rstan")

#  ..... Compile Stan Code .....
fit <- stan_model(model_code = tilt_code, model_name = "Tilt Fit");

# .... Take Subset of Data ....
Max <- 360;
subset <- which(long < Max & long > 0 & i100 > 0 & i100 < 4);
J <- 72;
K <- 9;
n.samp <- 200;
bin <- ceiling(long[subset]*J/Max-1)%%J+1;

# .... Stratified sampling to sample equally from each bin ....
sub <- as.vector(apply(as.matrix(unique(bin)), 1, function(x) sample(which(bin==x), n.samp)));
bin <- bin[sub]; 
x <- i100[subset][sub]; 
y <- FUV[subset][sub];
sub.data <- data.frame(x = x, y = y, bin = bin);

#  ..... Priors .....
pri.mu.cav <- log(c(200, 2.5, 2.5, 250, 5, .5, 50, 7, 1.6, 1.5, 1.6, 1.6, 1.3, 1.3, 1.15, 1.2, 200))
init.data <- list();
for(i in 1:4) { init.data[[i]] <- list(eta = matrix(0, nrow=8, ncol=J), phi = pri.mu.cav);}

#  ..... MCMC: Actual Data .....
actual.data <- list(N = length(y), M = 17, B = J, x = sub.data$x, y = sub.data$y, bin = sub.data$bin, Mu_Cav = pri.mu.cav, Sig_Cav = diag(17) * 4);
system.time(actual.fit <- sampling(fit, data = actual.data, iter = 100, chains = 4, init = init.data)); 
colMeans(extract(actual.fit)$phi)
#159 seconds, N = 1200, N/J = 20, J = 60, Max = 60
#2676 seconds, N = 7200, N/J = 20, J = 360, Max = 360
#2536 seconds, N = 12000, N/J = 200, J = 60, Max = 60
#5113 seconds, N = 72000, N/J = 200, J = 360, Max = 360
#2007 seconds, N = 14400, N/J = 200, J = 72, Max = 360

#  ..... EP: Actual Data .....
astro.fit <- astro.EP(fit = fit, data = sub.data, J = J, K = K, S = 4, prior.Mu = pri.mu.cav, SMOOTH = 0.3, parallel = TRUE);
astro.fit$Post.Mu
#2639 seconds per S, N = 7200, N/J = 20, K = 6, J = 360, Max = 360, S = 2
#21384 seconds per S, N = 72000, N/J = 200, K = 6, J = 360, Max = 360, S = 2
#10660 seconds, N = 14400, N/J = 200, K = 72, J = 72, Max = 360, S = 4

save.image("Exponential.RData")
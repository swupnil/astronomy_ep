#script to import astronomy data and run EP on it

require("FITSio");
require("rstan");

# ..... Read in the Data ......
setwd("~/Dropbox/Documents/School/Columbia/Research/Astronomy");
data <- readFITS("allskyfuv_all.fits");
data2 <- readFITS("stars_fuv_glon.fits");
#load hierarichical_data.RData

# ..... Parse/Clean the Variables ......
FUV <- data$col[[1]];
i100 <- data$col[[2]];
i100 <- i100[-which(is.na(FUV))];
long <- data$col[[4]];
long <- long[-which(is.na(FUV))];
FUV <- FUV[-which(is.na(FUV))];

# .... Take Subset of Data ....
J <- 3600;
K <- 6;
n.samp <- 20;
subset <- which(long < 360 & long > 0 & i100 > 0 & i100 < 4);
bin <- ceiling(long[subset]*J/360-1)%%J+1;
# .... Stratified sampling to sample equally from each bin ....
sub <- as.vector(apply(as.matrix(unique(bin)), 1, function(x) sample(which(bin==x), n.samp)));
bin <- bin[sub]; 
x <- i100[subset][sub]; 
y <- FUV[subset][sub];
sub.data <- data.frame(x = x, y = y, bin = bin);
sub.data2 <- data.frame(x = x, y = y, bin = ceiling(bin/10));

#  ..... Compile the Stan Code .....
fit <- stan_model(model_code = tilt_code, model_name = "Tilt Fit");

#  ..... MCMC: Actual Data .....
actual.data <- list(N=length(y), M=4, B=J, x=sub.data$x, y=sub.data$y, bin=sub.data$bin, Mu_Cav=rep(5.5,4), Sig_Cav=diag(4)*4);
system.time(actual.fit <- sampling(fit, data = actual.data, iter = 1000, chains = 4)); 
colMeans(extract(actual.fit)$phi) #2577 seconds, N = 72000, J = 3600

actual.data2 <- list(N=length(y), M=4, B=J/10, x=sub.data$x, y=sub.data$y, bin=sub.data2$bin, Mu_Cav=rep(5.5,4), Sig_Cav=diag(4)*4);
system.time(actual.fit2 <- sampling(fit, data = actual.data2, iter = 1000, chains = 4)); 
colMeans(extract(actual.fit2)$phi) #2055 seconds, N = 72000, J = 360

#  ..... EP: Simulated Data .....
sim.data <- astro.sim(phi = colMeans(extract(actual.fit)$phi), J = J, n.samp = n.samp);

sim.fit <- astro.EP(fit = fit, data = sim.data, J = J, K = K, S = 4, prior.Mu = rep(5.5,4), SMOOTH = 0.3, parallel = TRUE);
sim.fit$Post.Mu #6646 seconds, N = 72000, J = 3600

sim.data2 <- data.frame(x = sim.data$x, y = sim.data$y, bin = ceiling(sim.data$bin/10));
sim.fit2 <- astro.EP(fit = fit, data = sim.data2, J = J/10, K = K, S = 4, prior.Mu = rep(5.5,4), SMOOTH = 0.3, paralle = TRUE); 
sim.fit2$Post.Mu #3513 seconds, N = 72000, J = 360

sim.fit4 <- astro.EP(fit = fit, data = sim.data, J = J, K = K, S = 4, prior.Mu = rep(5.5,4), SMOOTH = 0.3, random = TRUE);
sim.fit4$Post.Mu

#  ..... EP: Actual Data .....
astro.fit <- astro.EP(fit = fit, data = sub.data, J = J, K = K, S = 2, prior.Mu = rep(5.5,4), SMOOTH = 0.3, parallel = TRUE);
astro.fit$Post.Mu

astro.fit2 <- astro.EP(fit = fit, data = sub.data2, J = J/10, K = K, S = 5, prior.Mu = rep(5.5,4), SMOOTH = 0.3);
astro.fit2$Post.Mu

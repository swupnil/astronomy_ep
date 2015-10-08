require("rstan")
require("scales")

# compile code in stan for exp models
exp.fit.normal <- stan_model(model_code = exp_code_normal, model_name = "Exp Fit");
exp.fit.normal2 <- stan_model(model_code = exp_code_normal2, model_name = "Exp Fit 2");

####################
#### EXP MODEL #####
####################

par(mfrow=c(2,3))
theta <- c(1:14)*25

#normal prior
test12 <- list()
pri.mu <- log(c(200, 250, 5, .5, 50, 7, 1, 50, .5))
pri.sig <- log(c(1.6, 1.6, 1.3, 1.3, 1.15, 1.2, 1.4, 1.15, 1.4))
par(mfrow=c(2,3))
for(i in 1:14){
  test12[[i]] <- exp.fitter(degrees = theta[i], inc = .1, fit = exp.fit.normal, mu = pri.mu, sig = pri.sig)
}

#######################
#### EXP MODEL 2 ######
#######################

# normal prior
test21 <- test22 <- list()
pri.mu <- log(c(200, 2.5, 2.5, 250, 5, .5, 50, 7))
pri.sig <- log(c(1.6, 1.5, 1.6, 1.6, 1.3, 1.3, 1.15, 1.2))
par(mfrow=c(2,3))

# 1 Degree Increments
for(i in 1:14){
  test21[[i]] <- exp.fitter(degrees = theta[i], inc = 1, fit = exp.fit.normal2, mu = pri.mu, sig = pri.sig)
}

# 0.1 Degree Increments
for(i in 1:14){
  test22[[i]] <- exp.fitter(degrees = theta[i], inc = .1, fit = exp.fit.normal2, mu = pri.mu, sig = pri.sig)
}

###########################
#### EXP MODEL 2 ALL ######
###########################
test31 <- list()

for(i in 0:359){
  cat("Degrees: ",i,"\n")
  if(i%%6==0){
    path <- file.path("PlotsBy0.1/",paste(i,"-",i+5,".png",sep=""));
    png(filename = path, width = 1500, height = 1000)
    par(mfrow=c(2,3))
  }
  test31[[i+1]] <- exp.fitter(degrees = i, inc = .1, fit = exp.fit.normal2, mu = pri.mu, sig = pri.sig)             
  a <- extract(test31[[i+1]])$a
  colnames(a) <- paste(rep("a",8),as.character(0:7),sep="")
  write.csv(a, paste("ParamsBy0.1/",i,"-",i+0.1,".csv",sep=""))
  if(i%%6==5){
    dev.off()
  }
}

for(i in 0:359){
  cat("Degrees: ",i,"\n")
  if(i%%6==0){
    path <- file.path("PlotsBy1/",paste(i,"-",i+5,".png",sep=""));
    png(filename = path, width = 1500, height = 1000)
    par(mfrow=c(2,3))
  }
  temp <- exp.fitter(degrees = i, inc = 1, fit = exp.fit.normal2, mu = pri.mu, sig = pri.sig)
  a <- extract(temp)$a
  colnames(a) <- paste(rep("a",8),as.character(0:7),sep="")
  write.csv(a, paste("ParamsBy1/",i,"-",i+1,".csv",sep=""))
  if(i%%6==5){
    dev.off()
  }
}

############################
#### DETERMINE PRIORS ######
############################
fitted_a <- matrix(0, nrow = 360, ncol = 8);

for(i in 0:359){
  cat("Degrees: ",i,"\n")
  current_data <- read.csv(paste("Output/ParamsBy0.1/",i,"-",i+0.1,".csv",sep=""));
  fitted_a[i+1,] <- apply(current_data,2,median)[-1];
}

rm(current_data)
colMeans(log(fitted_a))
apply(log(fitted_a), 2, sd)

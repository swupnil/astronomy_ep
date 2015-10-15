# Model 0: Uknown Mean; Known Variance.
model0_sim <- function(phi = 200, n.samp = 1000){
  y <- rnorm(n.samp, phi, 1);
  data.frame(y = y);
}

# Model 1: Uknown Global Mean; Known Global Variance.
model1_sim <- function(phi = 200, J = 50, n.samp = 1000){
  y <- bin <- c();
  theta <- rnorm(J, phi, 10);
  
  for(j in 1:J){
    y_j <- rnorm(n.samp, theta[j], 1);
    y <- c(y, y_j);
    bin <- c(bin, rep(j, n.samp));
  }
  x <- y;
  return(data.frame(x = x, y = y, bin = bin));
}

# Model 2: Unknown Global Mean; Unknown Global Variance.
model2_sim <- function(phi = c(200, log(5)), J = 50, n.samp = 1000){
  y <- bin <- c();
  theta <- rnorm(J, phi[1], exp(phi[2]));
  
  for(j in 1:J){
    y_j <- rnorm(n.samp, theta[j], 1);
    y <- c(y, y_j);
    bin <- c(bin, rep(j, n.samp));
  }
  x <- y;
  return(data.frame(x = x, y = y, bin = bin));
}

# Model 3: Unknown Global Mean; Unknown Global Variance. Two Parameters.
model3_sim <- function(phi = c(200, -100, log(5), log(4), log(20)), J = 50, n.samp = 1000){
  x <- y <- bin <- c();
  beta0 <- rnorm(J, phi[1], exp(phi[4]));
  beta1 <- rnorm(J, phi[2], exp(phi[5]));
  log_sd <- rnorm(J, phi[3], exp(phi[6]));
  
  print(log_sd);
  print(exp(log_sd));
  for(j in 1:J){
    x_j <- runif(n.samp, 0, 50);
    y_j <- beta0[j] + beta1[j] * x_j + rnorm(n.samp, 0, exp(log_sd[j]));
    x <- c(x, x_j);
    y <- c(y, y_j);
    bin <- c(bin, rep(j, n.samp));
  }
  
  return(data.frame(x = x, y = y, bin = bin));
}
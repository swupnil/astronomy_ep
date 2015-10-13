#function to simulate astronomy data
model4_sim <- function(phi = log(c(200, 2.5, 2.5)), J = 360, n.samp = 20){
  x <- y <- bin <- c();
  mu_log_a <- tau_log_a <- rep(0, 4);
  
  for(i in 1:4){
    mu_log_a[i] <- phi[i];
    tau_log_a[i] <- exp(phi[i + 4]);
  }
  
  mu_log_sigma <- phi[9];
  tau_log_sigma <- exp(phi[10]);
  
  for(j in 1:J){
    a_j <- exp(mu_log_a + rnorm(4, 0, 1) * tau_log_a);
    sigma_j <- exp(mu_log_sigma + rnorm(1, 0, 1) * tau_log_sigma);
    
    x_j <- runif(n.samp, 0, 50);
    y_j <- logit_curve(x_j, a_j) + rnorm(n.samp, 0, 1) * sigma_j;
    
    x <- c(x, x_j);
    y <- c(y, y_j);
    bin <- c(bin, rep(j, n.samp));
  }
  
  return(data.frame(x = x, y = y, bin = bin));
}

logit_curve <- function(x = 2, a) {
  result <- a[1] + inv_logit((x - a[2])/a[3]) * a[4];
  return(result)
}

inv_logit <- function(x){ return(1/(1 + exp(-x))) }
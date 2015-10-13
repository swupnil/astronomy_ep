#function to simulate astronomy data
astro_sim <- function(phi = log(c(200, 2.5, 2.5, 250, 5, .5, 50, 7)), J = 360, n.samp = 20){
  x <- y <- bin <- c();
  mu_log_a <- tau_log_a <- rep(0, 8);
  exp_phi <- exp(phi);
  
  for(i in 1:8){
    mu_log_a[i] <- phi[i];
    tau_log_a[i] <- exp(phi[i + 8]);
  }
  
  sigma <- exp_phi[17];
  
  for(j in 1:J){
    a_j <- exp(mu_log_a + rnorm(8, 0, 1) * tau_log_a);
      
    x_j <- runif(n.samp, 0, 50);
    y_j <- exp_curve(x_j, a_j) + rnorm(n.samp, 0, 1) * sigma;
    
    x <- c(x, x_j);
    y <- c(y, y_j);
    bin <- c(bin, rep(j, n.samp));
  }
  
  return(data.frame(x = x, y = y, bin = bin));
}

exp_curve <- function(x = 2, a) {
  result <- a[4] * a[5] * (1 - exp(-x/(a[5])) - a[6] * exp(-.5*((x - a[7])/a[8])^2));
  result <- result * inv_logit((x - a[2])/a[3]) + a[1];
  return(result)
}

inv_logit <- function(x){ return(1/(1 + exp(-x))) }
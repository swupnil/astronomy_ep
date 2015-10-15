#function to simulate astronomy data
model4_sim <- function(pri_mu = log(c(200, 2.5, 2.5)), pri_sig = c(0.5,0.5,0.5),  J = 360, n.samp = 20){
  x <- y <- bin <- c();
  
  phi <- pri_mu + rnorm(10, 0, 1) * diag(pri_sig);
  mu_log_a <- phi[1:5];
  tau_log_a <- exp(phi[6:10]);
  
  print(phi)
  
  for(j in 1:J){
    a_j <- exp(mu_log_a + rnorm(5, 0, 1) * tau_log_a);
    
    x_j <- runif(n.samp, 0, 50);
    y_j <- logit_curve(x_j, a_j) + rnorm(n.samp, 0, 1) * a_j[5];
    
    x <- c(x, x_j);
    y <- c(y, y_j);
    bin <- c(bin, rep(j, n.samp));
  }
  
  return(data.frame(x = x, y = y, bin = bin, phi = phi));
}

logit_curve <- function(x = 2, a) {
  result <- a[1] + inv_logit((x - a[2])/a[3]) * a[4];
  return(result)
}

inv_logit <- function(x){ return(1/(1 + exp(-x))) }
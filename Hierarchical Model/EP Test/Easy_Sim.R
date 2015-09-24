# Model 0
model0_sim <- function(phi = 200, n.samp = 1000){
  y <- rnorm(n.samp, phi, 1);
  data.frame(y = y);
}

# Model 1
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

# Model 2
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
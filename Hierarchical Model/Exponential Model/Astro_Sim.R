#function to simulate astronomy data
astro.sim <- function(phi = log(c(200, 2.5, 2.5, 250, 5, .5, 50, 7)), J = 360, n.samp = 20){
  x <- y <- bin <- c();
  a <- a_tau <- rep(0, 8);
  for(i in 1:8){
    a[i] <- exp(phi)[i];
  }
  for(i in 1:8){
    a_tau[i] <- exp(phi)[i+8];
  }
  sigma <- exp(phi)[17];
  
  for(j in 1:J){
    x <- c(x, runif(n.samp, 0, 25));
    a_j <- a + rnorm(8,0,1)*a_tau;
    eta <- rnorm(n.samp,0,1)*sigma;
    y <- c(y, exp.curve2(x[((j-1)*n.samp+1):(j*n.samp)], a_j) + eta);
    bin <- c(bin, rep(j,n.samp));
  }
  return(data.frame(x=x, y=y, bin=bin));
}
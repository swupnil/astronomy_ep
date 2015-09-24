#function to simulate astronomy data
astro.sim <- function(phi = c(5.41,4.23,5.52,5.58), J = 360, n.samp = 20){
  x <- y <- bin <- c();
  mu.alpha <- exp(phi)[1];
  tau.alpha <- exp(phi)[2];
  beta <- exp(phi)[3];
  sigma <- exp(phi)[4];
  for(j in 1:J){
    x <- c(x, runif(n.samp,0,4));
    alpha <- mu.alpha + rnorm(1,0,1)*tau.alpha;
    eta <- rnorm(n.samp,0,1)*sigma;
    y <- c(y, alpha+beta*x[((j-1)*n.samp+1):(j*n.samp)]+eta);
    bin <- c(bin, rep(j,n.samp));
  }
  return(data.frame(x=x, y=y, bin=bin));
}
poly.fitter <- function(degrees = 36, inc = 1, P = 1, x = i100[i100*FUV>0], y = FUV[i100*FUV>0], z = long[i100*FUV>0], bin.factor = 1, zoomed = TRUE){
  x <- x[z < degrees + inc & z >= degrees]
  y <- y[z < degrees + inc & z >= degrees]
  X <- x;
  for(p in 2:P){
    X <- cbind(X, x^p);
  }
  
  model <- lm(y~X);
  xax <- c(0:(max(x)*100))/100;
  fitted <- rep(model$coeff[1], length(xax));
  for(p in 1:P){
    fitted <- fitted + model$coeff[p+1]*xax^p;
  }
  
  degrees <- degrees / bin.factor
  if(zoomed){
    plot(x, y, xlab = "i100", ylab = "FUV", xlim = c(0, 45), main = bquote(~ P == .(P) ~ "," ~ theta == .(degrees) *degree ~ "-" ~ .(degrees+inc) *degree), cex = 1, bty = "l", col = "gray")
  } else{
    plot(x, y, xlab = "i100", ylab = "FUV", xlim = c(0, 45), ylim = c(0,10000), main = bquote(~ P == .(P) ~ "," ~ theta == .(degrees) *degree ~ "-" ~ .(degrees+inc) *degree), cex = 1, bty = "l", col = "gray")
  }
  lines(xax, fitted, col = 1, lwd = log(length(x))/2, lty = 1)
}

poly.plot <- function(degrees = 36, inc = .1, x = i100[i100*FUV>0], y = FUV[i100*FUV>0], z = long[i100*FUV>0], bin.factor = 1){
  par(mfrow = c(2,3))
  for (p in 1:6)
    poly.fitter(degrees, inc, p, x, y, z, bin.factor)
  par(mfrow = c(1,1))
}

sim.fitter <- function(degrees = 36, inc = 1, P = 5, x = sims$x[sims$x*sims$y>0 & sims$x<200], y = sims$y[sims$x*sims$y>0 & sims$x < 200], z = sims$bin[sims$x*sims$y>0 & sims$x < 200], bin.factor = 1){
  x <- x[z < degrees + inc & z >= degrees]
  y <- y[z < degrees + inc & z >= degrees]
  X <- x;
  for(p in 2:P){
    X <- cbind(X, x^p);
  }
  
  model <- lm(y~X);
  xax <- c(0:(max(x)*100))/100;
  fitted <- rep(model$coeff[1], length(xax));
  for(p in 1:P){
    fitted <- fitted + model$coeff[p+1]*xax^p;
  }
  
  degrees <- degrees / bin.factor
  plot(x, y, xlab = "i100", ylab = "FUV", xlim = c(0, 200), main = bquote(~ P == .(P) ~ "," ~ theta == .(degrees) *degree ~ "-" ~ .(degrees+inc) *degree), cex = 1, bty = "l", col = "gray")
  lines(xax, fitted, col = 1, lwd = log(length(x))/2, lty = 1)
}
exp.fitter <- function(degrees = 150, inc = 1, bin.factor = 1, is.plot = TRUE, fit = exp.fit, x = i100[i100*FUV>0], y = FUV[i100*FUV>0], z = long[i100*FUV>0], mu = rep(1,9), sig = rep(1,9)){
  x <- x[z < degrees + inc & z >= degrees];
  y <- y[z < degrees + inc & z >= degrees];
  
  #Fit the Model
  exp.data <- list(N = length(x), x = x, y = y, priorMu = mu, priorSig = sig);
  init.data <- list();
  for(i in 1:4) { init.data[[i]] <- list(log_a = mu, log_sigma = 5.3);}
  model.fit <- sampling(fit, data = exp.data, iter = 100, chains = 4, init = init.data);
  
  if(is.plot){
    # Plot Raw Data
    degrees <- degrees / bin.factor;
    plot(x, y, xlab = "i100", ylab = "FUV", main = bquote(theta == .(degrees) *degree ~ "-" ~ .(degrees+inc) *degree), cex = 1, bty = "l", col = "gray");
    
    # Find and Plot Fitted Values of Model
    xax <- c(0:(max(x)*100))/100;
    A <- extract(model.fit)$a
    
    for(i in 1:nrow(A)){
      if(fit@model_name == "Exp Fit"){
        fitted <- exp.curve(xax,A[i,]);
      } else {
        fitted <- exp.curve2(xax,A[i,]);
      }
      lines(xax, fitted, lwd = .5, lty = 1, col = alpha("red",.1));
    }
  }
  
  return (model.fit)
}

inv.logit <- function(x){ return(1/(1+exp(-x))) }

exp.curve = function(x = 2, a) {
  result <- a[1] + a[2]*a[3]*(1-exp(-x/(a[3]))-a[4]*exp(-.5*((x-a[5])/a[6])^2)) + a[7]*a[9]*x*log(1+exp((x-a[8])/a[9]))
  return(result)
}

exp.curve2 = function(x = 2, a) {
  result <- a[4]*a[5]*(1-exp(-x/(a[5]))-a[6]*exp(-.5*((x-a[7])/a[8])^2))
  result <- result * inv.logit((x-a[2])/a[3]) + a[1]
  return(result)
}

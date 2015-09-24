#function to run EP on the astronomy data
EP = function(fit = NULL, data = NULL, J = 360, K = 6, prior.Mu = log(c(200, 250, 5, .5, 50, 7, 1, 50, .5)), S = 10, 
                    SMOOTH = 0.9, randomSites = FALSE, parallel = FALSE, mc_iter = 100){

  # Extract data
  x <- data$x;
  y <- data$y;
  bin <- data$bin;
  P <- length(prior.Mu);
  
  # Group the Sites randomly
  bin.samp <- sample(unique(bin), replace = FALSE)
  
  # Initialize global and local natural parameters
  Sigma.k.inv.Mu <- Sigma.k.inv <- Sigma.inv <- Sigma.inv.Mu <- Post.Sigma <- list();
  Post.Mu <- matrix(0, ncol = P, nrow = S + 1);
  for(k in 1:K) Sigma.k.inv.Mu[[k]] <- rep(0,P);
  for(k in 1:K) Sigma.k.inv[[k]] <- diag(P)*0;
  Sigma.inv[[1]] <- Sigma.inv[[2]] <- diag(P)*K/10;
  Sigma.inv.Mu[[1]] <- Sigma.inv.Mu[[2]] <- Sigma.inv[[1]]%*%prior.Mu;
  Post.Sigma[[1]] <- solve(Sigma.inv[[1]]);
  Post.Mu[1,] <- as.vector(Post.Sigma[[1]]%*%Sigma.inv.Mu[[1]]);
  init.data <- list()

  # timer
  init.time <- proc.time();
  times <- matrix(0,nrow=S,ncol=3)

  # The actual algorithm...
  for(s in 1:S){
    #update natural parameters from previous iteration for parallel/serial
    if (parallel) {
      Sigma.inv[[s+1]] <- Sigma.inv[[1]];
      Sigma.inv.Mu[[s+1]] <- Sigma.inv.Mu[[1]];
    } else if (s > 1) {
      Sigma.inv[[s]] <- Sigma.inv[[s-1]];
      Sigma.inv.Mu[[s]] <- Sigma.inv.Mu[[s-1]];
    }
    
    for(k in 1:K){
      cat(" Current Iteration Status: [", s, " out of ", S, "] \n");
      cat(" Current Partition Status: [", k, " out of ", K, "] \n");
      
      # 1. Update the Cavity Distribution...
      Sigma.inv.cav <- Sigma.inv[[s]] - Sigma.k.inv[[k]];
      Sigma.inv.Mu.cav <- Sigma.inv.Mu[[s]] - Sigma.k.inv.Mu[[k]];
      Sigma.cav <- solve(Sigma.inv.cav);
      Sigma.cav <- (Sigma.cav+t(Sigma.cav))/2; #preserve symmetry
      Mu.cav <- as.vector(Sigma.cav%*%Sigma.inv.Mu.cav);
    
      # 2. Find tilted distribution in Stan...
      # Extract current partition of the data...
      if (randomSites) { subset <- which(bin %in% bin.samp[((k-1)*J/K+1) : (k*J/K)]); }
      else { subset <- which(bin <= k*J/K & bin > (k-1)*J/K); }
      y.k <- y[subset];
      x.k <- x[subset];
      bin.k <- (ceiling(bin[subset] - 1))%%(J/K) + 1;
      tilt.data <- list(N=length(y.k), M=P, B=J/K, x=x.k, y=y.k, bin=bin.k, Mu_Cav=Mu.cav, Sig_Cav=Sigma.cav);
    
      # Fit tilted distribution in Stan....
      #for(i in 1:4) { init.data[[i]] <- list(eta = matrix(0, nrow=P/2, ncol=J/K), phi = Mu.cav);}
      for(i in 1:4) { init.data[[i]] <- list(theta = rep(0, J/K), phi = Mu.cav);}
      times[s,] <- system.time(tilt.fit <- sampling(fit, data = tilt.data, iter = mc_iter, chains = 4, init = init.data))[1:3];
    
      # Extract mean and covariance matrix....
      Mu.tilt <- colMeans(extract(tilt.fit)$phi);
      Sigma.tilt <- matrix(0, nrow = P, ncol = P);
      for(i in 1:P)
        for(j in 1:P)
          Sigma.tilt[i,j] <- cov(extract(tilt.fit)$phi[,i], extract(tilt.fit)$phi[,j]);
      # Bias correction on covariance matrix....
      n <- length(extract(tilt.fit)$phi[,1]);
      Sigma.inv.tilt <- solve(Sigma.tilt)*(n-P-2)/(n-1);
      Sigma.inv.Mu.tilt <- Sigma.inv.tilt%*%Mu.tilt;
    
      # 3. Update the Site Distribution
      Sigma.k.inv[[k]] <- Sigma.inv.tilt - Sigma.inv.cav; 
      Sigma.k.inv.Mu[[k]] <- Sigma.inv.Mu.tilt - Sigma.inv.Mu.cav;
      
      # 4. Update g(phi) in Parallel/Serial
      if (parallel) {
        Sigma.inv[[s+1]] <- Sigma.inv.tilt + Sigma.inv[[s+1]];
        Sigma.inv.Mu[[s+1]] <- Sigma.inv.Mu.tilt + Sigma.inv.Mu[[s+1]];
      } else {
        Sigma.inv[[s]] <- SMOOTH*Sigma.inv.tilt + (1-SMOOTH)*Sigma.inv[[s]];
        Sigma.inv.Mu[[s]] <- SMOOTH*Sigma.inv.Mu.tilt + (1-SMOOTH)*Sigma.inv.Mu[[s]];
      }
    }
    
    # Convert Natural Parameters Back Into Usual Framework
    if (parallel) {
      Post.Sigma[[s+1]] = solve(Sigma.inv[[s+1]]);
      Post.Mu[s+1,] = as.vector(Post.Sigma[[s+1]]%*%Sigma.inv.Mu[[s+1]]);
    } else if (s > 1) {
      Post.Sigma[[s+1]] = solve(Sigma.inv[[s]]);
      Post.Mu[s+1,] = as.vector(Post.Sigma[[s+1]]%*%Sigma.inv.Mu[[s]]);
    }
  }

  #save the total time elapsed
  final.time <- proc.time();
  
  #plot phi v. s, and return all values
  EP.plot(exp(Post.Mu), J, K, S);
  return(list(Post.Sigma = Post.Sigma, Post.Mu = Post.Mu,
              time = final.time - init.time, times = times));
}

EP.plot = function(data=NULL, J=360, K=6, S=10, params=expression(mu[alpha],tau[alpha],beta,sigma)){
  P <- ncol(data);
  plot(c(0:S), data[,1], xlab="Iteration", ylab=expression(phi), 
       ylim=c(min(data),max(data)), type='l', main=paste("EP: J =",J,", K =",K));
  for(p in 2:P)
    lines(c(0:S), data[,p], col=p);
  legend("bottomleft", col=c(1:P), legend=params, lty=rep(1,P));
}

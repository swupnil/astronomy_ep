Model0_EP = function(fit = NULL, data = NULL, prior.Mu = 0, prior.Sigma.inv = 0.01, S = 10, 
              K = 5, SMOOTH = 0.9, parallel = TRUE, mc_iter = 100){
  
  # Extract data
  y <- data$y;
  sites <- apply(as.matrix(1:length(y)), 1, function(x) sample(1:K,1));
  prior.Sigma.inv.Mu <- prior.Sigma.inv * prior.Mu;
  
  # Initialize global and local natural parameters
  Sigma.k.inv.Mu <- Sigma.k.inv <- Sigma.inv <- Sigma.inv.Mu <- Post.Sigma <- Post.Mu <- list();
  
  for(k in 1:K) Sigma.k.inv.Mu[[k]] <- 0;
  for(k in 1:K) Sigma.k.inv[[k]] <- 0;
  Sigma.inv[[1]] <- prior.Sigma.inv;
  Sigma.inv.Mu[[1]] <- prior.Sigma.inv.Mu;
  
  Post.Sigma[[1]] <- 1/(Sigma.inv[[1]]);
  Post.Mu[[1]] <- as.vector(Post.Sigma[[1]] * Sigma.inv.Mu[[1]]);
  init.data <- list();
  
  # timer
  init.time <- proc.time();
  times <- matrix(0,nrow=S,ncol=3)
  
  # The actual algorithm...
  for(s in 1:S){
    if (parallel) {
      Sigma.inv[[s+1]] <- prior.Sigma.inv;
      Sigma.inv.Mu[[s+1]] <- prior.Sigma.inv.Mu;
    } else if (s > 1) {
      Sigma.inv[[s]] <- Sigma.inv[[s-1]];
      Sigma.inv.Mu[[s]] <- Sigma.inv.Mu[[s-1]];
    }
    
    for(k in 1:K){
      cat(" Current Iteration Status: [", s, " out of ", S, "] \n");
      cat(" Current Partition Status: [", k, " out of ", K, "] \n");
      
      # 1. Update the Cavity Distribution...
      Sigma.inv.cav <- Sigma.inv[[s]] - Sigma.k.inv[[k]];
      Sigma.inv.Mu.cav <- Sigma.inv.Mu[[s]] - Sigma.k.inv.Mu[[k]];;
      Sigma.cav <- 1/(Sigma.inv.cav);
      Mu.cav <- as.vector(Sigma.cav * Sigma.inv.Mu.cav);
      
      # 2. Find tilted distribution in Stan...
      # Extract current partition of the data...
      y.k <- y[which(sites == k)];
      tilt.data <- list(N = length(y.k), y = y.k, Mu_Cav = Mu.cav, Sig_Cav = Sigma.cav);
      
      # Fit tilted distribution in Stan....
      for(i in 1:4) { init.data[[i]] <- list(phi = Mu.cav);}
      times[s,] <- system.time(tilt.fit <- sampling(fit, data = tilt.data, iter = mc_iter, chains = 4, init = init.data))[1:3];
      
      # Extract mean and covariance matrix....
      Mu.tilt <- mean(extract(tilt.fit)$phi);
      Sigma.tilt <- var(extract(tilt.fit)$phi);
      
      # Bias correction on covariance matrix....
      n <- length(extract(tilt.fit)$phi);
      Sigma.inv.tilt <- 1/(Sigma.tilt) * (n-1-2)/(n-1);
      Sigma.inv.Mu.tilt <- Sigma.inv.tilt * Mu.tilt;
      
      # 3. Update the Site Distribution
      Sigma.k.inv[[k]] <- Sigma.inv.tilt - Sigma.inv.cav; 
      Sigma.k.inv.Mu[[k]] <- Sigma.inv.Mu.tilt - Sigma.inv.Mu.cav;
      
      # 4. Update g(phi) in Parallel/Serial
      if (parallel) {
        Sigma.inv[[s+1]] <- Sigma.inv[[s+1]] + Sigma.k.inv[[k]];
        Sigma.inv.Mu[[s+1]] <- Sigma.inv.Mu[[s+1]] + Sigma.k.inv.Mu[[k]];
      } else {
        Sigma.inv[[s]] <- SMOOTH * Sigma.k.inv[[k]] + (1-SMOOTH) * Sigma.inv.cav;
        Sigma.inv.Mu[[s]] <- SMOOTH * Sigma.k.inv.Mu[[k]] + (1-SMOOTH) * Sigma.inv.Mu.cav;
      }
    }
    
    # Convert Natural Parameters Back Into Usual Framework
    if (parallel) {
      Post.Sigma[[s+1]] = 1/(Sigma.inv[[s+1]]);
      Post.Mu[[s+1]] = Post.Sigma[[s+1]] * Sigma.inv.Mu[[s+1]];
    } else if (s > 1) {
      Post.Sigma[[s+1]] = 1/(Sigma.inv[[s]]);
      Post.Mu[[s+1]] = Post.Sigma[[s+1]] * Sigma.inv.Mu[[s]];
    }
  }
  
  #save the total time elapsed
  final.time <- proc.time();
  
  #plot phi v. s, and return all values
  return(list(Post.Sigma = Post.Sigma, Post.Mu = Post.Mu,
              time = final.time - init.time, times = times));
}
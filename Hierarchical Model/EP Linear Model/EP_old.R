#function to run EP on the astronomy data
astro.EP.old = function(fit = NULL, J = 360, K = 6, n.samp = 20,
                    prior.Mu = rep(5.5,4), S = 10, P = 4, SMOOTH = 0.1){

  # Initialize global and local natural parameters
  Sigma.k.inv.Mu <- Sigma.k.inv <- Sigma.inv <- Sigma.inv.Mu <- Post.Sigma <- list();
  Post.Mu <- matrix(0, ncol = P, nrow = S + 1);
  for(k in 1:K) Sigma.k.inv.Mu[[k]] <- rep(0,P);
  for(k in 1:K) Sigma.k.inv[[k]] <- diag(P)*0;
  Sigma.inv[[1]] <- diag(P)*K/4;
  Sigma.inv.Mu[[1]] <- Sigma.inv[[1]]%*%prior.Mu;
  Post.Sigma[[1]] <- solve(Sigma.inv[[1]]);
  Post.Mu[1,] <- as.vector(Post.Sigma[[1]]%*%Sigma.inv.Mu[[1]]);

  #timer
  init.time <- proc.time();

  # The actual algorithm...
  for(s in 1:S){
    #update natural parameters from previous iteration
    if(s>1){
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
      subset <- which(long <= k*360/K & long > (k-1)*360/K);
      y <- FUV[subset];
      x <- i100[subset];
      bin <- (ceiling(long[subset]/(360/J)) - 1)%%(J/K) + 1;
      #stratified sampling to sample equally from each bin...
      subset <- as.vector(apply(as.matrix(unique(bin)),1,
                        function(x) sample(which(bin==x), min(n.samp,length(which(bin==x))))));
      bin <- bin[subset]; x = x[subset]; y = y[subset];
      tilt.data <- list(N=length(y), M=P, B=J/K, x=x, y=y, bin=bin, Mu_Cav=Mu.cav, Sig_Cav=Sigma.cav);
    
      # Fit tilted distribution in Stan....
      tilt.fit <- sampling(fit, data = tilt.data, iter = 100, chains = 4);
    
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
    
      # 4. Update g(phi)
      Sigma.inv[[s]] <- (1-SMOOTH)*Sigma.inv.tilt + SMOOTH*Sigma.inv[[s]];
      Sigma.inv.Mu[[s]] <- (1-SMOOTH)*Sigma.inv.Mu.tilt + SMOOTH*Sigma.inv.Mu[[s]];
    }
    
    # Convert Natural Parameters Back Into Usual Framework
    Post.Sigma[[s+1]] = solve(Sigma.inv[[s]]);
    Post.Mu[s+1,] = as.vector(Post.Sigma[[s+1]]%*%Sigma.inv.Mu[[s]]);
  }

  #save the total time elapsed
  final.time <- proc.time();
  
  #plot phi v. s, and return all values
  EP.plot(exp(Post.Mu), J, K, S);
  return(list(Post.Sigma = Post.Sigma, Post.Mu = Post.Mu,
              time = final.time - init.time));
}

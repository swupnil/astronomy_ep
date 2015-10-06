#MODEL
model2_code = '
  data{
    int<lower=1> N;
    int<lower=1> B;
    real y[N];
    int<lower=1> bin[N];
    vector[2] Mu_Cav;
    matrix[2,2] Sig_Cav;
  }
  transformed data{
    matrix[2,2] L;
    L <- cholesky_decompose(Sig_Cav);
  }
  parameters {
    vector[2] phi;
    vector[B] theta;
  }
  transformed parameters {
    real<lower=0> sigma;
    sigma <- exp(phi[2]);
  }
  model{
    real theta_hat[N];
    phi ~ multi_normal_cholesky(Mu_Cav, L);
    
    for(b in 1:B){
      theta[b] ~ normal(phi[1], sigma);
    }

    for(n in 1:N){
      theta_hat[n] <- theta[bin[n]];
    }

    y ~ normal(theta_hat, 1);
  }
'
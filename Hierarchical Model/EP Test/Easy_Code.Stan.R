#MODEL
easy_code = '
  data{
    int<lower=1> N;
    int<lower=1> M;
    int<lower=1> B;
    real y[N];
    int<lower=1> bin[N];
    vector[M] Mu_Cav;
    matrix[M,M] Sig_Cav;
  }
  transformed data{
    matrix[M,M] L;
    L <- cholesky_decompose(Sig_Cav);
  }
  parameters{
    matrix[2,B] eta;
    vector[M] phi;
  }
  transformed parameters{
    vector[B] mu;
    real<lower=0> sigma[B];

    for(b in 1:B){
      mu[b] <- phi[1] + eta[1,b] * phi[2];
      sigma[b] <- exp(phi[3] + eta[2,b] * phi[4]);
    }
  }
  model{
    phi ~ multi_normal_cholesky(Mu_Cav,L);
    for(b in 1:B){
      eta[1,b] ~ normal(0,1);
      eta[2,b] ~ normal(0,1);
    }
    for(n in 1:N){
      y[n] ~ normal(mu[bin[n]], sigma[bin[n]]);
    }
  }
'
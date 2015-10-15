#MODEL
model4_code = '
  data{
    int<lower=1> N;
    int<lower=1> B;
    real<lower=0> x[N];
    real y[N];
    int<lower=1> bin[N];
    vector[10] Mu_Cav;
    matrix[10,10] Sig_Cav;
  }

  transformed data{
    matrix[10,10] L;
    L <- cholesky_decompose(Sig_Cav);
  }

  parameters{
    vector[10] phi;
    matrix[5,B] eta;
  }

  transformed parameters{
    matrix[4,B] a;
    vector[B] log_sd;

    for(b in 1:B){
      for(i in 1:4){
        a[i,b] <- phi[i] + exp(phi[i + 5]) * eta[i,b];
      }
      log_sd[b] <- phi[5] + exp(phi[10]) * eta[5,b];
    }
  }

  model{
    vector[N] y_hat;
    vector[N] sigma_hat;

    phi ~ multi_normal_cholesky(Mu_Cav, L);

    for(i in 1:5){
      for(b in 1:B){
        eta[i,b] ~ normal(0,1);
      }
    }

    for(n in 1:N){
      y_hat[n] <- a[1,bin[n]] + inv_logit( (x[n] - a[2,bin[n]]) / a[3,bin[n]] ) * a[4,bin[n]]; 
      sigma_hat[n] <- exp(log_sd[bin[n]]);
    }

    y ~ normal(y_hat, sigma_hat);
  }
'
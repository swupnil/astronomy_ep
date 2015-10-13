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
    matrix[5,B] eta;
    vector[10] phi;
  }

  transformed parameters{
    matrix[4,B] log_a;
    matrix[4,B] a;
    real log_sigma[B];
    real<lower=0, upper=positive_infinity()> sigma[B];

    for(b in 1:B){
      for(i in 1:4){
        log_a[i,b] <- phi[i] + exp(phi[i + 4]) * eta[i,b];
        a[i,b] <- exp(log_a[i,b]);
      }

      log_sigma[b] <- phi[9] + exp(phi[10]) * eta[5,b];
      sigma[b] <- exp(log_sigma[b]);
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
      y_hat[n] <- a[1,bin[n]] + inv_logit((x[n]-a[2,bin[n]])/a[3,bin[n]]); 
      sigma_hat[n] <- sigma[bin[n]];
    }

    y ~ normal(y_hat, sigma_hat);
  }
'
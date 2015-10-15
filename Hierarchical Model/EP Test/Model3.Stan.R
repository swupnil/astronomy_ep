#MODEL
model3_code = '
  data {
    int<lower=1> N;
    int<lower=1> B;
    real y[N];
    real x[N];
    int<lower=1> bin[N];
    vector[6] Mu_Cav;
    matrix[6,6] Sig_Cav;
  }

  transformed data {
    matrix[6,6] L;
    L <- cholesky_decompose(Sig_Cav);
  }

  parameters {
    vector[6] phi;
    matrix[3,B] eta;
  }

  transformed parameters {
    vector[B] beta0;
    vector[B] beta1;
    vector[B] log_sigma;

    for(b in 1:B) {
      beta0[b] <- phi[1] + exp(phi[4]) * eta[1,b];
      beta1[b] <- phi[2] + exp(phi[5]) * eta[2,b];
      log_sigma[b] <- phi[3] + exp(phi[6]) * eta[3,b];
    }
  }

  model {
    real y_hat[N];
    real sigma_hat[N];

    phi ~ multi_normal_cholesky(Mu_Cav, L);

    for(b in 1:B) {
      for(p in 1:3) {
        eta[p,b] ~ normal(0,1);
      }
    }

    for(n in 1:N){
      y_hat[n] <- beta0[bin[n]] + beta1[bin[n]] * x[n];
      sigma_hat[n] <- exp(log_sigma[bin[n]]);
    }

    y ~ normal(y_hat, sigma_hat);
  }
'
#MODEL
model3_code = '
  data {
    int<lower=1> N;
    int<lower=1> B;
    real y[N];
    real x[N];
    int<lower=1> bin[N];
    vector[5] Mu_Cav;
    matrix[5,5] Sig_Cav;
  }

  transformed data {
    matrix[5,5] L;
    L <- cholesky_decompose(Sig_Cav);
  }

  parameters {
    vector[5] phi;
    vector[B] beta0;
    vector[B] beta1;
  }

  transformed parameters {
    vector[2] theta;
    vector<lower=0>[2] tau;
    real<lower=0> sigma;

    theta[1] <- phi[1];
    theta[2] <- phi[2];

    tau[1] <- exp(phi[3]);
    tau[2] <- exp(phi[4]);

    sigma <- exp(phi[5]);
  }

  model {
    real y_hat[N];

    phi ~ multi_normal_cholesky(Mu_Cav, L);

    for(b in 1:B){
      beta0[b] ~ normal(theta[1], tau[1]);
      beta1[b] ~ normal(theta[2], tau[2]);
    }

    for(n in 1:N){
      y_hat[n] <- beta0[bin[n]] + beta1[bin[n]] * x[n];
    }

    y ~ normal(y_hat, sigma);
  }
'
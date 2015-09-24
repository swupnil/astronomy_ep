#MODEL
model1_code = '
  data{
    int<lower=1> N;
    int<lower=1> B;
    real y[N];
    int<lower=1> bin[N];
    real Mu_Cav;
    real Sig_Cav;
  }
  parameters {
    real phi;
    vector[B] theta;
  }
  model{
    real theta_hat[N];
    phi ~ normal(Mu_Cav, sqrt(Sig_Cav));
    theta ~ normal(phi, 10);
    for(n in 1:N){
      theta_hat[n] <- theta[bin[n]];
    }
    y ~ normal(theta_hat, 1);
  }
'
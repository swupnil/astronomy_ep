#MODEL
model0_code = '
  data{
    int<lower=1> N;
    real y[N];
    real Mu_Cav;
    real Sig_Cav;
  }
  parameters{
    real phi;
  }
  model{
    phi ~ normal(Mu_Cav, sqrt(Sig_Cav));
    y ~ normal(phi, 1);
  }
'
#MODEL
tilt_code = '
  data{
    int<lower=1> N;
    int<lower=1> M;
    int<lower=1> B;
    real<lower=0> x[N];
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
    vector[B] eta;
    vector[M] phi;
  }
  transformed parameters{
    vector[B] alpha;
    real<lower=0, upper=positive_infinity()> mu;
    real<lower=0, upper=positive_infinity()> tau;
    real<lower=0, upper=positive_infinity()> sigma;
    real<lower=0, upper=positive_infinity()> beta;
    mu <- exp(phi[1]);
    tau <- exp(phi[2]);
    beta <- exp(phi[3]);
    sigma <- exp(phi[4]);
    alpha <- mu + tau*eta;
  }
  model{
    vector[N] y_hat;
    phi ~ multi_normal_cholesky(Mu_Cav,L);
    eta ~ normal(0,1);
    for(n in 1:N){
      y_hat[n] <- alpha[bin[n]] + x[n]*beta;
    }
    y ~ normal(y_hat, sigma);
  }
'
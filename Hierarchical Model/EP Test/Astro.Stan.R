#MODEL
astro_code = '
  data{
    int<lower=1> N;
    int<lower=1> B;
    real<lower=0> x[N];
    real y[N];
    int<lower=1> bin[N];
    vector[17] Mu_Cav;
    matrix[17,17] Sig_Cav;
  }

  transformed data{
    matrix[17,17] L;
    L <- cholesky_decompose(Sig_Cav);
  }

  parameters{
    matrix[8,B] eta;
    vector[17] phi;
  }

  transformed parameters{
    matrix[8,B] log_a;
    matrix[8,B] a;
    real<lower=0, upper=positive_infinity()> sigma;
    
    sigma <- exp(phi[17]);

    for(b in 1:B){
      for(i in 1:8){
        log_a[i,b] <- phi[i] + exp(phi[i + 8]) * eta[i,b];
        a[i,b] <- exp(log_a[i,b]);
      }
    }
  }

  model{
    vector[N] y_hat;

    phi ~ multi_normal_cholesky(Mu_Cav, L);

    for(i in 1:8){
      for(b in 1:B){
        eta[i,b] ~ normal(0,1);
      }
    }

    for(n in 1:N){
      y_hat[n] <- a[1,bin[n]] + inv_logit((x[n]-a[2,bin[n]])/a[3,bin[n]])*(a[4,bin[n]]*a[5,bin[n]]*(1-exp(-x[n]/a[5,bin[n]])-a[6,bin[n]]*exp(-.5*pow((x[n]-a[7,bin[n]])/a[8,bin[n]],2))));   
    }

    y ~ normal(y_hat, sigma);
  }
'
#exp log normal model
exp_code_normal2 = "
    data {
        int<lower=0> N;
        vector[N] y;
        vector[N] x;
        vector[8] priorMu;
        vector[8] priorSig;
    }
    parameters {
        vector[8] log_a;
        real log_sigma;
    }
    transformed parameters {
        vector[8] a;
        real<lower=0> sigma;
        a <- exp(log_a);
        sigma <- exp(log_sigma);
    }
    model {
        vector[N] y_hat;
        log_a ~ normal(priorMu, priorSig);
        log_sigma ~ normal(5, 0.25);

        for(n in 1:N){
            y_hat[n] <- a[1] + inv_logit((x[n]-a[2])/a[3])*(a[4]*a[5]*(1-exp(-x[n]/a[5])-a[6]*exp(-.5*pow((x[n]-a[7])/a[8],2))));
        }
        y ~ normal(y_hat, sigma);
    }
"
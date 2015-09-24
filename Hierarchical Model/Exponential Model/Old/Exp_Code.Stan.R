#exp log normal model
exp_code_normal = "
    data {
        int<lower=0> N;
        vector[N] y;
        vector[N] x;
        vector[9] priorMu;
        vector[9] priorSig;
    }
    parameters {
        vector[9] log_a;
        real log_sigma;
    }
    transformed parameters {
        vector[9] a;
        real<lower=0> sigma;
        a <- exp(log_a);
        sigma <- exp(log_sigma);
    }
    model {
        vector[N] y_hat;
        log_a ~ normal(priorMu, priorSig);
        log_sigma ~ normal(5, 0.25);

        for(n in 1:N){
            y_hat[n] <- a[1]+a[2]*a[3]*(1-exp(-x[n]/a[3])-a[4]*exp(-.5*pow((x[n]-a[5])/a[6],2)))+a[7]*a[9]*log(1+exp((x[n]-a[8])/a[9]));
        }
        y ~ normal(y_hat, sigma);
    }
"
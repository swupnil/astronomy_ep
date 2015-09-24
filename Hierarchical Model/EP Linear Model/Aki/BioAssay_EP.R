x = c(-.86, -.3, -.05, -.73);
m = c(5, 5, 5, 5);
y = c(0, 1, 3, 5);

x = cbind(rep(1,4),x);
inv_logit = function(a) { 1/(1+exp(-a))};

Sig.Mu.i = Sig.i.inv = list();
for(i in 1:5){
  Sig.Mu.i[[i]] = rep(0,2);
  Sig.i.inv[[i]] = diag(2);
} 

Sig.inv = Reduce("+",Sig.i.inv);
Sig.Mu = rep(0,2);

for(n in 1:10){
  for(i in 1:4){
    #parameters of cavity distribution
    cav.SM = Sig.Mu - Sig.Mu.i[[i]];
    cav.S = Sig.inv - Sig.i.inv[[i]];
    
    #project to one dimensional subspace
    M_i = t(x[i,])%*%(solve(cav.S))%*%cav.SM;
    V_i = t(x[i,])%*%(solve(cav.S))%*%x[i,];
    
    #unnormalized tilted distribution of eta
    q = function(eta){ dnorm(eta,M_i,V_i)*dbinom(y[i],m[i],inv_logit(eta))}
    E0 = integrate(function(x) q(x),M_i-10*sqrt(V_i),M_i+10*sqrt(V_i))$value;
    E1 = integrate(function(x) x*q(x),M_i-10*sqrt(V_i),M_i+10*sqrt(V_i))$value;
    E2 = integrate(function(x) x^2*q(x),M_i-10*sqrt(V_i),M_i+10*sqrt(V_i))$value;
    M = E1/E0;
    V = E2/E0-M^2;
    
    #update moments
    MV.i = M/V-M_i/V_i;
    V.i.inv = 1/V - 1/V_i;
    
    #transform back to full space
    Sig.Mu.i[[i]] = x[i,]*MV.i;
    Sig.i.inv[[i]] = x[i,]%*%t(x[i,])*as.numeric(V.i.inv);
    
    #update g(theta)
    Sig.Mu = cav.SM + Sig.Mu.i[[i]];
    Sig.inv = cav.S + Sig.i.inv[[i]];
  }
}
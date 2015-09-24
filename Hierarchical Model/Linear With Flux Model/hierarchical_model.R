require("FITSio")
require("rstan")
#setwd("~/Dropbox/Documents/School/Columbia/Research/Astronomy")
#data = readFITS("allskyfuv_all.fits")
#data2 = readFITS("stars_fuv_glon.fits")
#load hierarichical_data.RData

#STAR LIGHT INTENSITY VARIABLES
#gl = data2$col[[1]]
#fluxN6 = data2$col[[5]]
#fluxS6 = data2$col[[6]]
#fluxN5 = data2$col[[7]]
#fluxS5 = data2$col[[8]]

fluxN6avg = fluxS6avg = fluxN5avg = fluxS5avg = rep(NA,length(b1))
for(i in 1:length(b1)){
  fluxN6avg[i] = mean(fluxN6[10*i-9:10*i])
  fluxS6avg[i] = mean(fluxS6[10*i-9:10*i])
  fluxN5avg[i] = mean(fluxN5[10*i-9:10*i])
  fluxS5avg[i] = mean(fluxS5[10*i-9:10*i])
} 

#ALL SKY VARIABLES
FUV = data$col[[1]]
i100 = data$col[[2]]
lat = data$col[[3]]
long = data$col[[4]]

x = matrix(NA,nrow=100,ncol=360)
y = x

for(i in 1:360){
  subset = which(long<=i & long>(i-1) & i100>0 & i100<4 & FUV>0)
  index = sample(1:length(subset),100,replace=F)
  x[,i] = i100[subset][index]
  y[,i] = FUV[subset][index]
  if(i%%10==0)
    cat("Longitude",i,"\n")
}
y[which(is.na(y))] = max(FUV[-which(is.na(FUV))])

#MODEL

model_code = '
  data{
    int<lower=1> I;
    int<lower=1> J;
    real<lower=0> f[I];
    real<lower=0> x[J,I];
    real<lower=0> y[J,I];
    int<lower=0,upper=1> squareRoot;
  }
  parameters{
    real<lower=0> alpha[I];
    real<lower=0> beta[I];
    real W[4];
    real<lower=0> sigma[I];
    real<lower=0> mu;
    real<lower=0> tau[3];
  }
  model{
    for(i in 1:I){
      alpha[i] ~ normal(W[1]+W[2]*f[i],tau[1]);
      beta[i] ~ normal(W[3]+W[4]*f[i],tau[2]);
      sigma[i] ~ lognormal(mu,tau[3]);
      for(j in 1:J){
        if(squareRoot==0)
          y[j,i] ~ normal(alpha[i]+beta[i]*x[j,i],sigma[i]);
        else
          y[j,i] ~ normal(sqrt(pow(alpha[i],2)+pow(beta[i]*x[j,i],2)),sigma[i]);
      }
    }
  }
'

#PRINT THE FITS
fit2 = stan(model_code = model_code, data=list(J=100,I=360,x=x,y=y,f=fluxS5, squareRoot=1),iter=1000,chains=4)
fit1 = stan(model_code = model_code, data=list(J=100,I=360,x=x/2,y=y,f=fluxS5/1000, squareRoot=0),iter=1000,chains=4)
fit = fit2

#assess the weights
par(mfrow=c(2,2))
plot(1:2000,extract(fit)$W[,1],xlab="Iteration",ylab=expression(A[1]))
plot(1:2000,extract(fit)$W[,3],xlab="Iteration",ylab=expression(B[1]))
plot(1:2000,extract(fit)$W[,2],xlab="Iteration",ylab=expression(A[2]))
plot(1:2000,extract(fit)$W[,4],xlab="Iteration",ylab=expression(B[2]))

b0 = apply(extract(fit)$W,2,mean)
b1 = apply(extract(fit)$b1,2,mean)
s = apply(extract(fit)$sigma,2,mean)
mu = apply(extract(fit)$mu,2,mean)
tau = apply(extract(fit)$tau,2,mean)

par(mfrow=c(3,1))

#Plot of fitted values
par(mfrow=c(3,3))
nsim = length(extract(fit)$mu)
xgrid = (0:40)/10
for(i in 1:36){
  plot(x[,i],y[,i],xlab="",ylab="",main=i)
  for(j in sample(1:nsim,10,rep=T))
    lines(xgrid,extract(fit)$W[j,3]+extract(fit)$W[j,4]*fluxS5avg[i]+extract(fit)$W[j,1]*xgrid+extract(fit)$W[j,2]*fluxS5avg[i]*xgrid,col=2)
}
plot(x[,1],y[,1])
lines(xgrid,mean(extract(fit)$W[,3])+mean(extract(fit)$W[j,4])*fluxS5[1]+mean(extract(fit)$W[j,1])*xgrid+mean(extract(fit)$W[j,2])*fluxS5[1]*xgrid,col=2)

#Plot of Intercept v. Longitude
plot(1:36*10,b0,xlab="Longitude",ylab=expression(paste(beta[0])),type='l')
lines(1:36*10,rep(mu[1],36),col="red",lty=1)
text(5,mu[1]-.01,expression(mu[0]),col="red",cex=1.1)
lines(1:36*10,rep(mu[1]+tau[1],36),col="red",lty=2)
text(5,mu[1]+tau[1]+.1,expression(mu[0]+tau[0]),col="red",cex=1.1)
lines(1:36*10,rep(mu[1]-tau[1],36),col="red",lty=2)
text(5,mu[1]-tau[1]-.1,expression(mu[0]-tau[0]),col="red",cex=1.1)

#Plot of Slope v. Longitude
plot(1:36*10,b1,xlab="Longitude",ylab=expression(paste(beta[1])),type='l')
lines(1:36*10,rep(mu[2],36),col="red",lty=1)
text(5,mu[2]-20,expression(mu[1]),col="red",cex=1.1)
lines(1:36*10,rep(mu[2]+tau[2],36),col="red",lty=2)
text(5,mu[2]+tau[2]+20,expression(mu[1]+tau[1]),col="red",cex=1.1)
lines(1:36*10,rep(mu[2]-tau[2],36),col="red",lty=2)
text(5,mu[2]-tau[2]-20,expression(mu[1]-tau[1]),col="red",cex=1.1)

#Plot of Log Sigma v. Longitude
plot(1:36*10,log(s),xlab="Longitude",ylab=expression(paste(log(sigma))),type='l')
lines(1:36*10,rep(mu[3],36),col="red",lty=1)
text(5,mu[3]-.1,expression(mu[sigma]),col="red",cex=1.1)
lines(1:36*10,rep(mu[3]+tau[3],36),col="red",lty=2)
text(5,mu[3]+tau[3]+.1,expression(mu[sigma]+tau[sigma]),col="red",cex=1.1)
lines(1:36*10,rep(mu[3]-tau[3],36),col="red",lty=2)
text(5,mu[3]-tau[3]-.1,expression(mu[sigma]-tau[sigma]),col="red",cex=1.1)

#Plot of Slope/Intercept v. Star Light Intensity
par(mfrow=c(2,2))
plot(fluxN6avg,b1,ylab=expression(beta[1]),main=paste("r =",round(cor(fluxN6avg,b1),3)))
plot(fluxS6avg,b1,ylab=expression(beta[1]),main=paste("r =",round(cor(fluxS6avg,b1),3)))
plot(fluxN5avg,b1,ylab=expression(beta[1]),main=paste("r =",round(cor(fluxN5avg,b1),3)))
plot(fluxS5avg,b1,ylab=expression(beta[1]),main=paste("r =",round(cor(fluxS5avg,b1),3)))
plot(fluxN6avg,b0,ylab=expression(beta[0]),main=paste("r =",round(cor(fluxN6avg,b0),3)))
plot(fluxS6avg,b0,ylab=expression(beta[0]),main=paste("r =",round(cor(fluxS6avg,b0),3)))
plot(fluxN5avg,b0,ylab=expression(beta[0]),main=paste("r =",round(cor(fluxN5avg,b0),3)))
plot(fluxS5avg,b0,ylab=expression(beta[0]),main=paste("r =",round(cor(fluxS5avg,b0),3)))

#Slope/Intercept Time Series v. Light Intesnity Time Series
par(mfrow=c(3,2))
plot(1:36*10,b1,xlab="Longitude",ylab=expression(beta[1]),type='l')
plot(1:36*10,b0,xlab="Longitude",ylab=expression(beta[0]),type='l')
plot(1:36*10,fluxN6avg,xlab="Longitude",type='l',col="red")
plot(1:36*10,fluxS6avg,xlab="Longitude",type='l',col="red")
plot(1:36*10,fluxN5avg,xlab="Longitude",type='l',col="red")
plot(1:36*10,fluxS5avg,xlab="Longitude",type='l',col="red")



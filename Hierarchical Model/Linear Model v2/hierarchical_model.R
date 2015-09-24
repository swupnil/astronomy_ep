require("FITSio")
require("rstan")
setwd("~/Dropbox/Documents/School/Columbia/Research/Astronomy")
data = readFITS("allskyfuv_all.fits")
data2 = readFITS("stars_fuv_glon.fits")

#STAR LIGHT INTENSITY VARIABLES
gl = data2$col[[1]]
fluxN6 = data2$col[[5]]
fluxS6 = data2$col[[6]]
fluxN5 = data2$col[[7]]
fluxS5 = data2$col[[8]]

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

x = matrix(NA,nrow=100,ncol=36)
y = x

for(i in 1:36){
  subset = which(long<=i*10 & long>(i-1)*10 & i100>0 & FUV>0)
  index = sample(1:length(subset),100,replace=F)
  x[,i] = i100[subset][index]
  y[,i] = FUV[subset][index]
  if(i%%10==0)
    cat("Longitude",i,"\n")
}
y[which(is.na(y))] = max(FUV[-which(is.na(FUV))])

#LINEAR MODEL

linear_code = '
data{
int<lower=1> I;
int<lower=1> J;
real x[J,I];
real y[J,I];
}
parameters{
real<lower=0> b0[I];
real<lower=0> b1[I];
real<lower=0> sigma[I];
real<lower=0> mu[3];
real<lower=0> tau[3];
}
model{
for(i in 1:I){
b0[i] ~ normal(mu[1],tau[1]);
b1[i] ~ normal(mu[2],tau[2]);
log(sigma[i]) ~ normal(mu[3],tau[3]);
for(j in 1:J){
y[j,i] ~ normal(b0[i]+b1[i]*x[j,i],sigma[i]);
}
}
}
'

#SQUARE ROOT MODEL WITH X^2

sqrt_code = '
data{
int<lower=1> I;
int<lower=1> J;
real x[J,I];
real y[J,I];
}
parameters{
real<lower=0> b0[I];
real<lower=0> b1[I];
real<lower=0> sigma[I];
real<lower=0> mu[3];
real<lower=0> tau[3];
}
model{
for(i in 1:I){
b0[i] ~ normal(mu[1],tau[1]);
b1[i] ~ normal(mu[2],tau[2]);
log(sigma[i]) ~ normal(mu[3],tau[3]);
for(j in 1:J){
y[j,i] ~ normal(sqrt(b0[i]*b0[i]+b1[i]*b1[i]*x[j,i]*x[j,i]),sigma[i]);
}
}
}
'

#PRINT THE FITS
fit2 = stan(model_code = sqrt_code,data=list(J=100,I=36,x=x,y=y),iter=1000,chains=4)
fit1 = stan(model_code = linear_code, data=list(J=100,I=36,x=x,y=y),iter=1000,chains=4)
fit = fit1

b0 = apply(extract(fit)$b0,2,mean)
b1 = apply(extract(fit)$b1,2,mean)
s = apply(extract(fit)$sigma,2,mean)
mu = apply(extract(fit)$mu,2,mean)
tau = apply(extract(fit)$tau,2,mean)

par(mfrow=c(3,1))

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

#########################
### NO INTERCEPT MODEL ##
#########################

#Plot of Slope v. Longitude
plot(1:36*10,b1,xlab="Longitude",ylab=expression(paste(beta[1])),type='l')
lines(1:36*10,rep(mu[1],36),col="red",lty=1)
text(10,mu[1]-20,expression(mu[1]),col="red",cex=1.1)
lines(1:36*10,rep(mu[1]+tau[1],36),col="red",lty=2)
text(10,mu[1]+tau[1]+20,expression(mu[1]+tau[1]),col="red",cex=1.1)
lines(1:36*10,rep(mu[1]-tau[1],36),col="red",lty=2)
text(10,mu[1]-tau[1]-20,expression(mu[1]-tau[1]),col="red",cex=1.1)

#Plot of Log Sigma v. Longitude
plot(1:36*10,log(s),xlab="Longitude",ylab=expression(paste(log(sigma))),type='l')
lines(1:36*10,rep(mu[2],36),col="red",lty=1)
text(10,mu[2]-.1,expression(mu[sigma]),col="red",cex=1.1)
lines(1:36*10,rep(mu[2]+tau[2],36),col="red",lty=2)
text(10,mu[2]+tau[2]+.1,expression(mu[sigma]+tau[sigma]),col="red",cex=1.1)
lines(1:36*10,rep(mu[2]-tau[2],36),col="red",lty=2)
text(10,mu[2]-tau[2]-.1,expression(mu[sigma]-tau[sigma]),col="red",cex=1.1)

#Plot of Slope/Intercept v. Star Light Intensity
par(mfrow=c(2,2))
plot(fluxN6avg,b1,ylab=expression(beta[1]),main=paste("r =",round(cor(fluxN6avg,b1),3)))
plot(fluxS6avg,b1,ylab=expression(beta[1]),main=paste("r =",round(cor(fluxS6avg,b1),3)))
plot(fluxN5avg,b1,ylab=expression(beta[1]),main=paste("r =",round(cor(fluxN5avg,b1),3)))
plot(fluxS5avg,b1,ylab=expression(beta[1]),main=paste("r =",round(cor(fluxS5avg,b1),3)))

#Slope/Intercept Time Series v. Light Intesnity Time Series
par(mfrow=c(3,2))
plot(1:36*10,b1,xlab="Longitude",ylab=expression(beta[1]),type='l')
plot(1:36*10,fluxN6avg,xlab="Longitude",type='l',col="red")
plot(1:36*10,fluxS6avg,xlab="Longitude",type='l',col="red")
plot(1:36*10,fluxN5avg,xlab="Longitude",type='l',col="red")
plot(1:36*10,fluxS5avg,xlab="Longitude",type='l',col="red")



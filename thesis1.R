rm(list=ls())
x1=c(101.3,102,102.7,103.5,104.2,104.9,105.6,106.9,107,107.7,108.5,109.1,109.9,110.6,111.3)
x2=c(25.3,25.5,25.7,25.9,26.1,26.2,26.4,26.6,26.8,27,27.1,27.3,27.5,27.7,27.8)
y=c(34.481,34.369,34.268,34.16,34.215,34.308,34.402,34.479,34.58,34.682,34.78,34.875,34.963,35.054,35.173)
z=data.frame(y,x1,x2)
z$x1_st<-scale(z$x1)
z$x2_st<-scale(z$x2)
z$y_st<-scale(z$y)
#-----------------------
model=lm(y_st~x1_st+x2_st,data=z)
summary(model)
#-----------------------
library(lmtest)
dwtest(model)
#-----------------------
x=pacf(residuals(model), 14)
x
#-----------------------
n= 15
F = numeric(n)
e = residuals(model)
for(t in 1:(n-1)){
F[t] = (e[t]*e[t+1])
}
rho = sum(F)/sum(e^2);rho
#-----------------------
X = as.matrix(z[,4:5])
y=c(1,x$acf)
m1 <- matrix(NA, 15, 15)
for (i in 1:15) {
  m1[i,i:15]= y[1:(15-i+1)]
}
m1
V=pmax(m1, t(m1), na.rm=TRUE)
H = t(X)%*%solve(V)%*%X
E = eigen(H)
E
E$values
E$vectors
#################################### the historical data
data=read.csv("historical data.csv")
data$x1_st<-scale(data$x1)  
data$x2_st<-scale(data$x2)
data$y_st<-scale(data$y)
##------------------
model2=lm(y_st~x1_st+x2_st,data=data)
summary(model2)
beta=as.numeric(model2$coefficients)[2:3]
beta=as.matrix(beta)
x=pacf (residuals(model2), 59)
#-----------------------
X = as.matrix(data[,4:5])
y=c(1,x$acf)
m1 <- matrix(NA, 60, 60)
for (i in 1:60) {
  m1[i,i:60]= y[1:(60-i+1)]
}
m1
V=pmax(m1, t(m1), na.rm=TRUE)
H = t(X)%*%solve(V)%*%X
#-----------------------
n= 60
F = numeric(n)
e = residuals(model2)
for(t in 1:(n-1)){
  F[t] = (e[t]*e[t+1])
}
rho = sum(F)/sum(e^2);rho
#-----------------------???????????? ??????????
X = as.matrix(data[,4:5])
y=as.matrix(data[,6])
g=c(1,x$acf)
m1 <- matrix(NA, 60, 60)
for (i in 1:60) {
  m1[i,i:60]= g[1:(60-i+1)]
}
m1
V=pmax(m1, t(m1), na.rm=TRUE)
H = t(X)%*%solve(V)%*%X
betag=solve(H)%*%t(X)%*%solve(V)%*%y
sigma_had=(t(y-X%*%betag)%*%solve(V)%*%(y-X%*%betag))/58
sigma_had=as.numeric(sigma_had)
####------variance
R=X[48:49,]
r=as.matrix(y[48:49,])
W=V[1:2,1:2]
k=2*sigma_had/(t(betag)%*%betag)
k=as.numeric(k)
M0=solve(H+t(R)%*%solve(W)%*%R)
Mk=solve(H+t(R)%*%solve(W)%*%R+k*diag(2))
var_betag=sigma_had*solve(H)
var_betam=sigma_had*M0
var_betamk=sigma_had*Mk%*%solve(M0)%*%Mk
##-----------
betam=M0%*%(t(X)%*%solve(V)%*%y+t(R)%*%solve(W)%*%r)
betamk=Mk%*%(t(X)%*%solve(V)%*%y+t(R)%*%solve(W)%*%r)
M_betam=sigma_had*M0
SM_betam=sum(diag(M_betam))
M_betamk=sigma_had*Mk%*%solve(M0)%*%Mk+k^2*Mk%*%beta%*%t(beta)%*%Mk
SM_betamk=sum(diag(M_betamk))
M_betag=var_betag
SM_betag=sum(diag(M_betag))  





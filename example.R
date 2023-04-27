# This R script simulates a panel dataset
# and then uses the conventional fixed effect estimation
# and the half-panel jackknife estimation 
# to estimate the coefficients and the 
# impulse response function


rm(list=ls())

# Step 1: Data Generating Process =======================

# Parameter settings

N = 30
T0 = 60
beta = -0.6
rho = 0

delta = 0.2
eta = 0.2
sigmamu = 1
sigmaep = 1
sigmachi = 1

lagX = 1
lagY = 0

# Data initialization

X=matrix(0,T0,N)
Y=matrix(0,T0,N)
mu=matrix(rnorm(N*T0,mean=0,sd=sigmamu),nrow=T0)
epsilon=matrix(rnorm(N*T0,mean=0,sd=sigmaep),nrow=T0)
Xbar=rep(0,1)
alpha=rep(0,N)
X0=rnorm(N,mean=0,sd=1)
chi=rnorm(n=1,mean=0,sd=sigmachi)

# Generate data
for (i in 1:N) {
  
  # x
  for(t in 1:T0){
    if(t==1){
      X[t,i]=delta+rho* X0[i]+mu[t,i]
    } else {
      X[t,i]=delta+rho* X[(t-1),i]+ mu[t,i]
    }
  }
  
  # Y
  Xbar[i]=mean(X[,i])
  alpha[i]=eta*sqrt(T0)*Xbar[i]+chi #fixed effect

  for (t in 1:T0){
    Y[t,i]= alpha[i]+ beta*X[t,i]+epsilon[t,i]
  }

}

YX <- cbind(c(Y), c(X))
data <- as.data.frame(cbind(rep(1:N,each = T0),rep(1:T0,N),YX))
colnames(data) <- c("id","time","Y","X")



# Step 2: Regression =====================================

source("LP_panel.R")

H = 10
Y.name <- c("Y")
X.name <- c("X")

# FE
fit.FE <- LP_panel(data, Y.name = Y.name, X.name = X.name,
                   method = "FE", lagX = lagX, lagY = lagY, H = H)

IRF.FE <- data.frame(est = fit.FE$IRF,
                     se = fit.FE$se,
                     lower = fit.FE$IRF  - 1.96*fit.FE$se, 
                     upper = fit.FE$IRF  - 1.96*fit.FE$se)

# HJ
fit.HJ <- LP_panel(data, Y.name = Y.name, X.name = X.name,
                   method = "HJ", lagX = lagX, lagY = lagY, H = H)

IRF.HJ <- data.frame(est = fit.HJ$IRF,
                     se = fit.HJ$se,
                     lower = fit.HJ$IRF  - 1.96*fit.HJ$se, 
                     upper = fit.HJ$IRF  - 1.96*fit.HJ$se)

# print results ##########
cat("estimated IRF by FE \n")
print(IRF.FE)

cat("estimated IRF by HJ \n")
print(IRF.HJ)



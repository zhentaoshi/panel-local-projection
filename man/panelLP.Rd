% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/panelLP.R
\name{panelLP}
\alias{panelLP}
\title{Panel Local Projection}
\usage{
panelLP(
  data,
  Y.name,
  X.name,
  c.name = NULL,
  id.name = NULL,
  time.name = NULL,
  method = "SPJ",
  te = F,
  lagX = NULL,
  lagY = 1,
  H = 5
)
}
\arguments{
\item{data}{A data frame, containing the panel data set.}

\item{Y.name}{Character. The dependent variable.}

\item{X.name}{Character. The shock variables.}

\item{c.name}{NULL or Character. The control variables. If NULL, the other columns of the
panel data set (except id.name, time.name, Y.name, X.name) are regarded the control variables.}

\item{id.name}{NULL or Character. The column to identify the cross sectional units. If NULL, the first column of the
panel data set will be the variable denoting the cross section.}

\item{time.name}{NULL or Character. The column to identify the time periods. If NULL, the second column of the panel
data set will be the variable denoting the time section.}

\item{method}{Character. Type of method, either "SPJ" (split-panel jackknife, default) or "FE" (conventional within-group demeaned estimation).}

\item{te}{Boolean. FALSE (default) to exclude the time effect from the model.}

\item{lagX}{NULL or Integer. The number of lagged shock variables included as regressors. If NULL, lagX = lagY.}

\item{lagY}{Integer. The number of lagged dependent variables included as regressors.}

\item{H}{Integer. The maximum horizon for impulse response function estimates.}
}
\value{
A list containing:

\item{IRF}{A matrix, containing the estimated impulse responses. Each row stands for one horizon.}

\item{se}{A matrix, containing the standard errors clustered by individuals. Each row stands for one horizon.}
}
\description{
It offers the main function `panelLP()` to implement the panel local projection that includes two methods: the conventional fixed effect (FE) estimator and the split-panel jackknife (SPJ) estimator that eliminates the asymptotical bias and delivers valid statistical inference.
}
\examples{


# Step 1: Data Generating Process =======================

# Sample size
N = 30
T0 = 60

# Parameter settings
beta = -0.6
rho = 0.5

delta = 0.2
eta = 0.2
sigmamu = 1
sigmaep = 1
sigmachi = 1

lagX = 1
lagY = 0

# Data initialization

set.seed(20230501)

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

H = 10
Y.name <- c("Y")
X.name <- c("X")

# FE
fit.FE <- panelLP(data, Y.name = Y.name, X.name = X.name,
                   method = "FE", lagX = lagX, lagY = lagY, H = H)

IRF.FE <- data.frame(est = fit.FE$IRF,
                     se = fit.FE$se,
                     lower = fit.FE$IRF  - 1.96*fit.FE$se,
                     upper = fit.FE$IRF  + 1.96*fit.FE$se)

# SPJ
fit.SPJ <- panelLP(data, Y.name = Y.name, X.name = X.name,
                   method = "SPJ", lagX = lagX, lagY = lagY, H = H)

IRF.SPJ <- data.frame(est = fit.SPJ$IRF,
                     se = fit.SPJ$se,
                     lower = fit.SPJ$IRF  - 1.96*fit.SPJ$se,
                     upper = fit.SPJ$IRF  + 1.96*fit.SPJ$se)

# print results ##########
cat("estimated IRF by FE \n")
print(IRF.FE)

cat("estimated IRF by SPJ \n")
print(IRF.SPJ)

}
\references{
Ziwei Mei, Liugang Sheng and Zhentao Shi (2023). "Implicit Nickell Bias in Panel Local Projection"
}

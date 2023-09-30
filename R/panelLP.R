#' @name panelLP
#' @title Panel Local Projection
#' @description It offers the main function `panelLP()` to implement the panel local projection that includes two methods: the conventional fixed effect (FE) estimator and the split-panel jackknife (SPJ) estimator that eliminates the asymptotical bias and delivers valid statistical inference.
#' @param data      A data frame, containing the panel data set.
#' @param Y.name    Character. The dependent variable.
#' @param X.name    Character. The shock variables.
#' @param c.name    NULL or Character. The control variables. If NULL, the other columns of the
#'                  panel data set (except id.name, time.name, Y.name, X.name) are regarded the control variables.
#' @param id.name   NULL or Character. The column to identify the cross sectional units. If NULL, the first column of the
#'                  panel data set will be the variable denoting the cross section.
#' @param time.name NULL or Character. The column to identify the time periods. If NULL, the second column of the panel
#'                  data set will be the variable denoting the time section.
#' @param method    Character. Type of method, either "SPJ" (split-panel jackknife, default) or "FE" (conventional within-group demeaned estimation).
#' @param te        Boolean. FALSE (default) to exclude the time effect from the model.
#' @param lagX      NULL or Integer. The number of lagged shock variables included as regressors. If NULL, lagX = lagY.
#' @param lagY      Integer. The number of lagged dependent variables included as regressors.
#' @param H         Integer. The maximum horizon for impulse response function estimates.
#'
#' @return A list containing:
#'
#'\item{IRF}{A matrix, containing the estimated impulse responses. Each row stands for one horizon.}
#'
#'\item{se}{A matrix, containing the standard errors clustered by individuals. Each row stands for one horizon.}
#'
#' @references
#'
#' Ziwei Mei, Liugang Sheng and Zhentao Shi (2023). "Implicit Nickell Bias in Panel Local Projection"
#'
#'
#' @export
#'
#'
#' @examples
#'
#'
#' # Step 1: Data Generating Process =======================
#'
#' # Sample size
#' N = 30
#' T0 = 60
#'
#' # Parameter settings
#' beta = -0.6
#' rho = 0.5
#'
#' delta = 0.2
#' eta = 0.2
#' sigmamu = 1
#' sigmaep = 1
#' sigmachi = 1
#'
#' lagX = 1
#' lagY = 0
#'
#' # Data initialization
#'
#' set.seed(20230501)
#'
#' X=matrix(0,T0,N)
#' Y=matrix(0,T0,N)
#' mu=matrix(rnorm(N*T0,mean=0,sd=sigmamu),nrow=T0)
#' epsilon=matrix(rnorm(N*T0,mean=0,sd=sigmaep),nrow=T0)
#' Xbar=rep(0,1)
#' alpha=rep(0,N)
#' X0=rnorm(N,mean=0,sd=1)
#' chi=rnorm(n=1,mean=0,sd=sigmachi)
#'
#' # Generate data
#' for (i in 1:N) {
#'
#'   # x
#'   for(t in 1:T0){
#'     if(t==1){
#'       X[t,i]=delta+rho* X0[i]+mu[t,i]
#'     } else {
#'       X[t,i]=delta+rho* X[(t-1),i]+ mu[t,i]
#'     }
#'   }
#'
#'   # Y
#'   Xbar[i]=mean(X[,i])
#'   alpha[i]=eta*sqrt(T0)*Xbar[i]+chi #fixed effect
#'
#'   for (t in 1:T0){
#'     Y[t,i]= alpha[i]+ beta*X[t,i]+epsilon[t,i]
#'   }
#'
#' }
#'
#' YX <- cbind(c(Y), c(X))
#' data <- as.data.frame(cbind(rep(1:N,each = T0),rep(1:T0,N),YX))
#' colnames(data) <- c("id","time","Y","X")
#'
#'
#'
#' # Step 2: Regression =====================================
#'
#' H = 10
#' Y.name <- c("Y")
#' X.name <- c("X")
#'
#' # FE
#' fit.FE <- panelLP(data, Y.name = Y.name, X.name = X.name,
#'                    method = "FE", lagX = lagX, lagY = lagY, H = H)
#'
#' IRF.FE <- data.frame(est = fit.FE$IRF,
#'                      se = fit.FE$se,
#'                      lower = fit.FE$IRF  - 1.96*fit.FE$se,
#'                      upper = fit.FE$IRF  + 1.96*fit.FE$se)
#'
#' # SPJ
#' fit.SPJ <- panelLP(data, Y.name = Y.name, X.name = X.name,
#'                    method = "SPJ", lagX = lagX, lagY = lagY, H = H)
#'
#' IRF.SPJ <- data.frame(est = fit.SPJ$IRF,
#'                      se = fit.SPJ$se,
#'                      lower = fit.SPJ$IRF  - 1.96*fit.SPJ$se,
#'                      upper = fit.SPJ$IRF  + 1.96*fit.SPJ$se)
#'
#' # print results ##########
#' cat("estimated IRF by FE \n")
#' print(IRF.FE)
#'
#' cat("estimated IRF by SPJ \n")
#' print(IRF.SPJ)
#'


panelLP = function(data,
                    Y.name,
                    X.name,
                    c.name = NULL,
                    id.name = NULL,
                    time.name = NULL,
                    method = "SPJ",
                    te = F,
                    lagX = NULL,
                    lagY = 1,
                    H = 5) {


  # Check NULL
  if(is.null(id.name)){
    id.name <- names(data)[1]
  }
  if(is.null(time.name)){
    time.name <- names(data)[2]
  }
  if (is.null(c.name)){
    c.name <- names(data)[which(!(colnames(data) %in% c(Y.name,X.name,id.name,time.name)))]
  }
  if(is.null(lagX)){
    lagX = max(lagY,1)
  }

  px = length(X.name)
  pc = length(c.name)
  N = length(unique(data[,id.name]))
  T0 <- length(sort(unique(data[,time.name])))
  lag.max <- max(lagX,lagY)
  h.seq = 0:H

  # Prepare matrices to store results
  IRF=matrix(NA, length(h.seq), px)
  se=matrix(NA, length(h.seq), px)

  # Add data to a balance panel
  for (nn in 1:N) {
    data[,id.name][data[,id.name] == unique(data[,id.name])[nn]] <- nn
  }
  for (tt in 1:T0) {
    data[,time.name][data[,time.name] == sort(unique(data[,time.name]))[tt]] <- tt
  }

  data.new <- NULL
  for (nn in 1:N){
    data_nn <- data[data[,id.name] == nn,]
    Missing <- !(1:T0 %in% data_nn[,time.name])
    if (sum(Missing) > 0){
      Missing.ind <- (1:T0)[Missing]
      add <- matrix(NA, nrow=length(Missing.ind), ncol = ncol(data))
      add[,which(names(data) == id.name)] <- rep(nn, length(Missing.ind))
      add[,which(names(data) == time.name)] <- Missing.ind
      colnames(add) <- colnames(data_nn)
      data_nn <- rbind(data_nn, add)
      order <- sort(data_nn[,which(names(data) == time.name)], index.return = T)$ix
      data_nn <- data_nn[order,]
    }
    data.new <- rbind(data.new, data_nn)
  }

  # Create data for regression
  data <- data.new

  y.long <- as.matrix(data[,Y.name])
  X.long <- as.matrix(data[,X.name])
  cc.long <- as.matrix(data[,c.name])

  y <- matrix(y.long, ncol = N)
  x <- NULL
  for (ix in 1:px){
    x = cbind(x, matrix(X.long[,ix], ncol = N))
  }
  cc <- NULL
  if (pc > 0){
    for (ic in 1:pc){
      cc = cbind(cc, matrix(cc.long[,ic], ncol = N))
    }
  }


  #===Function to Demean Data===
  demean = function(type_i, type_want, yt_X_CC_ylag, type="t") {

    # type="t": demean t; type="n": demean n

    yt_X_CC_ylag.dm <- yt_X_CC_ylag
    for (i in 1:type_i) {
      want <- 1:type_want + (i-1)*type_want
      want.mean <- rowSums(is.na(yt_X_CC_ylag[want,])) == 0
      for (iv in 1:ncol(yt_X_CC_ylag)) {
        yt_X_CC_ylag.dm[want,iv] <- yt_X_CC_ylag.dm[want,iv] - mean(yt_X_CC_ylag[want[want.mean],iv])
      }
    }

    yt_X_CC_ylag.dm
  }

  #===Function to Transform Data Arrangements===
  transf = function(N, T0h, yt_X_CC_ylag, type="tnl_ntl") {

    # type="tnl_ntl": T*N to N*T; type="ntl_tnl": N*T to T*N

    if (type == "tnl_ntl") { col=N } else if (type=="ntl_tnl") { col=T0h }

    yt_X_CC_ylag.tr = NULL
    for (iv in 1:ncol(yt_X_CC_ylag)) {
      yt_X_CC_ylag.tr = cbind(yt_X_CC_ylag.tr, matrix(t(matrix(yt_X_CC_ylag[,iv], ncol=col)), nrow=N*T0h))
    }

    yt_X_CC_ylag.tr
  }

  #===Function to Prepare Data for FE Estimation===
  ols_within_dataprepare = function(N, T0h, y_h, x_h, cc_h, pc, lagY, te) {

    yt <- c(y_h[,1:N])
    X <- matrix(x_h, nrow=length(yt))
    CC = NULL
    if (pc >0) {
      CC <- matrix(cc_h, nrow=length(yt))
    }
    ylag <- NULL
    if (lagY >= 1) {
      ylag <- matrix(y_h[,-(1:N)], nrow=length(yt))
    }

    yt_X_CC_ylag = cbind(yt, X, CC, ylag)

    # function--demean_t
    dt=demean(N, T0h, yt_X_CC_ylag, type="t")

    if (te == F) {

      dataprepare = dt

    } else if (te == T) {

      # function--transf_nt+demean_n+transfer_tn
      dataprepare=transf(N, T0h, demean(T0h, N, transf(N, T0h, dt, type="tnl_ntl"), type="n"), type="ntl_tnl")
    }

    dataprepare
  }

  #===Function to Compute the Variance of the Estimated Coefficients===
  var_hat = function(N, T0h, indep_var, res.vec, dd.mat=indep_var) {

    W_N = matrix(0, ncol(indep_var), ncol(indep_var))
    res_N <- matrix(res.vec, ncol = N)

    for (iN in 1:N){
      dd_iN <- as.matrix(dd.mat[(iN-1)*T0h+(1:T0h),])
      want_iN <- !( rowSums(is.na(dd_iN))>0 | is.na(res_N[,iN]))
      if (sum(want_iN) == 1) {
        temp <- as.matrix(dd_iN[want_iN,])%*%res_N[want_iN,iN]
        W_N = W_N + temp%*%t(temp)
      } else {
        W_N = W_N + t(dd_iN[want_iN,])%*%res_N[want_iN,iN]%*%t(res_N[want_iN,iN])%*%dd_iN[want_iN,]
      }
    }

    want_NT <- !( rowSums(is.na(indep_var))>0 | is.na(res.vec))
    smp <- length(res.vec[want_NT,])
    W = W_N * (N/(N-1)*(smp-1)/(smp-ncol(indep_var)))
    temp <- t(indep_var[want_NT,])%*%indep_var[want_NT,]

    var0_mat <- solve(temp) %*% W %*% solve(temp)
  }

  #===Function to Split Data===
  split_var = function(cut, N, var, fhalf=NULL, shalf=NULL) {

    var_half = matrix(NA, nrow(var), ncol(var))
    for (ixx in 1:ncol(var_half)) {
      if (ixx%%N !=0) {
        if (is.na(cut[1,ixx%%N])==0) {
          if(!is.null(fhalf)) {
            var_half[(1:cut[1,ixx%%N]),ixx] = fhalf[(1:cut[1,ixx%%N]),ixx]
          }
          if(!is.null(shalf)) {
            var_half[-(1:cut[1,ixx%%N]),ixx] = shalf[-(1:cut[1,ixx%%N]),ixx]
          }
        }
      } else if (ixx%%N ==0) {
        if (is.na(cut[1,N])==0) {
          if(!is.null(fhalf)) {
            var_half[(1:cut[1,N]),ixx] = fhalf[(1:cut[1,N]),ixx]
          }
          if(!is.null(shalf)) {
            var_half[-(1:cut[1,N]),ixx] = shalf[-(1:cut[1,N]),ixx]
          }
        }
      }
    }
    var_half
  }


  for(ih in 1:length(h.seq)){

    h = h.seq[ih]

    # To avoid an error when h = 0
    if (lagY == 0){
      hh=h
    }else if (lagY > 0){
      hh=max(h,1)
    }

    T0h=T0-hh-lag.max+1

    y_h <- y[(hh+lag.max):(T0),]
    if (lagY >= 1){
      for (ilag in 0:(lagY-1)){
        if (h == 0){
          y.temp <- y[(hh+lag.max-ilag-1):(T0-ilag-1),]
        }else{
          y.temp <- y[(lag.max-ilag):(T0-hh-ilag),]
        }
        y_h <- cbind(y_h, y.temp)
      }
    }

    x_h <- NULL
    for (ipx in 1:px){
      x_h_ipx <- NULL
      if (h == 0){
        x_h_ipx <- x[(hh+lag.max):(T0),(ipx-1)*N+1:N]
      }else{
        for (ilag in 0:(lagX-1)){
          x.temp <- x[(lag.max-ilag):(T0-hh-ilag),(ipx-1)*N+1:N]
          x_h_ipx <- cbind(x_h_ipx, x.temp)
        }
      }
      x_h <- cbind(x_h, x_h_ipx)
    }

    if (pc > 0){
      cc_h=cc[(hh+lag.max):T0,]
    }



    if(method == "FE") {

      # function--ols_within_dataprepare
      FEdata = ols_within_dataprepare(N, T0h, y_h, x_h, cc_h, pc, lagY, te)

      dep_var = FEdata[,1]
      indep_var = as.matrix(FEdata[,-1])
      fit_h=lm(dep_var ~ indep_var + 0)
      res.vec = dep_var - indep_var %*% fit_h$coefficients

      if (h == 0){
        IRF[ih,] = fit_h$coefficients[1:px]
      }else{
        IRF[ih,] = fit_h$coefficients[(1:px-1)*lagX + 1]
      }

      # function--var_hat
      var.hat=var_hat(N, T0h, indep_var, res.vec, dd.mat=indep_var)

      if (h == 0){
        se[ih,] = sqrt(diag( as.matrix(var.hat) ))[1:px]
      }else{
        se[ih,] = sqrt(diag( as.matrix(var.hat) ))[(1:px-1)*lagX + 1]
      }



    } else if (method == "SPJ") {

      # function--ols_within_dataprepare
      SPJdata = ols_within_dataprepare(N, T0h, y_h, x_h, cc_h, pc, lagY, te)

      dep_var = SPJdata[,1]
      indep_var = as.matrix(SPJdata[,-1])
      fit_h=lm(dep_var ~ indep_var + 0)
      beta_all=fit_h$coefficients

      # Compute cut point
      cut = matrix(NA,1,ncol(matrix(complete.cases(cbind(dep_var,indep_var)), ncol = N)))
      for (iN in 1:ncol(cut)) {
        cut[1,iN] = floor(median(which(matrix(complete.cases(cbind(dep_var,indep_var)), ncol = N)[,iN]==TRUE)))
      }


      ### Regression for sample a
      # function--split_var
      y_h_a = split_var(cut, N, y_h, fhalf=y_h, shalf=NULL)
      x_h_a = split_var(cut, N, x_h, fhalf=x_h, shalf=NULL)
      if (pc >0){
        cc_h_a = split_var(cut, N, cc_h, fhalf=cc_h, shalf=NULL)
      }

      # function--ols_within_dataprepare
      SPJdata_a = ols_within_dataprepare(N, T0h, y_h_a, x_h_a, cc_h_a, pc, lagY, te)

      dep_var_a = SPJdata_a[,1]
      indep_var_a = as.matrix(SPJdata_a[,-1])
      fit_h_a=lm(dep_var_a ~ indep_var_a + 0)
      beta_a = fit_h_a$coefficients


      ### Regression for sample b
      # function--split_var
      y_h_b = split_var(cut, N, y_h, fhalf=NULL, shalf=y_h)
      x_h_b = split_var(cut, N, x_h, fhalf=NULL, shalf=x_h)
      if (pc >0){
        cc_h_b = split_var(cut, N, cc_h, fhalf=NULL, shalf=cc_h)
      }

      # function--ols_within_dataprepare
      SPJdata_b = ols_within_dataprepare(N, T0h, y_h_b, x_h_b, cc_h_b, pc, lagY, te)

      dep_var_b = SPJdata_b[,1]
      indep_var_b = as.matrix(SPJdata_b[,-1])
      fit_h_b=lm(dep_var_b ~ indep_var_b + 0)
      beta_b = fit_h_b$coefficients


      ### Estimate IRF
      beta.hat = 2*beta_all-0.5*(beta_a+beta_b)

      if (h == 0){
        IRF[ih,] = beta.hat[1:px]
      }else{
        IRF[ih,] = beta.hat[(1:px-1)*lagX + 1]
      }

      ### S.E.
      res.vec <- dep_var - indep_var %*% beta.hat

      XX.mat <- matrix(indep_var, nrow = T0h)
      XX.mat.a <- matrix(indep_var_a, nrow = T0h)
      XX.mat.b <- matrix(indep_var_b, nrow = T0h)
      # function--split_var
      XX.mat.sub = split_var(cut, N, XX.mat, fhalf=XX.mat.a, shalf=XX.mat.b)
      dd.mat <- 2*XX.mat - XX.mat.sub
      dd.mat <- matrix(dd.mat, ncol = ncol(indep_var))

      # function--var_hat
      var.hat=var_hat(N, T0h, indep_var, res.vec, dd.mat=dd.mat)

      if (h == 0){
        se[ih,] = sqrt( diag(var.hat ))[1:px]
      }else{
        se[ih,] = sqrt( diag(var.hat ))[(1:px-1)*lagX + 1]
      }



    } else {
      cat("pls specify the method. ")
    }
  }

  return(list(IRF=IRF, se=se))
}

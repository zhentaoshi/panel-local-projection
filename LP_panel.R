
# this file contains one function: LP_panel()
# this function is developed as a companion of 

# Ziwei Mei, Liugang Sheng and Zhentao Shi (2023): "Implicit Nickell Bias in Panel Local Projection"

# Qiyu Dai, Ziwei Mei, and Shu Shen contribute to the code.


LP_panel = function(data,Y.name,X.name,
                    c.name = NULL,
                    id.name = NULL,
                    time.name = NULL,
                    H = 10,
                    lagY = 1, 
                    lagX = NULL,
                    method = "HJ") {

# INPUTS
  # data: a data frame for all variables
  # Y.name: the dependent variable in the panel LP regression
  # X.name: the shock variables in the panel LP regression
  # c.name: the control variables in the panel LP regression
  # id.name: the column to identify the cross sectional units
  # time.name: the column to identify the time periods
  # H: the maximum horizon for prediction
  # lagY: how many lagged dependent variables are to be included as regressors
  # lagX: how many lagged shock variables are to be included as regressors
  # method: either "HJ" (half-panel jackknife) or "FE" (conventional within-group estimation)
  
  
# OUTPUTS
  # IRF: the estimated coefficients of the shock variables
  # se: the standard errors of the estimated coefficients allowing for the within-group correlation
  
  h.seq=0:H
  px = length(X.name)
  IRF=matrix(NA,length(h.seq),px)
  se=matrix(NA,length(h.seq),px)
  
  if (is.null(id.name)){
    id.name <- names(data)[1]
  }
  if (is.null(time.name)){
    time.name <- names(data)[2]
  }
  if (is.null(lagX)){
    lagX = max(lagY,1)
  }
  lag.max <- max(lagX,lagY)
  id.temp <- data[,id.name]
  N = length(unique(id.temp))
  for (nn in 1:N){
    data[,id.name][id.temp == unique(id.temp)[nn]] <- nn 
  }
  time.temp <- data[,time.name]
  T0 <- length( sort(unique(time.temp)) )
  for (tt in 1:T0){
    data[,time.name][time.temp == sort(unique(time.temp))[tt]] <- tt
  }
  
  id <- data[,id.name]
  time <- data[,time.name]
  
  data.new <- NULL 
  for (nn in 1:N){
    data_nn <- data[data[,id.name] == nn,]
    Missing <- !(1:T0 %in% data_nn[,time.name])
    if (sum(Missing) > 0){
      Missing.ind <- (1:T0)[Missing]
      add <- matrix(NA,nrow=length(Missing.ind),ncol = ncol(data) )
      add[,which(names(data) == id.name)] <- rep(nn,length(Missing.ind))
      add[,which(names(data) == time.name)] <- Missing.ind
      colnames(add) <- colnames(data_nn)
      data_nn <- rbind(data_nn,add)
      order <- sort(data_nn[,which(names(data) == time.name)],index.return = T)$ix 
      data_nn <- data_nn[order,] 
    }
    
    data.new <- rbind(data.new, data_nn)
  }
  
  data <- data.new 
  
  id <- data[,id.name]
  time <- data[,time.name]
  
  y.long <- as.matrix(data[,Y.name])
  X.long <- as.matrix(data[,X.name]) 
  
  
  if (is.null(c.name)){
    cc.ind <- which(!(colnames(data) %in% c(Y.name,X.name,id.name,time.name)))
    c.name <- names(data)[cc.ind]
  }
  pc = length(c.name)
  cc.long <- as.matrix(data[ , c.name])
  
  y <- matrix(y.long, ncol = N)
  
  x <- NULL
  
  for (ix in 1:px){
    x.ix.long <- X.long[,ix]
    x.ix <- matrix(x.ix.long, ncol = N)
    x = cbind(x,x.ix)
  }
  
  cc <- NULL
  # cc.demean <- NULL 
  if (pc > 0){
    for (ic in 1:pc){
      cc.ic.long <- cc.long[,ic]
      cc.ic <- matrix(cc.ic.long, ncol = N)
      cc = cbind(cc,cc.ic)
    }
  }
  
  ### function: demean_t =================================================== ###
  demean_t = function(N, y_h, yt, X, CC, ylag, pc, lagY) {
    
    yt.dm <- yt
    X.dm <- X
    CC.dm <- CC
    ylag.dm <- ylag 
    
    for (iN in 1:N){
      want <- 1:nrow(y_h) + (iN-1)*nrow(y_h)
      want.mean <- rowSums(is.na(cbind(yt[want],X[want,],CC[want,],ylag[want,]))) == 0
      yt.dm[want] <-  yt.dm[want] - mean(yt[want[want.mean]])
      for (ix in 1:ncol(X)){
        X.dm[want,ix] <-  X.dm[want,ix] - mean(X[want[want.mean],ix])
      }
      
      if (pc >0){
        for (ic in 1:ncol(CC)){
          CC.dm[want,ic] <-  CC.dm[want,ic] -  mean(CC[want[want.mean],ic])
        }
      }
      if (lagY >= 1){
        for (iylag in 1:ncol(ylag)){
          ylag.dm[want,iylag] <-  ylag.dm[want,iylag] - mean(ylag[want[want.mean],iylag]) 
        }
      }
    }
    
    return(list(yt.dm=yt.dm, X.dm=X.dm, CC.dm=CC.dm, ylag.dm=ylag.dm))
  }
  ### ====================================================================== ###
  
  for(ih in 1:length(h.seq)){
    
    h = h.seq[ih]
    
    # When h=0, it is necessary to reduce one period of samples to ensure that lagy does not report an error
    if (lagY == 0){
      hh=h
    }else if (lagY > 0){
      hh=max(h,1)
    }
    
    lag_T=T0-hh
    
    y_h <- y[(hh+lag.max):(T0),]
    
    
    if (lagY >= 1){
      for (ilag in 0:(lagY-1)){
        
        if (h == 0){
          y.temp <- y[(hh+lag.max-ilag-1):(T0-ilag-1),]
        }else{
          y.temp <- y[(lag.max-ilag):(lag_T-ilag),]
        }
        
        
        y_h <- cbind(y_h,y.temp)
        
      }
    }
    
    
    x_h <- NULL 
    
    
    for (ipx in 1:px){
      x_h_ipx <- NULL
      
      if (h == 0){
        x_h_ipx <- x[(hh+lag.max):(T0),(ipx-1)*N+1:N]
      }else{
        for (ilag in 0:(lagX-1)){
          x.temp <- x[(lag.max-ilag):(lag_T-ilag),(ipx-1)*N+1:N]
          x_h_ipx <- cbind(x_h_ipx,x.temp)
        }
      }
      
      x_h <- cbind(x_h,x_h_ipx)
      
    }
    
    if (pc > 0){
      cc_h=cc[(hh+lag.max):T0,]
    }
    
    if(method == "FE") {
      
      yt <- c(y_h[,1:N])
      
      X <- matrix(x_h, nrow = length(yt) )
      
      CC = NULL 
      if (pc >0){
        CC <- matrix(cc_h, nrow = length(yt) )
      }
      
      ylag <- NULL
      if (lagY >= 1){
        ylag <- matrix(y_h[,-(1:N)], nrow = length(yt) )
      }
      
      # function:demean_t
      dt.FE=demean_t(N, y_h, yt, X, CC, ylag, pc, lagY)
      
      dep_var = dt.FE$yt.dm
      indep_var <- cbind(dt.FE$X.dm,dt.FE$CC.dm,dt.FE$ylag.dm)
      
      fit_h=lm(dep_var ~ indep_var + 0)
      res.vec = dep_var - indep_var %*% fit_h$coefficients
      
      if (h == 0){
        IRF[ih,] = fit_h$coefficients[1:px]
      }else{
        IRF[ih,] = fit_h$coefficients[(1:px-1)*lagX + 1]
      }
      
      W_N = matrix(0,ncol(indep_var),ncol(indep_var))
      res_N <- matrix(res.vec,ncol = N)
      for (iN in 1:N){ 
        indep_iN <- as.matrix(indep_var[(iN-1)*nrow(y_h)+(1:nrow(y_h)),])
        want_iN <- !( rowSums(is.na(indep_iN)) > 0 | is.na(res_N[,iN]) )
        
        if (sum(want_iN) == 1){
          temp_N <- as.matrix(indep_iN[want_iN,])%*%res_N[want_iN,iN]
          W_N = W_N + temp_N %*% t(temp_N)
        }else{
          W_N = W_N + t(indep_iN[want_iN,])%*%res_N[want_iN,iN]%*%t(res_N[want_iN,iN])%*%indep_iN[want_iN,]
        }
      }
      
      want_NT <- !( rowSums(is.na(indep_var)) > 0 | is.na(res.vec) )
      smp <- length(res.vec[want_NT,])
      
      W = W_N *(N/(N-1)*(smp-1)/(smp-ncol(indep_var)))
      
      temp <- t(indep_var[want_NT,])%*%indep_var[want_NT,]
      var0_mat <- solve(temp)%*%W%*%solve(temp)
      
      var_mat <- var0_mat
      
      if (h == 0){
        se[ih,] = sqrt(diag( as.matrix(var_mat) ))[1:px]
      }else{
        se[ih,] = sqrt(diag( as.matrix(var_mat) ))[(1:px-1)*lagX + 1]
      }
      
    } else if (method == "HJ") {
      
      ### function: split_var ============================================== ###
      split_var = function(cut, N, var, half=1) {
        #half: "1" first half, "-1" second half
        
        var_half = matrix(NA,nrow(var),ncol(var))
        for (ixx in 1:ncol(var_half)) {
          if (ixx%%N !=0) {
            if (is.na(cut[1,ixx%%N])==0) {
              var_half[half*(1:cut[1,ixx%%N]),ixx]=var[half*(1:cut[1,ixx%%N]),ixx]
            }
          } else if (ixx%%N ==0) {
            if (is.na(cut[1,N])==0) {
              var_half[half*(1:cut[1,N]),ixx]=var[half*(1:cut[1,N]),ixx]
            }
          }
        }
        
        return(list(var_half=var_half))
      }
      ### ================================================================== ###
      
      yt <- c(y_h[,1:N])
      
      X <- matrix(x_h, nrow = length(yt) )
      
      CC = NULL 
      if (pc >0){
        CC <- matrix(cc_h, nrow = length(yt) )
      }
      
      ylag <- NULL
      if (lagY >= 1){
        ylag <- matrix(y_h[,-(1:N)], nrow = length(yt) )
      }
      
      # function:demean_t
      dt.HJ=demean_t(N, y_h, yt, X, CC, ylag, pc, lagY)
      
      dep_var = dt.HJ$yt.dm
      indep_var <- cbind(dt.HJ$X.dm,dt.HJ$CC.dm,dt.HJ$ylag.dm)
      
      fit_h=lm(dep_var ~ indep_var + 0)
      beta_all=fit_h$coefficients 
      
      ######################
      T0h <- nrow(y_h)
      
      smp=cbind(dep_var,indep_var)
      smpna = complete.cases(smp)
      smpna.ti <- matrix(smpna, ncol = N)
      cut = matrix(NA,1,ncol(smpna.ti))
      for (iN in 1:ncol(cut)) {
        cut[1,iN] = floor(median(which(smpna.ti[,iN]==TRUE)))
      }
      
      ################## Sample a 
      
      #function:split_var
      y_a=split_var(cut, N, y_h, half=1)
      y_h_a = y_a$var_half
      yt_a <- c(y_h_a[,1:N])
      
      x_a=split_var(cut, N, x_h, half=1)
      x_h_a = x_a$var_half
      X_a <- matrix(x_h_a, nrow = length(yt_a))
      
      CC_a = NULL
      if (pc >0){
        cc_a=split_var(cut, N, cc_h, half=1)
        cc_h_a = cc_a$var_half
        CC_a <- matrix(cc_h_a, nrow = length(yt_a) )
      }
      
      ylag_a <- NULL
      if (lagY >= 1){
        ylag_a <- matrix(y_h_a[,-(1:N)], nrow = length(yt_a) )
      }
      
      # function:demean_t
      dt.HJ_a=demean_t(N, y_h_a, yt_a, X_a, CC_a, ylag_a, pc, lagY)
      
      dep_var_a = dt.HJ_a$yt.dm
      indep_var_a <- cbind(dt.HJ_a$X.dm,dt.HJ_a$CC.dm,dt.HJ_a$ylag.dm)
      
      fit_h_a=lm(dep_var_a ~ indep_var_a + 0)
      beta_a = fit_h_a$coefficients 
      
      ################## Sample b
      
      # function:split_var
      y_b=split_var(cut, N, y_h, half=-1)
      y_h_b = y_b$var_half
      yt_b <- c(y_h_b[,1:N])
      
      x_b=split_var(cut, N, x_h, half=-1)
      x_h_b = x_b$var_half
      X_b <- matrix(x_h_b, nrow = length(yt_b))
      
      CC_b = NULL
      if (pc >0){
        cc_b=split_var(cut, N, cc_h, half=-1)
        cc_h_b = cc_b$var_half
        CC_b <- matrix(cc_h_b, nrow = length(yt_b) )
      }
      
      ylag_b <- NULL
      if (lagY >= 1){
        ylag_b <- matrix(y_h_b[,-(1:N)], nrow = length(yt_b) )
      }
      
      # function:demean_t
      dt.HJ_b=demean_t(N, y_h_b, yt_b, X_b, CC_b, ylag_b, pc, lagY)
      
      dep_var_b = dt.HJ_b$yt.dm
      indep_var_b <- cbind(dt.HJ_b$X.dm,dt.HJ_b$CC.dm,dt.HJ_b$ylag.dm)
      
      fit_h_b=lm(dep_var_b ~ indep_var_b + 0)
      beta_b = fit_h_b$coefficients 
      
      ######### estimate IRF 
      beta.hat = 2*beta_all-0.5*(beta_a+beta_b)
      
      if (h == 0){
        IRF[ih,] = beta.hat[1:px]
      }else{
        IRF[ih,] = beta.hat[(1:px-1)*lagX + 1]
      }
      ################### s.e. 
      
      res.vec <- dep_var - indep_var %*% beta.hat
      
      XX.mat <- matrix(indep_var, nrow = T0h)
      # apply(x_h,MARGIN=2,FUN=mean)
      XX.mat.a <- matrix(indep_var_a, nrow = T0h)
      XX.mat.b <- matrix(indep_var_b, nrow = T0h)
      
      XX.mat.sub = matrix(NA,nrow(XX.mat),ncol(XX.mat))
      for (ixx in 1:ncol(XX.mat.sub)) {
        if (ixx%%N !=0) {
          if (is.na(cut[1,ixx%%N])==0) {
            XX.mat.sub[(1:cut[1,ixx%%N]),ixx]=XX.mat.a[(1:cut[1,ixx%%N]),ixx]
            XX.mat.sub[-(1:cut[1,ixx%%N]),ixx]=XX.mat.b[-(1:cut[1,ixx%%N]),ixx]
          }
        } else if (ixx%%N ==0) {
          if (is.na(cut[1,N])==0) {
            XX.mat.sub[(1:cut[1,N]),ixx]=XX.mat.a[(1:cut[1,N]),ixx]
            XX.mat.sub[-(1:cut[1,N]),ixx]=XX.mat.b[-(1:cut[1,N]),ixx]
          }
        }
      }
      
      dd.mat <- 2*XX.mat - XX.mat.sub
      dd.mat <- matrix(dd.mat, ncol = ncol(indep_var) )
      
      R_N.hat = matrix(0,ncol(indep_var),ncol(indep_var))
      res_N <- matrix(res.vec, ncol = N)
      for (iN in 1:N){ 
        dd_iN <- as.matrix(dd.mat[(iN-1)*nrow(y_h)+(1:nrow(y_h)),])
        want_iN <- !( rowSums(is.na(dd_iN)) > 0 | is.na(res_N[,iN]) )
        
        if (sum(want_iN) == 1){
          temp <- as.matrix(dd_iN[want_iN,])%*%res_N[want_iN,iN]
          R_N.hat = R_N.hat + temp %*% t(temp)
        }else{
          R_N.hat = R_N.hat+t(dd_iN[want_iN,])%*%res_N[want_iN,iN]%*%t(res_N[want_iN,iN])%*%dd_iN[want_iN,]
        }
      }
      
      want_NT <- !( rowSums(is.na(indep_var)) > 0 | is.na(res.vec) )
      smp <- length(res.vec[want_NT,])
      
      R.hat = R_N.hat *(N/(N-1)*(smp-1)/(smp-ncol(indep_var)))
      
      Q.hat = t(indep_var[want_NT,]) %*% indep_var[want_NT,]
      var0.hat <- solve(Q.hat) %*% R.hat %*% solve(Q.hat)
      
      var.hat <- var0.hat
      
      if (h == 0){
        se[ih,] = sqrt( diag(var.hat ))[1:px]
      }else{
        se[ih,] = sqrt( diag(var.hat ))[(1:px-1)*lagX + 1]
      }
      
    } else{
      cat("pls speficy the method")
    }
  }
  
  results=list(IRF,se)
  names(results)=c("IRF","se")
  return(results)
}
#' @name panelLP
#' @title Panel Local Projection
#' @description It offers the main function `panelLP()` to implement the panel local projection that includes two methods: the
#'              conventional fixed effect (FE) estimator and the split-panel jackknife (SPJ) estimator that eliminates the
#'              asymptotical bias and delivers valid statistical inference.
#' @param data      A data frame, containing the panel data set sorted by the cross sectional units and the time periods.
#' @param Y.name    Character. The dependent variable.
#' @param X.name    Character. The shock variables.
#' @param c.name    NULL (default) or Character. The control variables.
#' @param id.name   NULL (default) or Character. The column to identify the cross sectional units. If NULL, the first column of
#'                  the panel data set will be the variable denoting the cross section.
#' @param time.name NULL (default) or Character. The column to identify the time periods. If NULL, the second column of the panel
#'                  data set will be the variable denoting the time section.
#' @param method    Character. "SPJ" (default) or "FE". Type of method, either "SPJ" (split-panel jackknife, default) or "FE"
#'                  (conventional within-group demeaned estimation).
#' @param te        Boolean. FALSE (default) to exclude the time effect from the model.
#' @param lagX      NULL (default) or Integer. The number of lagged shock variables included as regressors. If NULL or 0, no
#'                  lagged shock variables will be included.
#' @param lagY      NULL (default) or Integer. The number of lagged dependent variables included as regressors. If NULL or 0, no
#'                  lagged dependent variables will be included.
#' @param H         Integer. The horizon for impulse response function estimates is 0:H.
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
#' rho = 0.8
#'
#' delta = 0.2
#' eta = 0.2
#' sigmamu = 1
#' sigmaep = 1
#' sigmachi = 1
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
#' fit.FE <- panelLP(data, Y.name = Y.name, X.name = X.name, method = "FE", H = H)
#'
#' results.FE <- data.frame(IRF = fit.FE$IRF,
#'                          se = fit.FE$se,
#'                          lower = fit.FE$IRF  - 1.96*fit.FE$se,
#'                          upper = fit.FE$IRF  + 1.96*fit.FE$se)
#'
#' # SPJ
#' fit.SPJ <- panelLP(data, Y.name = Y.name, X.name = X.name, method = "SPJ", H = H)
#'
#' results.SPJ <- data.frame(IRF = fit.SPJ$IRF,
#'                           se = fit.SPJ$se,
#'                           lower = fit.SPJ$IRF  - 1.96*fit.SPJ$se,
#'                           upper = fit.SPJ$IRF  + 1.96*fit.SPJ$se)
#'
#' # print results ##########
#' cat("estimated IRF by FE \n")
#' print(results.FE)
#'
#' cat("estimated IRF by SPJ \n")
#' print(results.SPJ)
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
                   lagY = NULL,
                   H ) {


  # Check NULL
  if(is.null(id.name)){
    id.name <- names(data)[1]
  }
  if(is.null(time.name)){
    time.name <- names(data)[2]
  }
  if(is.null(lagX)){
    lagX = 0
  }
  if(is.null(lagY)){
    lagY = 0
  }

  # Obtain basic parameters
  px = length(X.name)
  pc = length(c.name)
  N = length(unique(data[,id.name]))
  T0 <- length(sort(unique(data[,time.name])))
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
  data <- data.new

  # Create data for regression
  y_tn <- as.matrix(data[,Y.name])
  x_tn_p <- as.matrix(data[,X.name])
  c_tn_p <- as.matrix(data[,c.name])

  # Create lagged matrix for Y and X

  #===Function to Lag Each Column of a Matrix===
  lag_matrix <- function(x, n) {
    # Initialize a lagged matrix with NA values
    lag_x <- matrix(NA, nrow = nrow(x), ncol = ncol(x))

    # Apply lag operation for each column
    for (i in 1:ncol(x)) {
      if (n >= nrow(x)) {
        # If lag period is greater than or equal to the number of rows, set entire column to NA
        lag_x[, i] <- rep(NA, nrow(x))
      } else {
        # Otherwise, perform the lag operation
        lag_x[, i] <- c(rep(NA, n), x[1:(nrow(x) - n), i])
      }
    }

    # Return the lagged matrix
    return(lag_x)
  }

  # Reshape y and x into matrix
  ymat_t_n <- matrix(y_tn, ncol = N) # T rows and N columns
  xmat_t_np <- NULL  # T rows and N*px columns
  for (ix in 1:px){
    xmat_t_np = cbind(xmat_t_np, matrix(x_tn_p[,ix], ncol = N))
  }

  # Create lagged matrix for y and x
  ylagmat_t_nl <- NULL
  if (lagY >= 1) {
    for (ilagY in 1:lagY) {
      ylagmat_t_nl = cbind(ylagmat_t_nl, lag_matrix(ymat_t_n, ilagY))
    }
  }
  xlagmat_t_npl <- NULL
  if (lagX >= 1) {
    for (ilagX in 1:lagX) {
      xlagmat_t_npl = cbind(xlagmat_t_npl, lag_matrix(xmat_t_np, ilagX))
    }
  }



  #===Function to Prepare Data for Within-group OLS Estimation===
  ols_within_dataprepare = function(N, T0, yf_x_ylag_xlag_c, te) {

    #===Function to Demean Data===
    demean = function(type_i, type_want, yf_x_ylag_xlag_c) {

      # type_i=t: demean t; type_i=n: demean n

      yf_x_ylag_xlag_c.dm <- yf_x_ylag_xlag_c
      for (i in 1:type_i) {
        want <- 1:type_want + (i-1)*type_want
        want.mean <- rowSums(is.na(yf_x_ylag_xlag_c[want,])) == 0
        for (iv in 1:ncol(yf_x_ylag_xlag_c)) {
          yf_x_ylag_xlag_c.dm[want,iv] <- yf_x_ylag_xlag_c.dm[want,iv] - mean(yf_x_ylag_xlag_c[want[want.mean],iv])
        }
      }

      yf_x_ylag_xlag_c.dm
    }

    # function--demean_t
    data.dn=demean(N, T0, yf_x_ylag_xlag_c)

    if (te == F) {

      dataprepare = data.dn

    } else if (te == T) {

      #===Function to Transform Data Arrangements===
      transf = function(N, T0, yf_x_ylag_xlag_c, type="tnl_ntl") {

        # type="tnl_ntl": T*N to N*T; type="ntl_tnl": N*T to T*N

        if (type == "tnl_ntl") { col=N } else if (type=="ntl_tnl") { col=T0 }

        yf_x_ylag_xlag_c.tr = NULL
        for (iv in 1:ncol(yf_x_ylag_xlag_c)) {
          yf_x_ylag_xlag_c.tr = cbind(yf_x_ylag_xlag_c.tr, matrix(t(matrix(yf_x_ylag_xlag_c[,iv], ncol=col)), nrow=N*T0))
        }

        yf_x_ylag_xlag_c.tr
      }

      # function--transf_nt+demean_n+transfer_tn
      dataprepare=transf(N, T0, demean(T0, N, transf(N, T0, data.dn, type="tnl_ntl")), type="ntl_tnl")
    }

    dataprepare
  }


  #===Function to Compute the Variance of the Estimated Coefficients===
  var_hat = function(N, T0, indep_var, res.vec, dd.mat=indep_var) {

    W_N = matrix(0, ncol(indep_var), ncol(indep_var))
    res_N <- matrix(res.vec, ncol = N)

    for (iN in 1:N){
      dd_iN <- as.matrix(dd.mat[(iN-1)*T0+(1:T0),])
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

    # Create leaded matrix for y

    #===Function to lead each column of a matrix by a specified number of periods===
    lead_matrix <- function(x, n) {
      # Initialize a leaded matrix with NA values
      lead_x <- matrix(NA, nrow = nrow(x), ncol = ncol(x))

      # Apply lead operation for each column
      for (i in 1:ncol(x)) {
        if (n >= nrow(x)) {
          # If lead period is greater than or equal to the number of rows, set entire column to NA
          lead_x[, i] <- rep(NA, nrow(x))
        } else {
          # Otherwise, perform the lead operation
          lead_x[, i] <- c(x[(n + 1):nrow(x), i], rep(NA, n))
        }
      }

      # Return the leaded matrix
      return(lead_x)
    }

    # Create leads for y
    yfmat_t_n = lead_matrix(ymat_t_n, h)

    # Create data for regression
    yf_tn <- matrix(yfmat_t_n, nrow=N*T0)

    ylag_tn_l <- NULL
    if (lagY >= 1) {
      ylag_tn_l <- matrix(ylagmat_t_nl, nrow=N*T0)
    }
    xlag_tn_pl <- NULL
    if (lagX >= 1) {
      xlag_tn_pl <- matrix(xlagmat_t_npl, nrow=N*T0)
    }

    yf_x_ylag_xlag_c = cbind(yf_tn, x_tn_p, ylag_tn_l, xlag_tn_pl, c_tn_p)



    if(method == "FE") {

      # function--ols_within_dataprepare
      FEdata = ols_within_dataprepare(N, T0, yf_x_ylag_xlag_c, te)

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
      var.hat=var_hat(N, T0, indep_var, res.vec, dd.mat=indep_var)

      if (h == 0){
        se[ih,] = sqrt(diag( as.matrix(var.hat) ))[1:px]
      }else{
        se[ih,] = sqrt(diag( as.matrix(var.hat) ))[(1:px-1)*lagX + 1]
      }



    } else if (method == "SPJ") {

      # function--ols_within_dataprepare
      SPJdata = ols_within_dataprepare(N, T0, yf_x_ylag_xlag_c, te)

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

      yf_x_ylag_xlag_c_a = NULL
      for (iv in 1:ncol(yf_x_ylag_xlag_c)) {
        var = matrix(yf_x_ylag_xlag_c[,iv], ncol=N)
        # function--split_var
        var_half = split_var(cut, N, var, fhalf=var, shalf=NULL)
        yf_x_ylag_xlag_c_a = cbind(yf_x_ylag_xlag_c_a, c(var_half))
      }

      # function--ols_within_dataprepare
      SPJdata_a = ols_within_dataprepare(N, T0, yf_x_ylag_xlag_c_a, te)

      dep_var_a = SPJdata_a[,1]
      indep_var_a = as.matrix(SPJdata_a[,-1])
      fit_h_a=lm(dep_var_a ~ indep_var_a + 0)
      beta_a = fit_h_a$coefficients


      ### Regression for sample b
      # function--split_var
      yf_x_ylag_xlag_c_b = NULL
      for (iv in 1:ncol(yf_x_ylag_xlag_c)) {
        var = matrix(yf_x_ylag_xlag_c[,iv], ncol=N)
        # function--split_var
        var_half = split_var(cut, N, var, fhalf=NULL, shalf=var)
        yf_x_ylag_xlag_c_b = cbind(yf_x_ylag_xlag_c_b, c(var_half))
      }


      # function--ols_within_dataprepare
      SPJdata_b = ols_within_dataprepare(N, T0, yf_x_ylag_xlag_c_b, te)

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

      XX.mat <- matrix(indep_var, nrow = T0)
      XX.mat.a <- matrix(indep_var_a, nrow = T0)
      XX.mat.b <- matrix(indep_var_b, nrow = T0)
      # function--split_var
      XX.mat.sub = split_var(cut, N, XX.mat, fhalf=XX.mat.a, shalf=XX.mat.b)
      dd.mat <- 2*XX.mat - XX.mat.sub
      dd.mat <- matrix(dd.mat, ncol = ncol(indep_var))

      # function--var_hat
      var.hat=var_hat(N, T0, indep_var, res.vec, dd.mat=dd.mat)

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

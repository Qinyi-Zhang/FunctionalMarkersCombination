## Functions for compute scores--------
get.phi <- function(Cov,t, FVE.th = 99.99){
  eig <- eigen(Cov)
  positiveInd <- eig[['values']] >= 0
  if (sum(positiveInd) == 0) {
    stop('All eigenvalues are negative. The covariance estimate is incorrect.')
  }
  d <- eig[['values']][positiveInd]
  eigenV <- eig[['vectors']][, positiveInd, drop=FALSE]
  
  maxK <- ncol(eigenV)
  d <- d[1:maxK]
  eigenV <- eigenV[, 1:maxK, drop=FALSE]
  
  # thresholding for corresponding FVE option 
  #(not before to avoid not being able to reach the FVEthreshold when pos eigenvalues > maxk)
  # i.e. default FVE 0.9999 outputs all components remained here.
  FVE <- cumsum(d) / sum(d) * 100  # cumulative FVE for all available eigenvalues from fitted cov
  no_opt <- min(which(FVE >= FVE.th)) # final number of component chosen based on FVE

  # normalization
  muWork = 1:dim(eigenV)[1]
  phi0 <- apply(eigenV, 2, function(x) {
    x <- x / sqrt(trapzRcpp(t, x^2)) 
    if ( 0 <= sum(x*muWork) )
      return(x)
    else
      return(-x)
  })
  return(list(phi=phi0,FVE=FVE,no_opt=no_opt,eigenV=eigenV,d=d))
}

get.xi <- function(X,eigenV,phi0,t){
  xi <- NULL
  for (i in 1:dim(eigenV)[2]){
    xi <- cbind(xi,apply(t(matrix(rep(phi0[,i],nrow(X)),length(t)))*X,1,function(x){return(trapzRcpp(t,x))}))
  }
  return(xi)
}

get.score <- function(X, Y, FVE.th = 99.99, T0 = 0, T1 = 1, X.new = NULL, Y.new = NULL){
  n0 <- dim(X)[1]
  n1 <- dim(Y)[1]
  p <- dim(X)[3]
  no.grid <- dim(X)[2]
  t.grid <- seq(T0,T1,length.out = no.grid)
  
  zeta.hat <- list()
  eta.hat <- list()
  score.var <- list()
  
  zeta.hat.test <- eta.hat.test <- list()
  
  for(j in 1:p){
    SX <- Sparsify(X[ , ,j],t.grid,sparsity = no.grid)
    SY <- Sparsify(Y[ , ,j],t.grid,sparsity = no.grid)
    
    resX <- FPCA(SX$Ly, SX$Lt)
    resY <- FPCA(SY$Ly, SY$Lt)
    
    cov.est <- (n0)/(n0+n1)*(resX$fittedCov) + (n1)/(n0+n1)*(resY$fittedCov)
  
    #get phi
    temp.phi <- get.phi(cov.est,t.grid, FVE.th)
    no_opt <- temp.phi$no_opt
    phi0 <- temp.phi$phi[ ,1:no_opt]
    #FVE <- temp.phi$FVE
    eigenV <- temp.phi$eigenV[ ,1:no_opt]
    score.var[[j]]  <- eigenValue <- temp.phi$d[1:no_opt]
    
    #get score
    zeta.hat[[j]] <- get.xi(X[ , ,j],eigenV,phi0,t.grid)
    eta.hat[[j]] <- get.xi(Y[ , ,j],eigenV,phi0,t.grid)
    
    if(!is.null(X.new)){
      zeta.hat.test[[j]] <- get.xi(X.new[,,j],eigenV,phi0,t.grid)
      eta.hat.test[[j]] <- get.xi(Y.new[,,j],eigenV,phi0,t.grid)
    }    
  }
  return(list(zeta.hat = zeta.hat, eta.hat = eta.hat, score.var = score.var, zeta.hat.test = zeta.hat.test, eta.hat.test = eta.hat.test))
}

get.square.dist <- function(zeta.hat, eta.hat, weight = NULL){
  n0 <- nrow(zeta.hat[[1]])
  n1 <- nrow(eta.hat[[1]])
  p <- length(zeta.hat)
  
  mu.hat <- lapply(zeta.hat, colMeans)
  nu.hat <- lapply(eta.hat, colMeans)
  
  D <- matrix(NA, n0+n1, p)
  
  if(is.null(weight)) {
    for(j in 1:p){
      D[1:n0,j] <- apply(
        ((zeta.hat[[j]] - matrix(rep(mu.hat[[j]], n0),n0,byrow = T))^2 -
           (zeta.hat[[j]] - matrix(rep(nu.hat[[j]], n0),n0,byrow = T))^2), 
        1, sum)
      D[(n0+1):(n1+n0),j] <- apply(
        ((eta.hat[[j]] - matrix(rep(mu.hat[[j]], n1),n1,byrow = T))^2 -
           (eta.hat[[j]] - matrix(rep(nu.hat[[j]], n1),n1,byrow = T))^2), 
        1, sum)
    }
  } else {
    for(j in 1:p){
      D[1:n0,j] <- apply(
        ((zeta.hat[[j]] - matrix(rep(mu.hat[[j]], n0),n0,byrow = T))^2 -
           (zeta.hat[[j]] - matrix(rep(nu.hat[[j]], n0),n0,byrow = T))^2) / matrix(rep(sqrt(weight[[j]]),n0),n0,byrow = T), 
        1, sum)
      D[(n0+1):(n1+n0),j] <- apply(
        ((eta.hat[[j]] - matrix(rep(mu.hat[[j]], n1),n1,byrow = T))^2 -
           (eta.hat[[j]] - matrix(rep(nu.hat[[j]], n1),n1,byrow = T))^2) / matrix(rep(sqrt(weight[[j]]),n1),n1,byrow = T), 
        1, sum)
    }
  }
  return(D)
}

get.square.dist.cv <- function(zeta.hat, eta.hat, weight = NULL, zeta.hat.test, eta.hat.test){
  
  n0 <- nrow(zeta.hat[[1]])
  n1 <- nrow(eta.hat[[1]])
  n0.test <- nrow(zeta.hat.test[[1]])
  n1.test <- nrow(eta.hat.test[[1]])
  p <- length(zeta.hat)
  
  mu.hat <- lapply(zeta.hat, colMeans)
  nu.hat <- lapply(eta.hat, colMeans)
  
  D <- matrix(NA, n0+n1, p)
  D.test <- matrix(NA, n0.test+n1.test, p)
  
  if(is.null(weight)) {
    for(j in 1:p){
      D[1:n0,j] <- apply(
        ((zeta.hat[[j]] - matrix(rep(mu.hat[[j]], n0),n0,byrow = T))^2 -
           (zeta.hat[[j]] - matrix(rep(nu.hat[[j]], n0),n0,byrow = T))^2), 
        1, sum)
      D[(n0+1):(n1+n0),j] <- apply(
        ((eta.hat[[j]] - matrix(rep(mu.hat[[j]], n1),n1,byrow = T))^2 -
           (eta.hat[[j]] - matrix(rep(nu.hat[[j]], n1),n1,byrow = T))^2), 
        1, sum)
    }
  } else {
    for(j in 1:p){
      D[1:n0,j] <- apply(
        ((zeta.hat[[j]] - matrix(rep(mu.hat[[j]], n0),n0,byrow = T))^2 -
           (zeta.hat[[j]] - matrix(rep(nu.hat[[j]], n0),n0,byrow = T))^2) / matrix(rep(sqrt(weight[[j]]),n0),n0,byrow = T), 
        1, sum)
      D[(n0+1):(n1+n0),j] <- apply(
        ((eta.hat[[j]] - matrix(rep(mu.hat[[j]], n1),n1,byrow = T))^2 -
           (eta.hat[[j]] - matrix(rep(nu.hat[[j]], n1),n1,byrow = T))^2) / matrix(rep(sqrt(weight[[j]]),n1),n1,byrow = T), 
        1, sum)
    }
  }
  
  if(is.null(weight)) {
    for(j in 1:p){
      D.test[1:n0.test,j] <- apply(
        ((zeta.hat.test[[j]] - matrix(rep(mu.hat[[j]], n0.test),n0.test,byrow = T))^2 -
           (zeta.hat.test[[j]] - matrix(rep(nu.hat[[j]], n0.test),n0.test,byrow = T))^2), 
        1, sum)
      D.test[(n0.test+1):(n1.test+n0.test),j] <- apply(
        ((eta.hat.test[[j]] - matrix(rep(mu.hat[[j]], n1.test),n1.test,byrow = T))^2 -
           (eta.hat.test[[j]] - matrix(rep(nu.hat[[j]], n1.test),n1.test,byrow = T))^2), 
        1, sum)
    }
  } else {
    for(j in 1:p){
      D.test[1:n0.test,j] <- apply(
        ((zeta.hat.test[[j]] - matrix(rep(mu.hat[[j]], n0.test),n0.test,byrow = T))^2 -
           (zeta.hat.test[[j]] - matrix(rep(nu.hat[[j]], n0.test),n0.test,byrow = T))^2) / matrix(rep(sqrt(weight[[j]]),n0.test),n0.test,byrow = T), 
        1, sum)
      D.test[(n0.test+1):(n1.test+n0.test),j] <- apply(
        ((eta.hat.test[[j]] - matrix(rep(mu.hat[[j]], n1.test),n1.test,byrow = T))^2 -
           (eta.hat.test[[j]] - matrix(rep(nu.hat[[j]], n1.test),n1.test,byrow = T))^2) / matrix(rep(sqrt(weight[[j]]),n1.test),n1.test,byrow = T), 
        1, sum)
    }
  }
  return(list(D=D,D.test=D.test))
}

## Mann-Whitney U stat for AUC: continuous data----------
nonp.auc <- function(u,v){
n1  =  length(u)
n2  =  length(v)
return(sum(sapply(u,function(x) sum(x<v)))/n1/n2)
} 

####### data.1, data.2 must be of two-column

nonpar.combine2.auc <- function(alpha,rate,data.1,data.2){ 
  n1  =  nrow(data.1)
  n2  =  nrow(data.2)
  new.1  =  data.1%*%c(alpha,rate)
  new.2  =  data.2%*%c(alpha,rate)
  nonp.auc(new.1,new.2)
}

nonpar.combine2.coef <- function(new.1,new.2,evalnum = 201) { 
  rate = seq(-1,1,length = evalnum)
  alpha = rev(rate)[-1]
  auc.rate_x  =  sapply(rate, nonpar.combine2.auc, alpha = 1,data.1 = new.1,data.2 = new.2) 
  auc.alpha_x  =  sapply(alpha, nonpar.combine2.auc, rate = 1,data.1 = new.1,data.2 = new.2) 
  auc.0  =  c(auc.rate_x,auc.alpha_x) 
  amax.idx  =  which.max(auc.0) 
  if(amax.idx <= evalnum) return(c(alpha = 1,rate = rate[amax.idx], 
                                   auc.max = auc.0[amax.idx]))
  if(amax.idx > evalnum) return(c(alpha = alpha[amax.idx-evalnum],rate = 1, 
                                  auc.max = auc.0[amax.idx])) 
}

nonp.auc.check <- function(health,middle) { 
  auc.i = numeric(ncol(health))
  for (i in 1:ncol(health)) { 
    new.1 = health[,i]
    new.2 = middle[,i]
    auc.i[i] = nonp.auc(new.1,new.2)
    }
    auc.i
}

### Step-wise method
step.coef <- function(new.1,new.2,design = 'step-down'){
  n1  =  nrow(new.1)
  n2  =  nrow(new.2)
  VARnum  =  ncol(new.1)
  combcoef  =  matrix(0,nrow = VARnum-1,ncol = 2)
  
  if (design == 'step-down') {
  auc.order  =  sort(nonp.auc.check(health = new.1,middle = new.2),index.return = T,decreasing = T)$ix
  } 
  else {
    auc.order  =  sort(nonp.auc.check(health = new.1,middle = new.2),index.return = T,decreasing = F)$ix
    } 
  combmarker.1 = new.1[,auc.order[1]] 
  combmarker.2 = new.2[,auc.order[1]] 
  nal.coef  =  1
  for (i in 2:VARnum){
    combmarker.1  =  cbind(combmarker.1,new.1[,auc.order[i]]) 
    combmarker.2  =  cbind(combmarker.2,new.2[,auc.order[i]]) 
    temp.info  =  nonpar.combine2.coef(combmarker.1,combmarker.2) 
    combcoef[i-1,]  =  temp.info[1:2]
    nal.coef  =  c(nal.coef*combcoef[i-1,1],combcoef[i-1,2]) 
    combmarker.1  =  combmarker.1%*%temp.info[1:2] 
    combmarker.2  =  combmarker.2%*%temp.info[1:2]
}
nal.coef  =  nal.coef[sort(auc.order,index.return = T)$ix] 
check.sign  =  nonp.auc(new.1%*%nal.coef,new.2%*%nal.coef) 
if(check.sign <= 0.5) nal.coef = -nal.coef 
return(list(coef = as.numeric(nal.coef),
            auc.combined = as.numeric(temp.info[3]),
            check = (max(check.sign,1-check.sign) == temp.info[3]) ))
}


## Min-Max method
liu.coef <- function(data.1,data.2) {
  max_min.1  =  cbind(apply(data.1,1,max),apply(data.1,1,min))
  max_min.2  =  cbind(apply(data.2,1,max),apply(data.2,1,min))
  est.coef  =  nonpar.combine2.coef(max_min.1,max_min.2)[1:2]
  check.sign  =  nonp.auc(max_min.1%*%est.coef,max_min.2%*%est.coef)
  if(check.sign >= 0.5) 
    return(list(coef = est.coef, auc.combined = check.sign)) 
  else 
    return(list(coef = -est.coef, auc.combined = 1-check.sign))
}
## Methods by Kang ei al. 2016, SIM END

get.score.step <- function(labels, D, D.test){
  time.start=Sys.time()
  step.comb <- step.coef(D[labels==0, ], D[labels==1, ])
  coef.step <- step.comb$coef
  score.step <- coef.step%*%t(D.test)
  time.end=Sys.time()
  time.step=time.end-time.start
  
  return(list(score.step = score.step, time.step = time.step, coef.step = coef.step))
}

get.score.min.max <- function(labels, D, D.test){
  time.start=Sys.time()
  min.max <- liu.coef(D[labels==0, ], D[labels==1, ])
  coef.min.max <- min.max$coef
  score.min.max <- coef.min.max%*%t(cbind(apply(D.test,1,max), apply(D.test,1,min)))
  time.end=Sys.time()
  time.min.max=time.end-time.start
  
  return(list(score.min.max = score.min.max, time.min.max = time.min.max, coef.min.max = coef.min.max))
}

get.score.lg <- function(ytr, xtr, xt){
  time.start=Sys.time()
  lam.lg <- cv.glmnet(x=xtr, y=ytr, 
                          family = "binomial" , alpha = 0.5)$lambda.min
  res.lg.reg <- glmnet(x=xtr, y=ytr, 
                           family = "binomial" , lambda = lam.lg, alpha = 0.5)
  score.lg = predict(res.lg.reg, newx=xt)
  time.end=Sys.time()
  time.lg = time.end-time.start
  return( list( score.lg = score.lg , time.lg=time.lg, coef.lg=res.lg.reg$beta) )
}


get.score.elr <- function(labels, D, D.test, delta = 1/5){
  data.x = D[labels == 0, ]
  data.y = D[labels == 1, ]
  time.start = Sys.time()
  n <- nrow(data.x)
  m <- nrow(data.y)
  h <- (m+n)^(-delta) # define bandwidth
  combine<-function(data,lambda) data%*%as.vector(lambda)
  # kernel functions where the normal kernel is considered
  # and can be replace by other kernels
  K.h <- function(x) pnorm (x/h,0,1)
  k.h <- function(x) dnorm (x/h,0,1)/h
  # the objective function to be optimized
  find.lam.kernel <- function(lambda){
    x.c <- combine(data.x,lambda)
    y.c <- combine(data.y,lambda)
    -mean(mapply(function(x) K.h(x-y.c),x.c) )
  }
  # the gradient function
  grn <- function(lambda){
    x.c <- combine(data.x,lambda)
    y.c <- combine(data.y,lambda)
    xy.diff.k <- mapply(function(x) k.h(x-y.c),x.c)
    -sapply(1:ncol(data.x),function(i) {
      mean(xy.diff.k*mapply(function(x) (x-data.y[,i]),data.x[,i]) )
    })
  }
  # obtain BLC coefficients
  initial = rep(1, ncol(data.x)) / ncol(data.x)
  opt.kernel <- optim(initial,find.lam.kernel,method="BFGS",gr=grn,
                      control=list(maxit=30000, ndeps=.1))
  A.kernel <- -opt.kernel$value
  a.kernel <- opt.kernel$par
  
  time.end = Sys.time()
  time.elr = time.end - time.start
  score.elr = D.test %*% a.kernel
  return(list(time.elr = time.elr, score.elr = score.elr, coef.elr = a.kernel))
}

get.score.pl <- function(labels, D, D.test){
  data.x = D[labels == 0, ]
  data.y = D[labels == 1, ]
  time.start = Sys.time()
  Sigma.x = var(data.x)
  Sigma.y = var(data.y)
  mu.x = colMeans(data.x)
  mu.y = colMeans(data.y)
  coef.pl = solve(Sigma.x + Sigma.y) %*% (mu.y - mu.x)
  time.end = Sys.time()
  time.pl = time.end - time.start
  score.pl = c( D.test %*% coef.pl )
  
  return(list(time.pl = time.pl, score.pl = score.pl, coef.pl = coef.pl))
}

## Functions for data generation-------
gen.data <- function(case, n0, n1, p, K = 50, T0 = 0, T1 = 1, no.grid = 101, j0 = 1, j1 = 2){
  # case: simulation setting
  # K: no. of bases
  # p: no. of biomarkers
  # n0: sample size of group 0
  # n1: sample size of group 1
  
  #return: 
  #X and Y: n*no.grid*p arrays
  
  phi <-  make.basis(K = K, T0 = T0, T1 = T1, no.grid = no.grid)
  
  zeta.eta.temp <- zeta.eta(case, p, n0, n1, K)
  zeta <- zeta.eta.temp$zeta
  eta <- zeta.eta.temp$eta
  
  V0 <- array(0, dim = c(n0,no.grid,p))
  V1 <- array(0, dim = c(n1,no.grid,p))
  for(i in 1:p){
    V0[ , ,i] <- zeta[,i,] %*% t(phi)
    V1[ , ,i] <- eta[,i,] %*% t(phi)
  }
  
  X <- array(rnorm(n0*no.grid*p), dim = c(n0,no.grid,p)) #define X and Y with measurement errors
  Y <- array(rnorm(n1*no.grid*p), dim = c(n1,no.grid,p))
  
  t.grid <- seq(T0, T1, length.out = no.grid)
  
  mu0 <- 0
  # mu1 <- array(rep(runif(p,min = 0,max = 1),each = n1*no.grid), dim = c(n1,no.grid,p))
  mu1 <- array(1, dim = c(n1,no.grid,p))
  tt <- array(rep(t.grid, each = n1*p),dim = c(n1,p,no.grid))
  tt <- aperm(tt, c(1,3,2))
  mu1 <- mu1 * tt #*sqrt(2)
  
  V0j0 <- array(V0[,,j0], dim = c(n0,no.grid,p))
  V1j0 <- array(V1[,,j0], dim = c(n1,no.grid,p))
  if (p >=2){
    V0j1 <- array(V0[,,j1], dim = c(n0,no.grid,p))
    V1j1 <- array(V1[,,j1], dim = c(n1,no.grid,p))
  }
  else{
    V0j1 <- array(0, dim = c(n0,no.grid,p))
    V1j1 <- array(0, dim = c(n1,no.grid,p))
  }
  
  X <- mu0 + V0 + 0.5 * (V0j0 + V0j1) + X
  Y <- mu1 + V1 + 0.5 * (V1j0 + V1j1) + Y
  
  return(list(X=X, Y=Y))
}

phi.fun=function(t,k){
  if(k%%2==0)  
    return(sqrt(2) * sin((k-1)*pi*t))  
  else  
    return(sqrt(2)*  cos(k*pi*t))
}


make.basis <- function(K = 50, T0 = 0, T1 = 1, no.grid = 101){
  t.grid <- seq(T0, T1, length.out = no.grid)
  kk <- matrix(rep(c(1:K),no.grid), no.grid, K, byrow = TRUE)
  tt <- matrix(rep(t.grid, K), no.grid, K, byrow = FALSE)
  
  phi <- matrix(
    mapply(phi.fun, tt, kk),no.grid, K) 
  
  return(phi)
}


zeta.eta <- function(case, p, n0, n1, K = 50){
  if(case == 1){
    sd.temp <-  array(rep(1/(1:K),each = n0*p), dim = c(n0,p,K))
    zeta <- array(sapply(sd.temp, function(x){rnorm(1,sd = x)}), dim = c(n0,p,K))
    sd.temp <-  array(rep(1/(1:K),each = n1*p), dim = c(n1,p,K))
    eta <- array(sapply(sd.temp, function(x){rnorm(1,sd = x)}), dim = c(n1,p,K))
  }
  
  if(case == 2){
    sd.temp <-  array(rep(1/(1:K),each = n0*p), dim = c(n0,p,K)) * 
      aperm(array(rep(sqrt(1+sample(seq(0.5,3,0.5),p,replace = T)),each = n0*K), dim = c(n0,K,p)), c(1,3,2))
    zeta <- array(sapply(sd.temp, function(x){rnorm(1,sd = x)}), dim = c(n0,p,K))
    
    sd.temp <-  array(rep(1/(1:K),each = n1*p), dim = c(n1,p,K)) 
    eta <- array(sapply(sd.temp, function(x){rnorm(1,sd = x)}), dim = c(n1,p,K))
    
  }
  
  if(case == 3){
    sd.temp <-  array(rep(0.5/(1:K),each = n0*p), dim = c(n0,p,K)) * 
      aperm(array(rep(sqrt(1+sample(seq(0.5,3,0.5),p,replace = T)),each = n0*K), dim = c(n0,K,p)), c(1,3,2))
    zeta <- exp(array(sapply(sd.temp, function(x){rnorm(1,sd = x)}), dim = c(n0,p,K))) 
    - array(sapply(sd.temp, function(x){ exp(x^2/2) }), dim = c(n0,p,K))
    
    sd.temp <-  array(rep(0.5/(1:K),each = n1*p), dim = c(n1,p,K))
    eta <- exp(array(sapply(sd.temp, function(x){rnorm(1,sd = x)}), dim = c(n1,p,K))) 
    - array(sapply(sd.temp, function(x){ exp(x^2/2) }), dim = c(n0,p,K))
  }
  
  if(case == 4){
    aj <- sample(seq(3,6,1),p,replace = T)
    sd.temp<-  array(rep(1/(1:K),each = n0*p), dim = c(n0,p,K))
    temp.norm <- array(rep(1/rep(1,K),each = n0*p), dim = c(n0,p,K))
    shape.temp <- aperm(array(rep((0.2 * aj),each = n0*K), dim = c(n1,K,p)), c(1,3,2))
    zeta1 <- array(sapply(temp.norm, function(x){rnorm(1,sd = x)}), dim = c(n0,p,K))
    zeta2 <- array(sapply(shape.temp, function(x){rgamma(1,shape = x,rate = 1) - x}), dim = c(n0,p,K))
    zeta <- 0.5*(zeta1 + zeta2)*sd.temp
    
    sd.temp<-  array(rep(1/(1:K),each = n1*p), dim = c(n1,p,K))
    temp.norm <- array(rep(1/rep(1,K),each = n1*p), dim = c(n1,p,K))
    shape.temp <- array(0.1,dim = c(n1,p,K))
    eta1 <- array(sapply(temp.norm, function(x){rnorm(1,sd = x)}), dim = c(n1,p,K))
    eta2 <- array(sapply(shape.temp, function(x){rgamma(1,shape = x,rate = 1) - x}), dim = c(n1,p,K))
    eta <- 0.5*(eta1 + eta2)*sd.temp
  }
  
  # changes for case 4
  # 4.1 variance of gamma increase to 2,3,4,5,6 and 
  if(case == 4.1){
    aj <- sample(seq(3,6,1),p,replace = T)
    sd.temp<-  array(rep(1/(1:K),each = n0*p), dim = c(n0,p,K))
    temp.norm <- array(rep(1/rep(1,K),each = n0*p), dim = c(n0,p,K))
    shape.temp <- aperm(array(rep((aj),each = n0*K), dim = c(n1,K,p)), c(1,3,2))
    zeta1 <- array(sapply(temp.norm, function(x){rnorm(1,sd = x)}), dim = c(n0,p,K))
    zeta2 <- array(sapply(shape.temp, function(x){rgamma(1,shape = x,rate = 1) - x}), dim = c(n0,p,K))
    zeta <- 0.5*(zeta1 + zeta2)*sd.temp
    
    sd.temp<-  array(rep(1/(1:K),each = n1*p), dim = c(n1,p,K))
    temp.norm <- array(rep(1/rep(1,K),each = n1*p), dim = c(n1,p,K))
    shape.temp <- array(1,dim = c(n1,p,K))
    eta1 <- array(sapply(temp.norm, function(x){rnorm(1,sd = x)}), dim = c(n1,p,K))
    eta2 <- array(sapply(shape.temp, function(x){rgamma(1,shape = x,rate = 1) - x}), dim = c(n1,p,K))
    eta <- 0.5*(eta1 + eta2)*sd.temp
  }
  
  # 4.2 mixture gamma shape = 2 and mean 2 and -2, and mixture normal
  if(case == 4.2){
    shape.temp1 <- aperm(array(rep((2),each = n0*K), dim = c(n1,K,p)), c(1,3,2))
    shape.temp2 <- aperm(array(rep((2),each = n0*K), dim = c(n1,K,p)), c(1,3,2))
    sd.temp<-  array(rep(1/(1:K),each = n0*p), dim = c(n0,p,K))
    
    zeta1 <- array(sapply(shape.temp1, function(x){rgamma(1,shape = x,rate = 1)}), dim = c(n0,p,K))
    zeta2 <- array(sapply(shape.temp2, function(x){rgamma(1,shape = x,rate = 1) - 4}), dim = c(n0,p,K))
    
    ber.coef <- array(rbinom(n0*p*K,1,0.5), dim = c(n0,p,K))
    zeta <- (ber.coef*zeta1 + (1-ber.coef)* zeta2)*sd.temp
    
    sd.temp<-  array(rep(1/(1:K),each = n1*p), dim = c(n1,p,K))
    temp.norm1 <- array(rep(2,n1*p*K), dim = c(n1,p,K))
    temp.norm2 <- array(rep(-2,n1*p*K), dim = c(n1,p,K))
    shape.temp <- array(1,dim = c(n1,p,K))
    eta1 <- array(sapply(temp.norm1, function(x){rnorm(1,mean = x)}), dim = c(n1,p,K))
    eta2 <- array(sapply(temp.norm2, function(x){rnorm(1,mean = x)}), dim = c(n1,p,K))
    ber.coef <- array(rbinom(n1*p*K,1,0.5), dim = c(n1,p,K))
    eta <- (ber.coef*eta1 + (1-ber.coef)* eta2)*sd.temp
  }
  
  # 4.3 gamma bimodal and unimodal
  if(case == 4.3){
    shape.temp1 <- aperm(array(rep((2),each = n0*K), dim = c(n1,K,p)), c(1,3,2))
    shape.temp2 <- aperm(array(rep((2),each = n0*K), dim = c(n1,K,p)), c(1,3,2))
    sd.temp<-  array(rep(1/(1:K),each = n0*p), dim = c(n0,p,K))
    
    zeta1 <- array(sapply(shape.temp1, function(x){rgamma(1,shape = x,rate = 1)}), dim = c(n0,p,K))
    zeta2 <- array(sapply(shape.temp2, function(x){rgamma(1,shape = x,rate = 1) - 4}), dim = c(n0,p,K))
    
    ber.coef <- array(rbinom(n0*p*K,1,0.5), dim = c(n0,p,K))
    zeta <- (ber.coef*zeta1 + (1-ber.coef)* zeta2)*sd.temp
    
    sd.temp<-  array(rep(1/(1:K),each = n1*p), dim = c(n1,p,K))
    temp.norm <- array(rep(1/rep(1,K),each = n1*p), dim = c(n1,p,K))
    shape.temp <- array(2,dim = c(n1,p,K))
    eta <- array(sapply(shape.temp, function(x){rgamma(1,shape = x,rate = 1) - x}), dim = c(n1,p,K))
    eta <- eta*sd.temp
  }
  
  if(case == 5){
    sd.temp <-  array(rep(1/(sqrt(3)*(1:K)),each = n0*p), dim = c(n0,p,K))
    t.temp <- array(1,dim = c(n0,p,K))
    zeta <- array(sapply(t.temp, function(x){rt(1,df = 3)}), dim = c(n0,p,K)) * sd.temp
    
    sd.temp <-  array(rep(1/(sqrt(3)*(1:K)),each = n1*p), dim = c(n1,p,K))
    t.temp <- array(1,dim = c(n1,p,K))
    eta <- array(sapply(t.temp, function(x){rt(1,df = 3)}), dim = c(n1,p,K)) * sd.temp
    
  }
  
  if(case == 6){
    sd.temp <-  array(rep(1/(sqrt(3)*(1:K)),each = n0*p), dim = c(n0,p,K))
    t.temp <- array(1,dim = c(n0,p,K))
    zeta <- array(sapply(t.temp, function(x){rt(1,df = 3)}), dim = c(n0,p,K)) * sd.temp * 
      aperm(array(rep(sqrt(1+sample(seq(0.5,3,0.5),p,replace = T)),each = n0*K), dim = c(n0,K,p)), c(1,3,2))
    
    sd.temp <-  array(rep(1/(sqrt(3)*(1:K)),each = n1*p), dim = c(n1,p,K))
    t.temp <- array(1,dim = c(n1,p,K))
    eta <- array(sapply(t.temp, function(x){rt(1,df = 3)}), dim = c(n1,p,K)) * sd.temp
    
  }
  
  return(list(zeta = zeta, eta = eta))
}

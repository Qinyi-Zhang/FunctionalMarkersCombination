rm(list=ls())

library(emplik)
library(MASS)
library(fdapace)
library(pROC)
library(glmnet)
library(mnormt)

set.seed(2019)
source("./fROC_functions.R")

######## data generation #########
source("./fROC_gendata.R")

#sample size
n0 <- 150
n1 <- 150
#number of functional markers
p <- 5
#number of observations generated on each curve
grid <- 101
#number of real observations on each curve
grid.obs <- 20
#case
case <- 1

data <- gen.data(case, n0, n1, p, no.grid = grid, K=50)
X <- data$X[, (grid-grid.obs):grid, ] # non-diseased group
Y <- data$Y[, (grid-grid.obs):grid, ] # diseased group

######## seperate training and testing data ########
index0 <- sample(1:5, n0, replace = T)
index1 <- sample(1:5, n1, replace = T)
kc <- 1

n0.train <- sum(index0!=kc)
n1.train <- sum(index1!=kc)
n0.test <- sum(index0==kc)
n1.test <- sum(index1==kc)

X.train <- X[index0!=kc, , ]
Y.train <- Y[index1!=kc, , ]
X.test <- X[index0==kc, , ]
Y.test <- Y[index1==kc, , ]
label.test <- c(rep(0, n0.test), rep(1, n1.test))

######## obtain the projection scores ########
fve <- 99.99 #FVE
res.temp <- get.score(X.train,Y.train,FVE.th = fve, X.new = X.test,Y.new = Y.test)
zeta.hat <- res.temp$zeta.hat
eta.hat <- res.temp$eta.hat
zeta.hat.test <- res.temp$zeta.hat.test
eta.hat.test <- res.temp$eta.hat.test

######## derive the one-dimension scalar features ########
dist.temp <- get.square.dist.cv(zeta.hat, eta.hat, zeta.hat.test = zeta.hat.test, eta.hat.test = eta.hat.test)
D <- dist.temp$D
D.test <- dist.temp$D.test

###### derive combinations of scalar features by existing methods ######
break.tie <- runif(n0.test + n1.test) * 1e-8

# Su and Liu's linear combination
res.pl.test <- get.score.pl(c(rep(0,n0.train), rep(1,n1.train)), D, D.test)
score.pl <- res.pl.test$score.pl
res.pl <- roc(label.test, score.pl + break.tie )
se.pl <- res.pl$sensitivities
sp.pl <- res.pl$specificities

auc.pl <- res.pl$auc
yi.pl <- max(se.pl + sp.pl) - 1
coef.pl <- res.pl.test$coef.pl

# Liu et al.'s min-max combination 
res.mm.test <- get.score.min.max(c(rep(0,n0.train), rep(1,n1.train)), D, D.test)
score.mm <- res.mm.test$score.min.max
res.mm <- roc(label.test, score.mm + break.tie)
se.mm <- res.mm$sensitivities
sp.mm <- res.mm$specificities

auc.mm <- res.mm$auc 
yi.mm <- max(se.mm + sp.mm) - 1 
coef.mm <- res.mm.test$coef.min.max

# Chen et al.'s ELR linear combination
res.elr.test <- get.score.elr(c(rep(0,n0.train), rep(1,n1.train)), D, D.test)
score.elr <- res.elr.test$score.elr
res.elr <- roc(label.test, score.elr + break.tie)
se.elr <- res.elr$sensitivities
sp.elr <- res.elr$specificities

auc.elr <- res.elr$auc
yi.elr <- max(se.elr + sp.elr) - 1
coef.elr <- res.elr.test$coef.elr

# Kang et al.'s stepwise method 
res.step.test <- get.score.step(c(rep(0,n0.train), rep(1,n1.train)), D, D.test)
score.step <- res.step.test$score.step
res.step <- roc(label.test, score.step + break.tie)
se.step <- res.step$sensitivities
sp.step <- res.step$specificities

auc.step <- res.step$auc 
yi.step <- max(se.step + sp.step) - 1 
coef.step <- res.step.test$coef.step

# logistic regression
res.lg.test <- get.score.lg(c(rep(0,n0.train),rep(1,n1.train)), D, D.test)
score.lg <- res.lg.test$score.lg
res.lg <- roc(label.test, score.lg + break.tie)
se.lg <- res.lg$sensitivities
sp.lg <- res.lg$specificities

auc.lg <- res.step$auc 
yi.lg <- max(se.lg + sp.lg) - 1 
coef.lg <- res.lg.test$coef.lg

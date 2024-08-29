rm(list=ls())
    



#library(survival)
#library(quantreg)
#library(rankreg)
library(micEcon)

#file.dir <- ""
#file.dir <- "C:/Documents and Settings/xsong/My Documents/WangCY/grant/code/results/MCAR_EstimatePi/"
file.dir <- "C:/Documents and Settings/xsong/My Documents/WangCY/grant/code/results/MAR/MAR_2covs/MAR_EstimatePi/"
#file.dir <- "h:\\aft\\code\\"
#file.dir <- "C:\\Documents and Settings\\xsong\\My Documents\\aft\\code\\"
#aft.file <- paste(file.dir,"aftsp.r",sep="")
#source(aft.file)




mean.cut <- 1e20
se.cut <- 1e20

#set.seed(342)
#set.seed(4729)

n.simul <- 1000
n.boot <- 50
iter.max <- 100

n <- 500
valid.rate.true <- 0.36
error.type <- 2
#sigma.var <- 0.5^2  # error variance
#sigma.ic <- 0.25^2   # error variance for the instrumental variable
sigma.var <- 0.1  # error variance
sigma.ic <- 0.1   # error variance for the instrumental variable

obj.dat <- c(0,0,1,1,2,2)    # 0 incomplete, 1 complete, 2 all
corr.dat <- c(2,1,0,2,0,1)

n.cov <- 2
n.cov.ec <- 1                        # number of error contaminated covariate.
n.cov.ef <- n.cov-n.cov.ec           # number of error prone covariates.
beta.true <- c(0,-1,1)

for(n in c(1000,2000))
#for(n in c(4000))
{
for(valid.rate.true in c(0.36,0.64))
{
#for(error.type in c(1,2))
for(error.type in c(1))
{
for(sigma.var in c(0.1,0.2))
{
for(sigma.ic in c(0.1,0.2))
{

set.seed(4729)
tstart <- proc.time()

cat("start time:   ")
print(tstart)
cat("\n\n")

if(valid.rate.true==0.36)
{
    #b.true <- c(1.5,-1.8,0.4)
    #b.true <- c(-1.5,1.8,-0.4)
    b.true <- c(-1.4,1.5,-0.3)
    #b.true <- c(2,-2,0.5)
    #b.true <- c(3,-2,0.5)
    #b.true <- c(3,-2,0.5)
    #b.true <- c(3,-3,0.5)         # logit(P(R=1))=b.true*(1,Y,W,Z), valid.rate=0.3
    #b.true <- c(3,-4,1)         # logit(P(R=1))=b.true*(1,Y,W,Z), valid.rate=0.3
    #b.true <- c(3,-3,1)         # logit(P(R=1))=b.true*(1,Y,W,Z), valid.rate=0.3
    #b.true <- c(-1,0.8,0)         # logit(P(R=1))=b.true*(1,Y,W,Z), valid.rate=0.3
}else if(valid.rate.true==0.64)
{
    b.true <- c(1.4,-1.5,0.3)
}else if(valid.rate.true==0.63)
{
    b.true <- c(1.5,-1.8,0.4)
}else if(valid.rate.true==0.33)
{
    b.true <- c(-1,0.5,0.1)         # logit(P(R=1))=b.true*(1,Y,W,Z), valid.rate=0.3
    #b.true <- c(-1,0.8,0)         # logit(P(R=1))=b.true*(1,Y,W,Z), valid.rate=0.3
}else if(valid.rate.true==0.67)
{
    b.true <- c(1,-0.5,-0.1)
}else if(valid.rate.true==0.3)
{
    #b.true <- c(-1,0.5,0.1)           # logit(P(R=1))=b.true*(1,Y,W,Z), valid.rate=0.3
    b.true <- c(-2.3,2,1) 
}else if(valid.rate.true==0.4)
{
    b.true <- c(-1.6,2,1)          # logit(P(R=1))=b.true*(1,Y,W,Z), valid.rate=0.5
}else if(valid.rate.true==0.5)
{
    b.true <- c(-1.06,2,1)          # logit(P(R=1))=b.true*(1,Y,W,Z), valid.rate=0.5
}else if(valid.rate.true==0.7)
{
    b.true <- c(0.29,2,1)          # logit(P(R=1))=b.true*(1,Y,W,Z), valid.rate=0.7
}
#alpha.true <- c(0.5,2,1)
#alpha.true <- c(0,2,-1)
#alpha.true <- c(0,1,0)
#alpha.true <- c(-0.2,0.8,0.5)
alpha.true <- c(-0.2,0.8,0.5,0.6)

if(n.cov==1)
{
   X.mu <- 0
   X.var <- matrix(c(1),nrow=1)
}else
{
   X.mu <- rep(0,n.cov)
   X.var <- matrix(c(1,0.3,0.3,1),nrow=n.cov)
}
X.var.half <- t(chol(X.var))                             
sigma.var.half <- sqrt(sigma.var)

MY_SEP_ERROR=2
MY_MIXPROP_ERROR=0.3

error.type.string <- ifelse(error.type==1,"NormalError_", "MixedNormalError_")
error.type.string <- paste(as.character(n),"_",error.type.string,as.character(sigma.var),"_IV",as.character(sigma.ic),
    "_ValidRate",as.character(valid.rate.true),sep="")       
result.file <- paste(file.dir,error.type.string,".R.out",sep="")
sink(result.file,split=TRUE)
print(paste("++++++++++++++++++++ ",error.type.string," ++++++++++++++++++++"))
cat("\n")



# order of the covariates: error free covariate, error prone covariates
GenerateData <- function(n,beta.true)
{   
    n.var <- length(X.mu)
    delta <- rep(0,n)      # failure indicator for Cox model, not used for linear model and logistic model
    
    #browser()

    if(n.cov==1) # 1 cov    
    {
        X <- t(X.mu+X.var.half%*%matrix(rnorm(n.var*n,0,1),n.var,n))
        if(error.type==1) # normal
        {
            e <- rnorm(n,0,sigma.var.half)
        }
        else      # mixed normal
        {
            tmp <- rbinom(n,1,MY_MIXPROP_ERROR)
            m_dMixedSD = sqrt(sigma.var/(1+(1-(2*MY_MIXPROP_ERROR-1)*(2*MY_MIXPROP_ERROR-1))*MY_SEP_ERROR*MY_SEP_ERROR/4))
            m_dMixedMean = MY_SEP_ERROR/2*m_dMixedSD
            e <- rnorm(n,0,1)
            e = e*m_dMixedSD
            e =e+(tmp==1)*m_dMixedMean-(tmp==0)*m_dMixedMean
            
            e = e-(2*MY_MIXPROP_ERROR-1)*m_dMixedMean
        }
        
        W <- X+e  
        M <- alpha.true[1]*X^2+alpha.true[2]*X+alpha.true[3]+rnorm(n,0,sqrt(sigma.ic))
    }
    else # 2 covs    
    {
        X <- t(X.mu+X.var.half%*%matrix(rnorm(n.var*n,0,1),n.var,n))
        X[1,] <- (X[1,]<0.4)
        X2 <- as.matrix(X[,2],ncol=1)
        W <- X2+rnorm(n,0,sigma.var.half)  
        M <- alpha.true[1]*X2^2+alpha.true[2]*X2+alpha.true[3]+alpha.true[4]*X[1,]+rnorm(n,0,sqrt(sigma.ic))
    }
    
    X <- cbind(1,X)
    
    Xbeta <- X%*%beta.true    
    Y <- rbinom(n,1,1/(1+exp(-Xbeta)))      


    #r <- rbinom(n,1,valid.rate)
    pi.true <- c(1/(1+exp(-b.true[1]-b.true[2]*Y-b.true[3]*W)))
    #pi.true <- valid.rate.true
    r <- rbinom(n,1,pi.true)
    
    if(n.cov.ef>0)
    {
        W.ext <- cbind(1,X[,1+(1:n.cov.ef)],W)
        M.ext <- cbind(1,X[,1+(1:n.cov.ef)],M)
    }
    else
    {
        W.ext <- cbind(1,W)
        M.ext <- cbind(1,M)        
    }
    S <- cbind(W.ext,Y)
    dat <- list(X=X[,-1],W=W,M=M,Y=Y,r=r,W.ext=W.ext,M.ext=M.ext,S=S,pi.true=pi.true)

    return(dat)
    
}

GetEstimate.ideal <- function(dat)
{
    fit <- glm(Y~X, family=binomial(link="logit"),data=dat)
    beta.est <- fit$coef
    beta.se <- sqrt(diag(vcov(fit)))
    return(list(est=beta.est,se=beta.se))
}

GetEstimate.naive <- function(dat)
{
    fit <- glm(Y~W.ext-1, family=binomial(link="logit"),data=dat)
    beta.est <- fit$coef
    beta.se <- sqrt(diag(vcov(fit)))
    return(list(est=beta.est,se=beta.se))
}

# Get the matrix Sigma for the simple correction
GetSigma.sc <- function(beta.ini,dat)
{
    beta.n <- beta.ini[-2]
    beta.p <- beta.ini[-1]
    Sigma <- matrix(NA,2*(1+n.cov),2*(1+n.cov))
    if(!is.na(beta.ini[1]))
    {     
        #Sigma <- matrix(0,2*(1+n.cov),2*(1+n.cov))
        #for(i in 1:n)
        #{
        #    U.c.n <- (dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])
        #    U.c.p <- (dat$r[i]*(dat$Y[i]+(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])
        #    U.c <- c(U.c.n,U.c.p)
        #    Sigma <- Sigma+U.c%*%t(U.c)
        #}
        tmp.n <- (dat$Y-1+dat$Y*exp(-dat$W.ext%*%beta.n))
        G0.n.i <- tmp.n*dat$r/dat$pi.est
        tmp.p <- (dat$Y+(dat$Y-1)*exp(dat$W.ext%*%beta.p))
        G0.p.i <- tmp.p*dat$r/dat$pi.est
        
        U.c.n <- apply(dat$M.ext,2,'*',G0.n.i)  # n*(n.cov+n.cov.ec)
        U.c.p <- apply(dat$M.ext,2,'*',G0.p.i)        
        MS.ext <- matrix(apply(dat$S,2,"*",dat$M.ext),nrow=n) # n X ((1+n.cov)*(2+n.cov)) for the ith subject, dat$M.ext*t(dat$S) stored in the ith row
        
        U.der.b <- apply(MS.ext,2,"*",dat$r*(1-1/dat$pi.est)*tmp.n)
        U.der.b <- matrix(apply(U.der.b,2,mean),nrow=1+n.cov)
        U.der.b.omega.b <- t(apply(dat$omega.b,1,"%*%",t(U.der.b)))
        U.c.n <- U.c.n-U.der.b.omega.b
        
        U.der.b <- apply(MS.ext,2,"*",dat$r*(1-1/dat$pi.est)*tmp.p)
        U.der.b <- matrix(apply(U.der.b,2,mean),nrow=1+n.cov)
        U.der.b.omega.b <- t(apply(dat$omega.b,1,"%*%",t(U.der.b)))
        U.c.p <- U.c.p-U.der.b.omega.b
        
        U.c <- cbind(U.c.n,U.c.p)                 
        Sigma <- apply(matrix(apply(U.c,2,"*",U.c),nrow=length(dat$Y)),2,mean)
        Sigma <- matrix(Sigma,nrow=2*(1+n.cov))
            
       
    }
    return(Sigma)
}

GetEstimate.initial <- function(beta.start,dat)
{
    #beta.old <- beta.start
    #beta.new <- beta.start
    #beta.n <- rep(NA,n.cov+1)
    #beta.p <- rep(NA,n.cov+1)
    #iter <- 0
    #bSolution <- FALSE
    #while(iter<iter.max && !bSolution)
    #{
    #    U.c <- obj.c.n(beta.old,dat)
    #    U.c.der <- obj.c.n.der(beta.old,dat)
    #    U.c.der.inv <- solve(U.c.der)
    #    beta.new <- beta.old - U.c.der.inv%*%U.c
    #    bSolution <- all(abs(beta.new-beta.old)<1e-3)
    #    beta.old <- beta.new
    #    iter <- iter+1
    #}
    #
    #if(!bSolution)
    #    beta.n <- rep(NA,n.cov+1)
    #else
    #    beta.n <- beta.new
    #
    #iter <- 0 
    #bSolution <- FALSE   
    #while(iter<iter.max && !bSolution)
    #{
    #    U.c <- obj.c.p(beta.old,dat)
    #    U.c.der <- obj.c.p.der(beta.old,dat)
    #    beta.new <- beta.old - solve(U.c.der)%*%U.c
    #    bSolution <- all(abs(beta.new-beta.old)<1e-3)
    #    beta.old <- beta.new
    #    iter <- iter+1
    #}
    #
    #
    #if(!bSolution)
    #    beta.p <- rep(NA,n.cov+1)    
    #else
    #    beta.p <- beta.new
    #
    #beta.est <- rep(NA,2+n.cov)
    #if(!is.na(beta.n[1]) && !is.na(beta.p[1]))
    #    beta.est <- c(beta.n[1],beta.p[1],beta.n[2:(1+n.cov)])
    #else if(!is.na(beta.n[1]))
    #    beta.est <- c(beta.n[1],0,beta.n[2:(1+n.cov)])
    #else if(!is.na(beta.p[1]))
    #    beta.est <- c(0,beta.p)
    #    
       
    beta.est <- rep(NA,2+n.cov)
    beta.n <- rep(NA,1+n.cov)
    beta.p <- rep(NA,1+n.cov)
    beta.se.n <- rep(NA,1+n.cov)
    beta.se.p <- rep(NA,1+n.cov)

#browser()    
   
    fit <- optim(beta.start,obj.c.n,gr=obj.c.n.der,dat=dat,method="Nelder-Mead")
    if(fit$convergence==0 && fit$value < 1e-3)
    {
        beta.n <- fit$par 
        beta.se.n <- GetSe.n(beta.n,dat)
    }
    
    fit <- optim(beta.start,obj.c.p,gr=obj.c.p.der,dat=dat,method="Nelder-Mead")
    if(fit$convergence==0 && fit$value < 1e-3)
    {
        beta.p <- fit$par 
        beta.se.p <- GetSe.p(beta.p,dat)
    }  
    
    if(!is.na(beta.n[1]) && !is.na(beta.p[1]))
        beta.est <- c(beta.n[1],beta.p[1],(beta.n+beta.p)[2:(1+n.cov)]/2)
    #else if(!is.na(beta.n[1]))
    #    beta.est <- c(beta.n[1],0,beta.n[2:(1+n.cov)])
    #else if(!is.na(beta.p[1]))
    #    beta.est <- c(0,beta.p)
    
       
    

    return(list(beta.n=beta.n, beta.p=beta.p,beta.est=beta.est,beta.se.n=beta.se.n,beta.se.p=beta.se.p))
}

GetSe.n <- function(beta,dat)
{
    tmp <- (dat$Y-1+dat$Y*exp(-dat$W.ext%*%beta))
    G0.n.i <- tmp*dat$r/dat$pi.est
    #G0.p.i <- dat$Y+(dat$Y-1)*exp(dat$W.ext%*%beta.ini)
    
    U.c <- apply(dat$M.ext,2,'*',G0.n.i)  # n*(n.cov+n.cov.ec)
    MS.ext <- matrix(apply(dat$S,2,"*",dat$M.ext),nrow=n) # n X ((1+n.cov)*(2+n.cov)) for the ith subject, dat$M.ext*t(dat$S) stored in the ith row
        
    U.der.b <- apply(MS.ext,2,"*",dat$r*(1-1/dat$pi.est)*tmp)
    U.der.b <- matrix(apply(U.der.b,2,mean),nrow=1+n.cov)
    U.der.b.omega.b <- t(apply(dat$omega.b,1,"%*%",t(U.der.b)))
    U.c <- U.c-U.der.b.omega.b
    
    Sigma <- apply(matrix(apply(U.c,2,"*",U.c),nrow=length(dat$Y)),2,mean)
    Sigma <- matrix(Sigma,nrow=(1+n.cov)) 
    
    U.c.der <- obj.c.n.der(beta,dat)
    U.c.der.inv <- solve(U.c.der)
    beta.var <- U.c.der.inv%*%Sigma%*%t(U.c.der.inv)/n
    beta.se <- sqrt(diag(beta.var))
    return(beta.se)
}

GetSe.p <- function(beta,dat)
{
    tmp <- (dat$Y+(dat$Y-1)*exp(dat$W.ext%*%beta))
    G0.p.i <- tmp*dat$r/dat$pi.est
    #G0.p.i <- dat$Y+(dat$Y-1)*exp(dat$W.ext%*%beta.ini)
    
    U.c <- apply(dat$M.ext,2,'*',G0.p.i)  # n*(n.cov+n.cov.ec)
    MS.ext <- matrix(apply(dat$S,2,"*",dat$M.ext),nrow=n) # n X ((1+n.cov)*(2+n.cov)) for the ith subject, dat$M.ext*t(dat$S) stored in the ith row
        
    U.der.b <- apply(MS.ext,2,"*",dat$r*(1-1/dat$pi.est)*tmp)
    U.der.b <- matrix(apply(U.der.b,2,mean),nrow=1+n.cov)
    U.der.b.omega.b <- t(apply(dat$omega.b,1,"%*%",t(U.der.b)))
    U.c <- U.c-U.der.b.omega.b
    
    Sigma <- apply(matrix(apply(U.c,2,"*",U.c),nrow=length(dat$Y)),2,mean)
    Sigma <- matrix(Sigma,nrow=(1+n.cov)) 
    
    U.c.der <- obj.c.p.der(beta,dat)
    U.c.der.inv <- solve(U.c.der)
    beta.var <- U.c.der.inv%*%Sigma%*%t(U.c.der.inv)/n
    beta.se <- sqrt(diag(beta.var))
    return(beta.se)
}

obj.c.n <- function(beta,dat)
{    
    #U.c <- rep(0,1+n.cov)
    #for(i in 1:n)
    #{
    #    U.c <- U.c+t(dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])     
    #}
    
    G0.n.i <- (dat$Y-1+dat$Y*exp(-dat$W.ext%*%beta))*dat$r/dat$pi.est
    #G0.n.i <- (dat$Y-1+dat$Y*exp(-dat$W.ext%*%beta))*dat$r
    #G0.p.i <- dat$Y+(dat$Y-1)*exp(dat$W.ext%*%beta.ini)
    
    U.c <- apply(dat$M.ext,2,'*',G0.n.i)  # n*(n.cov+n.cov.ec)
    U.c <- apply(U.c,2,mean)  
    obj <- sum(U.c^2) 
    
    return(obj)
}

obj.c.n.der <- function(beta,dat)
{   
    #U.c.der <- matrix(0,1+n.cov,1+n.cov)
    #for(i in 1:n)
    #{
    #    U.c.der <- U.c.der-t(dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])%*%dat$W.ext[i,]
    #}
    G0.n.i <- dat$Y*exp(-dat$W.ext%*%beta)*dat$r/dat$pi.est
    WM.ext <- matrix(apply(dat$M.ext,2,"*",dat$W.ext),nrow=length(dat$Y))        # each row is W.ext.i%*%t(M.ext.i) stored in rows
    U.c.der <- -apply(WM.ext,2,"*",G0.n.i)    
    U.c.der <- apply(U.c.der,2,mean)
    U.c.der <- t(matrix(U.c.der,nrow=1+n.cov))
    
    return(U.c.der)    
}


obj.c.p <- function(beta,dat)
{    
    #U.c <- rep(0,1+n.cov)
    #for(i in 1:n)
    #{
    #    U.c <- U.c+t(dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])     
    #}
    
    G0.p.i <- (dat$Y+(dat$Y-1)*exp(dat$W.ext%*%beta))*dat$r/dat$pi.est
    
    U.c <- apply(dat$M.ext,2,'*',G0.p.i)  # n*(n.cov+n.cov.ec)
    U.c <- apply(U.c,2,mean)   
    obj <- sum(U.c^2) 
       
    return(obj)
    
}

obj.c.p.der <- function(beta,dat)
{   
    #U.c.der <- matrix(0,1+n.cov,1+n.cov)
    #for(i in 1:n)
    #{
    #    U.c.der <- U.c.der-t(dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])%*%dat$W.ext[i,]
    #}
    
    G0.p.i <- (dat$Y-1)*exp(dat$W.ext%*%beta)*dat$r/dat$pi.est
    WM.ext <- matrix(apply(dat$M.ext,2,"*",dat$W.ext),nrow=length(dat$Y))        # each row is W.ext.i%*%t(M.ext.i) stored in rows
    U.c.der <- apply(WM.ext,2,"*",G0.p.i)
    U.c.der <- apply(U.c.der,2,mean)
    U.c.der <- t(matrix(U.c.der,nrow=1+n.cov))
    
    return(U.c.der)    
}

# objective function for simple correction
obj.sc <- function(beta,dat,A)
{   
    #U.c <- rep(0,2*(1+n.cov))
    #for(i in 1:n)
    #{
    #    U.c.n <- (dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])
    #    U.c.p <- (dat$r[i]*(dat$Y[i]+(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])
    #    U.c <- U.c+c(U.c.n,U.c.p)
    #}    
      
    beta.n <- beta[-2]
    beta.p <- beta[-1]
    
    G0.n.i <- (dat$Y-1+dat$Y*exp(-dat$W.ext%*%beta.n))*dat$r/dat$pi.est
    G0.p.i <- (dat$Y+(dat$Y-1)*exp(dat$W.ext%*%beta.p))*dat$r/dat$pi.est
    
    tmp.n <- apply(dat$M.ext,2,'*',G0.n.i)  # n*(n.cov+n.cov.ec)
    tmp.p <- apply(dat$M.ext,2,'*',G0.p.i)
    U.c <- cbind(tmp.n,tmp.p)                 
    U.c <- apply(U.c,2,mean)   
    
    obj <- t(U.c)%*%A%*%U.c
    return(obj)
}

obj.sc.der <- function(beta,dat,A)
{  
    #U.c <- rep(0,2*(1+n.cov))
    #U.c.der <- matrix(0,2*(1+n.cov),1+n.cov)
    #for(i in 1:n)
    #{
    #    U.c.n <- (dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])
    #    U.c.p <- (dat$r[i]*(dat$Y[i]+(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])
    #    U.c <- U.c+c(U.c.n,U.c.p)
    #    
    #    U.c.n.der <- -(dat$r[i]*dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini)*dat$M.ext[i,])%*%dat$W.ext[i,]
    #    U.c.p.der <- (dat$r[i]*(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])%*%dat$W.ext[i,]
    #    U.c.der <- U.c.der+c(U.c.n.der,U.c.p.der)
    #}
    #
    #
    
    beta.n <- beta[-2]
    beta.p <- beta[-1]
    
    G0.n.beta.i <- dat$Y*exp(-dat$W.ext%*%beta.n)*dat$r/dat$pi.est
    G0.n.i <- (dat$Y-1)*dat$r/dat$pi.est+G0.n.beta.i
    G0.p.beta.i <- (dat$Y-1)*exp(dat$W.ext%*%beta.p)*dat$r/dat$pi.est
    G0.p.i <- dat$Y*dat$r/dat$pi.est+G0.p.beta.i
    
    tmp.n <- apply(dat$M.ext,2,'*',G0.n.i)  # n*(n.cov+n.cov.ec)
    tmp.p <- apply(dat$M.ext,2,'*',G0.p.i)
    U.c <- cbind(tmp.n,tmp.p)                 
    U.c <- apply(U.c,2,mean)
    
    WM.ext.n <- matrix(apply(dat$M.ext,2,"*",insertCol(dat$W.ext,2,0)),nrow=length(dat$Y))        # each row is W.ext.i%*%t(M.ext.i) stored in rows
    WM.ext.p <- matrix(apply(dat$M.ext,2,"*",insertCol(dat$W.ext,1,0)),nrow=length(dat$Y))        # each row is W.ext.i%*%t(M.ext.i) stored in rows
    G2.WM.n.ext.i <- -apply(WM.ext.n,2,"*",G0.n.beta.i)
    G2.WM.p.ext.i <- apply(WM.ext.p,2,"*",G0.p.beta.i)
        
    U.c.der <- cbind(G2.WM.n.ext.i,G2.WM.p.ext.i)
    U.c.der <- apply(U.c.der,2,mean)
    U.c.der <- t(matrix(U.c.der,nrow=2+n.cov))
    
    obj.der <- 2*t(U.c)%*%A%*%U.c.der
    return(obj.der)
}


# simple correction
GetEstimate.sc <- function(beta.start,dat)
{
#browser()
    beta.est <- rep(NA,2+n.cov) 
    beta.se <- rep(NA,2+n.cov) 
    
    fit <- GetEstimate.initial(beta.start,dat)
    beta.n <- fit$beta.n
    beta.p <- fit$beta.p
    beta.se.n <- fit$beta.se.n
    beta.se.p <- fit$beta.se.p
    beta.initial <- fit$beta.est
    if(!is.na(beta.initial[1]))
    {
        Sigma <- GetSigma.sc(beta.initial,dat)
        A <- try(solve(Sigma))
    
        
        
        if(!is.character(A[1])) 
        {
            #beta.est <- optim(beta.initial,obj.sc,gr=obj.sc.der,dat=dat,A=Sigma,method="BFGS")
            fit <- try(optim(beta.initial,obj.sc,gr=obj.sc.der,dat=dat,A=A,method="Nelder-Mead"))
            if(!is.character(fit[1]))
            {
                if(fit$convergence==0)
                {
                    beta.est <- fit$par
                    if(!is.na(beta.est[1]))
                    {
                        beta.se <- GetSe.sc(beta.est,dat,A)
                        if(is.na(beta.se[1])) beta.est <- rep(NA,2+n.cov) 
                    }
                   
                }
            }
        }
    }
    return(list(beta.n=beta.n, beta.p=beta.p, est=beta.est,se=beta.se,se.n=beta.se.n,se.p=beta.se.p))
    
}


GetSe.sc <- function(beta.est,dat,A)
{
    #B <- matrix(0,2*(1+n.cov),2*(1+n.cov))
    #U.c <- rep(0,2*(1+n.cov))
    #U.c.der <- matrix(0,2*(1+n.cov),1+n.cov)
    #for(i in 1:n)
    #{
    #    U.c.n <- t(dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])
    #    U.c.p <- t(dat$r[i]*(dat$Y[i]+(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])
    #    U.c <- U.c+rbind(U.c.n,U.c.p)
        
    #    U.c.n.der <- -t(dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])%*%dat$W.ext[i,]
    #    U.c.p.der <- t(dat$r[i]*(dat$Y[i]+(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])%*%dat$W.ext[i,]
    #    U.c.der <- U.c.der+rbind(U.c.n.der,U.c.p.der)
    #    B <- B+U.c%*%t(U.c)
    #}    
   
    beta.n <- beta.est[-2]
    beta.p <- beta.est[-1]
    
    G0.n.beta.i <- (dat$Y*exp(-dat$W.ext%*%beta.n))*dat$r/dat$pi.est
    G0.p.beta.i <- ((dat$Y-1)*exp(dat$W.ext%*%beta.p))*dat$r/dat$pi.est
    
    WM.ext.n <- matrix(apply(dat$M.ext,2,"*",insertCol(dat$W.ext,2,0)),nrow=length(dat$Y))        # each row is W.ext.i%*%t(M.ext.i) stored in rows
    WM.ext.p <- matrix(apply(dat$M.ext,2,"*",insertCol(dat$W.ext,1,0)),nrow=length(dat$Y))        # each row is W.ext.i%*%t(M.ext.i) stored in rows
    G2.WM.n.ext.i <- -apply(WM.ext.n,2,"*",G0.n.beta.i)
    G2.WM.p.ext.i <- apply(WM.ext.p,2,"*",G0.p.beta.i)
        
    U.c.der <- cbind(G2.WM.n.ext.i,G2.WM.p.ext.i)
    U.c.der <- apply(U.c.der,2,mean)
    U.c.der <- t(matrix(U.c.der,nrow=2+n.cov))
     
    B <- GetSigma.sc(beta.est,dat)   
   
    
    tmp <- try(solve(t(U.c.der)%*%A%*%U.c.der)) 
    beta.se <- rep(NA,2+n.cov)
    if(!is.character(tmp[1])) 
    {
        beta.var <- tmp%*%t(U.c.der)%*%A%*%B%*%A%*%U.c.der%*%tmp
        beta.se <- sqrt(diag(beta.var)/n)
    }
            
    
    
    return(beta.se)
}

GetCorrectionTerm.n <- function(beta.est,dat,corr.dat)
{
    #num <- rep(0,n.cov.ec)
    #for(i in 1:n)
    #{
    #    num <- num+(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$W[i,])    
    #}
    #den <- sum(dat$Y-1+dat$Y*exp(-dat$W.ext%*%beta.ini))
    #c.n <- num/den
    
    
    #beta.n <- beta.est[-2]
    beta.n <- beta.est
    

    weight <- rep(1,length(dat$Y))
    if(corr.dat==0) # incomplete
    {
        weight <- (1-dat$r)/(1-dat$pi.est)   
    }
    else if(corr.dat==1) # complete
    {
        weight <- dat$r/dat$pi.est
    }   
        
    G0.beta.n.i <- dat$Y*exp(-dat$W.ext%*%beta.n)*weight   
    G0.n.i <- (dat$Y-1)*weight+G0.beta.n.i 
    
    
    U.c <- apply(dat$W,2,'*',G0.n.i)  # n*(n.cov+n.cov.ec)
    U.c <- apply(U.c,2,mean)
    G0.beta <- apply(G0.beta.n.i,2,mean)
    c.n <- U.c/G0.beta
    
    return(c.n)
    
    
}


GetCorrectionTerm.p <- function(beta.est,dat,corr.dat)
{
    #num <- rep(0,n.cov.ec)
    #for(i in 1:n)
    #{
    #    num <- (dat$Y[i]+(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*dat$W[i,])        
    #}
    #den <- sum(dat$Y+(dat$Y-1)*exp(dat$W.ext%*%beta.ini))
    #c.p <- num/den
    
    #beta.p <- beta.est[-1]
    beta.p <- beta.est
    
    weight <- rep(1,length(dat$Y))
    if(corr.dat==0) # incomplete
    {
        weight <- (1-dat$r)/(1-dat$pi.est)   
    }
    else if(corr.dat==1) # complete
    {
        weight <- dat$r/dat$pi.est
    }     
    
    G0.beta.p.i <- (dat$Y-1)*exp(dat$W.ext%*%beta.p)*weight  
    G0.p.i <- dat$Y*weight+G0.beta.p.i
    
    U.c <- apply(dat$W,2,'*',G0.p.i)  # n*(n.cov+n.cov.ec)
    U.c <- apply(U.c,2,mean)
    G0.beta <- apply(G0.beta.p.i,2,mean)
    c.p <- U.c/G0.beta
    
    return(c.p)  
    
}


obj.ic <- function(beta,dat,A,c.n,c.p,obj.dat)
{   
    #U.c <- rep(0,2*(1+n.cov+n.cov.ec))
    #for(i in 1:n)
    #{
    #    tmp.n <- (dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*c(dat$M.ext[i,]),dat$W[i,]-c.n)
    #    tmp.p <- (dat$Y[i]+(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*c(dat$M.ext[i,]),dat$W[i,]-c.p)
    #    U.c <- U.c+c(tmp.n[1:(1+n.cov.ef)],
    #             tmp.p[1:(1+n.cov.ef)],
    #             r[i]*tmp.n[(2+n.cov.ef):(1+n.cov)]
    #             r[i]*tmp.p[(2+n.cov.ef):(1+n.cov)],
    #             (1-r[i])*tmp.n[(2+n.cov):(1+n.cov+n.cov.ec)],
    #             (1-r[i])*tmp.p[(2+n.cov):(1+n.cov+n.cov.ec)])        
    #}
    
    beta.n <- beta[-2]
    beta.p <- beta[-1]
    
    weight <- rep(1,length(dat$Y))
    if(obj.dat==0) # incomplete
    {
        weight <- (1-dat$r)/(1-dat$pi.est)
    }
    else if(obj.dat==1) # complete
    {
        weight <- dat$r/dat$pi.est
    }    
    
    G0.n.beta.i <- dat$Y*exp(-dat$W.ext%*%beta.n)
    G0.n.i <- dat$Y-1+G0.n.beta.i
    G0.p.beta.i <- (dat$Y-1)*exp(dat$W.ext%*%beta.p)
    G0.p.i <- dat$Y+G0.p.beta.i
    
    
    tmp.n <- apply(cbind(dat$M.ext,dat$W),2,'*',G0.n.i)  # n*(n.cov+n.cov.ec)
    tmp.p <- apply(cbind(dat$M.ext,dat$W),2,'*',G0.p.i)
    tmp.n.cr <- apply(t(c.n),2,'*',G0.n.beta.i)
    tmp.p.cr <- apply(t(c.p),2,'*',G0.p.beta.i)
    rho <- cbind(tmp.n[,1:(1+n.cov.ef)],
             tmp.p[,1:(1+n.cov.ef)],
             (tmp.n*dat$r/dat$pi.est)[,(2+n.cov.ef):(1+n.cov)],
             (tmp.p*dat$r/dat$pi.est)[,(2+n.cov.ef):(1+n.cov)],
             (tmp.n[,(2+n.cov):(1+n.cov+n.cov.ec)]-tmp.n.cr)*weight,
             (tmp.p[,(2+n.cov):(1+n.cov+n.cov.ec)]-tmp.p.cr)*weight
             )           
    
             
    U.c <- apply(rho,2,mean)
    
    
    obj <- t(U.c)%*%A%*%U.c
    return(obj)
}


obj.ic.der <- function(beta,dat,A,c.n,c.p)
{    
    #U.c <- rep(0,2*(1+n.cov+n.cov.ec))
    #U.c.der <- matrix(0,2*(1+n.cov+n.cov.ec),1+n.cov)
    #for(i in 1:n)
    #{
    #    tmp.n <- (dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*c(dat$M.ext[i,]),dat$W[i,]-c.n)
    #    tmp.p <- (dat$Y[i]+(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*c(dat$M.ext[i,]),dat$W[i,]-c.p)
    #    U.c <- U.c+c(tmp.n[1:(1+n.cov.ef)],
    #             tmp.p[1:(1+n.cov.ef)],
    #             r[i]*tmp.n[(2+n.cov.ef):(1+n.cov)]
    #             r[i]*tmp.p[(2+n.cov.ef):(1+n.cov)],
    #             (1-r[i])*tmp.n[(2+n.cov):(1+n.cov+n.cov.ec)],
    #             (1-r[i])*tmp.p[(2+n.cov):(1+n.cov+n.cov.ec)])        
    #             
    #    tmp.n.der <- -(dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*c(dat$M.ext[i,]),dat$W[i,]-c.n)%*%dat$W.ext[i,]
    #    tmp.p.der <- (dat$r[i]*(dat$Y[i]+(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*c(dat$M.ext[i,]),dat$W[i,]-c.p)%*%dat$W.ext[i,]
    #    U.c.der <- U.c.der+rbind(tmp.n[1:(1+n.cov.ef),],
    #             tmp.p[1:(1+n.cov.ef),],
    #             r[i]*tmp.n[(2+n.cov.ef):(1+n.cov),]
    #             r[i]*tmp.p[(2+n.cov.ef):(1+n.cov),],
    #             (1-r[i])*tmp.n[(2+n.cov):(1+n.cov+n.cov.ec),],
    #             (1-r[i])*tmp.p[(2+n.cov):(1+n.cov+n.cov.ec),]) 
    #    
    #}
    
    beta.n <- beta[-2]
    beta.p <- beta[-1]
    
    weight <- rep(1,length(dat$Y))
    if(obj.dat==0) # incomplete
    {
        weight <- (1-dat$r)/(1-dat$pi.est)   
    }
    else if(obj.dat==1) # complete
    {
        weight <- dat$r/dat$pi.est
    }     
    
   
    G0.n.beta.i <- dat$Y*exp(-dat$W.ext%*%beta.n)
    G0.n.i <- dat$Y-1+G0.n.beta.i
    G0.p.beta.i <- (dat$Y-1)*exp(dat$W.ext%*%beta.p)
    G0.p.i <- dat$Y+G0.p.beta.i
    
    tmp.n <- apply(cbind(dat$M.ext,dat$W),2,'*',G0.n.i)  # n*(n.cov+n.ec)
    tmp.p <- apply(cbind(dat$M.ext,dat$W),2,'*',G0.p.i)
    tmp.n.cr <- apply(t(c.n),2,'*',G0.n.beta.i)
    tmp.p.cr <- apply(t(c.p),2,'*',G0.p.beta.i)
    U.c <- cbind(tmp.n[,1:(1+n.cov.ef)],
             tmp.p[,1:(1+n.cov.ef)],
             (tmp.n*dat$r/dat$pi.est)[,(2+n.cov.ef):(1+n.cov)],
             (tmp.p*dat$r/dat$pi.est)[,(2+n.cov.ef):(1+n.cov)],
             (tmp.n[,(2+n.cov):(1+n.cov+n.cov.ec)]-tmp.n.cr)*weight,
             (tmp.p[,(2+n.cov):(1+n.cov+n.cov.ec)]-tmp.p.cr)*weight
             )           
       
    
    
    WM.ext.n <- matrix(apply(dat$M.ext,2,"*",insertCol(dat$W.ext,2,0)),nrow=length(dat$Y))        # each row is W.ext.i%*%t(M.ext.i) stored in rows
    WM.ext.p <- matrix(apply(dat$M.ext,2,"*",insertCol(dat$W.ext,1,0)),nrow=length(dat$Y))        # each row is W.ext.i%*%t(M.ext.i) stored in rows
    G2.WM.n.ext.i <- -apply(WM.ext.n,2,"*",G0.n.beta.i)
    G2.WM.p.ext.i <- apply(WM.ext.p,2,"*",G0.p.beta.i)
    
    WW.ext.n <- matrix(apply((t(t(dat$W)-c.n))*weight,2,"*",insertCol(dat$W.ext,2,0)),nrow=length(dat$Y))        # each row is W.ext.i%*%t(W-c) stored in rows
    WW.ext.p <- matrix(apply((t(t(dat$W)-c.p))*weight,2,"*",insertCol(dat$W.ext,1,0)),nrow=length(dat$Y))        # each row is W.ext.i%*%t(W-c) stored in rows
    G2.WW.n.ext.i <- -apply(WW.ext.n,2,"*",G0.n.beta.i)
    G2.WW.p.ext.i <- apply(WW.ext.p,2,"*",G0.p.beta.i)
        
    U.c.der <- cbind(G2.WM.n.ext.i[,1:((2+n.cov)*(1+n.cov.ef))],
                 G2.WM.p.ext.i[,1:((2+n.cov)*(1+n.cov.ef))],
                 (G2.WM.n.ext.i*dat$r/dat$pi.est)[,((2+n.cov)*(1+n.cov.ef)+1):((2+n.cov)*(1+n.cov))],
                 (G2.WM.p.ext.i*dat$r/dat$pi.est)[,((2+n.cov)*(1+n.cov.ef)+1):((2+n.cov)*(1+n.cov))],
                 (G2.WW.n.ext.i),
                 (G2.WW.p.ext.i))      
    
    
    
    U.c <- apply(U.c,2,mean)         
    
    U.c.der <- apply(U.c.der,2,mean)
    U.c.der <- t(matrix(U.c.der,nrow=(2+n.cov)))
    
    obj.der <-2*t(U.c)%*%A%*%U.c.der
    return(obj.der)
}

# improved correction
GetEstimate.ic <- function(beta.sc,beta.n,beta.p,dat,obj.dat,corr.dat)
{
#browser()
    c.n <- GetCorrectionTerm.n(beta.n,dat,corr.dat)
    c.p <- GetCorrectionTerm.p(beta.p,dat,corr.dat)
    A <- GetSigma.ic(beta.sc,dat,c.n,c.p,obj.dat,corr.dat)
    
    A <- try(solve(A))
    #if(is.character(A)) browser()

    
    beta.est <- rep(NA,2+n.cov) 
    beta.se <- rep(NA,2+n.cov) 
    #if(!is.na(A[1,1]))    
    if(!is.character(A))
    {
        #beta.est <- optim(beta.initial,obj.sc,gr=obj.sc.der,dat=dat,A=Sigma,method="BFGS")
        fit <- optim(beta.sc,obj.ic,gr=obj.ic.der,dat=dat,A=A,c.n=c.n,c.p=c.p,obj.dat=obj.dat,method="Nelder-Mead")
        if(fit$convergence==0)
        {
            beta.est <- fit$par
            if(!is.na(beta.est[1]))
            {
                beta.se <- GetSe.ic(beta.est,dat,A,c.n,c.p,obj.dat,corr.dat)
                if(is.na(beta.se[1])) beta.est <- rep(NA,2+n.cov) 
            }            
        }
    }
    return(list(est=beta.est,se=beta.se))
    
}


GetSe.ic <- function(beta.est,dat,A,c.n,c.p,obj.dat,corr.dat)
{
    #B <- matrix(0,2*(1+n.cov),2*(1+n.cov))
    #U.c <- rep(0,2*(1+n.cov))
    #U.c.der <- matrix(0,2*(1+n.cov),1+n.cov)
    #for(i in 1:n)
    #{
    #    U.c.n <- t(dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])
    #    U.c.p <- t(dat$r[i]*(dat$Y[i]+(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])
    #    U.c <- U.c+rbind(U.c.n,U.c.p)
        
    #    U.c.n.der <- -t(dat$r[i]*(dat$Y[i]-1+dat$Y[i]*exp(-dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])%*%dat$W.ext[i,]
    #    U.c.p.der <- t(dat$r[i]*(dat$Y[i]+(dat$Y[i]-1)*exp(dat$W.ext[i,]%*%beta.ini))*dat$M.ext[i,])%*%dat$W.ext[i,]
    #    U.c.der <- U.c.der+rbind(U.c.n.der,U.c.p.der)
    #    B <- B+U.c%*%t(U.c)
    #}
    
    
    beta.n <- beta.est[-2]
    beta.p <- beta.est[-1]  
    
    weight.obj <- rep(1,length(dat$Y))
    if(obj.dat==0) # incomplete
    {
        weight.obj <- (1-dat$r)/(1-dat$pi.est)   
    }
    else if(obj.dat==1) # complete
    {
        weight.obj <- dat$r/dat$pi.est
    } 
    
    
     
    G0.n.beta.i <- c(dat$Y*exp(-dat$W.ext%*%beta.n))
    G0.p.beta.i <- c((dat$Y-1)*exp(dat$W.ext%*%beta.p))
    WM.ext.n <- matrix(apply(dat$M.ext,2,"*",insertCol(dat$W.ext,2,0)),nrow=length(dat$Y))        # each row is W.ext.i%*%t(M.ext.i) stored in rows
    WM.ext.p <- matrix(apply(dat$M.ext,2,"*",insertCol(dat$W.ext,1,0)),nrow=length(dat$Y))        # each row is W.ext.i%*%t(M.ext.i) stored in rows
    G2.WM.n.ext.i <- -apply(WM.ext.n,2,"*",G0.n.beta.i)
    G2.WM.p.ext.i <- apply(WM.ext.p,2,"*",G0.p.beta.i)
    
    WW.ext.n <- matrix(apply(t(t(dat$W)-c.n),2,"*",insertCol(dat$W.ext,2,0)),nrow=length(dat$Y))        # each row is W.i%*%t(M.ext.i) stored in rows
    WW.ext.p <- matrix(apply(t(t(dat$W)-c.p),2,"*",insertCol(dat$W.ext,1,0)),nrow=length(dat$Y))        # each row is W.i%*%t(M.ext.i) stored in rows
    G2.WW.n.ext.i <- -apply(WW.ext.n,2,"*",G0.n.beta.i)
    G2.WW.p.ext.i <- apply(WW.ext.p,2,"*",G0.p.beta.i)
        
    U.c.der <- cbind(G2.WM.n.ext.i[,1:((2+n.cov)*(1+n.cov.ef))],
                 G2.WM.p.ext.i[,1:((2+n.cov)*(1+n.cov.ef))],
                 (G2.WM.n.ext.i*dat$r/dat$pi.est)[,((2+n.cov)*(1+n.cov.ef)+1):((2+n.cov)*(1+n.cov))],
                 (G2.WM.p.ext.i*dat$r/dat$pi.est)[,((2+n.cov)*(1+n.cov.ef)+1):((2+n.cov)*(1+n.cov))],
                 (G2.WW.n.ext.i*weight.obj),
                 (G2.WW.p.ext.i*weight.obj)) 
                 
    U.c.der <- apply(U.c.der,2,mean)
    U.c.der <- t(matrix(U.c.der,nrow=(2+n.cov)))
    
    B <- GetSigma.ic(beta.est,dat,c.n,c.p,obj.dat,corr.dat)
    
    tmp <- try(solve(t(U.c.der)%*%A%*%U.c.der)) 
    beta.se <- rep(NA,2+n.cov)
    if(!is.na(B) && !is.character(tmp[1])) 
    {
        beta.var <- tmp%*%t(U.c.der)%*%A%*%B%*%A%*%U.c.der%*%tmp
        beta.se <- sqrt(diag(beta.var)/n)
    }         
   
    
    return(beta.se)
}


# Get the matrix A for the improved correction
GetSigma.ic <- function(beta.sc,dat,c.n,c.p,obj.dat,corr.dat)
{
#browser()
    beta.n <- beta.sc[-2]
    beta.p <- beta.sc[-1]
    
    weight.obj <- rep(1,length(dat$Y))
    weight.obj.der <- 0
    if(obj.dat==0) # incomplete
    {
        weight.obj <- (1-dat$r)/(1-dat$pi.est)
        weight.obj.der <- weight.obj/(1-dat$pi.est)
    }
    else if(obj.dat==1) # complete
    {
        weight.obj <- dat$r/dat$pi.est
        weight.obj.der <- -weight.obj/dat$pi.est
    } 
    
    weight.corr <- rep(1,length(dat$Y))
    weight.corr.der <- 0
    if(corr.dat==0) # incomplete
    {
        weight.corr <- (1-dat$r)/(1-dat$pi.est)
        weight.corr.der <- weight.corr/(1-dat$pi.est)
    }
    else if(corr.dat==1) # complete
    {
        weight.corr <- dat$r/dat$pi.est
        weight.corr.der <- -weight.corr/dat$pi.est
    } 
    
    WW <- matrix(apply(dat$W,2,"*",dat$W),nrow=length(dat$Y))                    # each row is W.i%*%t(W.i) stored in rows
    WW.ext <- matrix(apply(dat$W,2,"*",dat$W.ext),nrow=length(dat$Y))            # each row is W.i%*%t(W.ext.i) stored in rows
    WM.ext <- matrix(apply(dat$M.ext,2,"*",dat$W.ext),nrow=length(dat$Y))        # each row is W.ext.i%*%t(M.ext.i) stored in rows
    
    # for the negative functions
    
    G0.n.beta.i <- c(dat$Y*exp(-dat$W.ext%*%beta.n))
    G0.n.i <- dat$Y-1+G0.n.beta.i
    G0.n.beta <- mean(G0.n.beta.i)
    G0.n.beta.obj <- mean(G0.n.beta.i*weight.obj)
    G0.n.beta.corr <- mean(G0.n.beta.i*weight.corr)
    G0.n <- mean(G0.n.i)    
    G0.n.corr <- mean(G0.n.i*weight.corr)    
    G1.n.beta.i <- apply(dat$W,2,"*",G0.n.beta.i)
    G1.n.i <- apply(dat$W,2,"*",G0.n.i)
    G1.n <- apply(G1.n.i,2,mean)
    G1.n.corr <- apply(G1.n.i*weight.corr,2,mean)
    G1.n.ext.i <- apply(dat$W.ext,2,"*",G0.n.i)
    G1.n.ext.beta.i <- apply(dat$W.ext,2,"*",G0.n.beta.i)
    # Calculation of G2    
    G2.n.i <- apply(WW,2,"*",G0.n.beta.i)    # each row is the elements of the matrix for the ith subject stored by rows
    G2.n.ext.i <- apply(WW.ext,2,"*",G0.n.beta.i) # each row is the elements of the matrix for the ith subject stored by rows
    G2.WM.n.ext.i <- apply(WM.ext,2,"*",dat$r/dat$pi.est*G0.n.beta.i)
    G2.WM.n.ext <- t(matrix(apply(G2.WM.n.ext.i,2,mean),nrow=1+n.cov))
    G2.WM.n.ext.inv <- try(solve(G2.WM.n.ext))
    
    Sigma <- NA
    
    if(!is.character(G2.WM.n.ext.inv[1]))
    {
        #G0.der.c.n <- -mean(G0.n.beta.i*weight.obj)
        G0.der.c.n <- -G0.n.beta.obj
        
        # for positive functions
        G0.p.beta.i <- c((dat$Y-1)*exp(dat$W.ext%*%beta.p))
        G0.p.i <- dat$Y+G0.p.beta.i
        G0.p.beta <- mean(G0.p.beta.i)
        G0.p.beta.obj <- mean(G0.p.beta.i*weight.obj)
        G0.p.beta.corr <- mean(G0.p.beta.i*weight.corr)
        G0.p <- mean(G0.p.i)    
        G0.p.corr <- mean(G0.p.i*weight.corr)   
        G1.p.beta.i <- apply(dat$W,2,"*",G0.p.beta.i)
        G1.p.i <- apply(dat$W,2,"*",G0.p.i)
        G1.p <- apply(G1.p.i,2,mean)
        G1.p.corr <- apply(G1.p.i*weight.corr,2,mean)
        G1.p.ext.i <- apply(dat$W.ext,2,"*",G0.p.i)
        G1.p.ext.beta.i <- apply(dat$W.ext,2,"*",G0.p.beta.i)
        # Calculation of G2    
        G2.p.i <- apply(WW,2,"*",G0.p.beta.i)    # each row is the elements of the matrix for the ith subject stored by rows
        G2.p.ext.i <- apply(WW.ext,2,"*",G0.p.beta.i) # each row is the elements of the matrix for the ith subject stored by rows
        G2.WM.p.ext.i <- apply(WM.ext,2,"*",dat$r/dat$pi.est*G0.p.beta.i)
        G2.WM.p.ext <- t(matrix(apply(G2.WM.p.ext.i,2,mean),nrow=1+n.cov))
        G2.WM.p.ext.inv <- try(solve(G2.WM.p.ext))
        if(!is.character(G2.WM.p.ext.inv[1]))
        {
        #G0.der.c.p <- -mean(G0.p.beta.i*weight.obj)
        G0.der.c.p <- -G0.p.beta.obj
    
    
    
        #for(i in 1:n)
        #{
        #    tmp.n <- G0.n.i[i]*c(dat$M.ext[i,]),dat$W[i,]-c.n)
        #    tmp.p <- G0.p.i[i]*c(dat$M.ext[i,]),dat$W[i,]-c.p)
        #    rho <- c(tmp.n[1:(1+n.cov.ef)],
        #             tmp.p[1:(1+n.cov.ef)],
        #             r[i]*tmp.n[(2+n.cov.ef):(1+n.cov)]
        #             r[i]*tmp.p[(2+n.cov.ef):(1+n.cov)],
        #             (1-r[i])*tmp.n[(2+n.cov):(1+n.cov+n.cov.ec)],
        #             (1-r[i])*tmp.p[(2+n.cov):(1+n.cov+n.cov.ec)])
           # more to input
            #ksi.n <- G1.n.i[i,]/G0.n-G1.n*G0.n.i[i]/G0.n^2
            #ita.n <- -G2.n.ext.i[i,]/G0.n.i-G1.n%*%t(t(G1.n.ext.i[i,]))
         #}
                     
        tmp.n <- apply(cbind(dat$M.ext,dat$W),2,'*',G0.n.i)  # n*(n.cov+n.cov.ec)
        tmp.p <- apply(cbind(dat$M.ext,dat$W),2,'*',G0.p.i)
        tmp.n.cr <- apply(t(c.n),2,'*',G0.n.beta.i)
        tmp.p.cr <- apply(t(c.p),2,'*',G0.p.beta.i)
        rho <- cbind(tmp.n[,1:(1+n.cov.ef)],
                 tmp.p[,1:(1+n.cov.ef)],
                 (tmp.n*dat$r/dat$pi.est)[,(2+n.cov.ef):(1+n.cov)],
                 (tmp.p*dat$r/dat$pi.est)[,(2+n.cov.ef):(1+n.cov)],
                 (tmp.n[,(2+n.cov):(1+n.cov+n.cov.ec)]-tmp.n.cr)*weight.obj,
                 (tmp.p[,(2+n.cov):(1+n.cov+n.cov.ec)]-tmp.p.cr)*weight.obj
                 )
                 
        
        MS.ext <- matrix(apply(dat$S,2,"*",dat$M.ext),nrow=n) # n X (n.cov.ec*(2+n.cov)) for the ith subject, dat$M.ext*t(dat$S) stored in the ith row
        WS <- matrix(apply(dat$S,2,"*",dat$W),nrow=n) # n X (n.cov.ec*(2+n.cov)) for the ith subject, dat$M.ext*t(dat$S) stored in the ith row    
        
        
        G1.b.n.2 <- apply(WS,2,"*",weight.corr.der*dat$pi.est*(1-dat$pi.est)*G0.n.i)
        G1.b.n.2 <- matrix(apply(G1.b.n.2,2,mean),nrow=n.cov.ec)
        G1.b.n.1 <- apply(weight.corr.der*dat$pi.est*(1-dat$pi.est)*apply(dat$S,2,"*",G0.n.i),2,mean)
        delta.n <- t(matrix(apply(dat$omega.b,1,"%*%",t(G1.b.n.2/G0.n.beta.corr-G1.n.corr%*%t(G1.b.n.1)/G0.n.beta.corr^2)),nrow=n.cov.ec))
        
        
        U.WM.der.b <- apply(MS.ext,2,"*",dat$r*(1-1/dat$pi.est)*G0.n.i)
        U.WM.der.b <- matrix(apply(U.WM.der.b,2,mean),nrow=1+n.cov)        
        U.WM.der.b.omega.b.n <- t(matrix(apply(dat$omega.b,1,"%*%",t(U.WM.der.b)),nrow=1+n.cov))
        omega.n <- apply(dat$M.ext,2,"*",G0.n.i)*dat$r/dat$pi.est-U.WM.der.b.omega.b.n     
        
        tmp.n.cr <- apply(t(c.n),2,'*',G0.n.beta.i)
        U.WW.der.b  <- apply(WS,2,"*",weight.obj.der*dat$pi.est*(1-dat$pi.est)*(G0.n.i-tmp.n.cr))
        U.WW.der.b <- matrix(apply(U.WW.der.b,2,mean),nrow=n.cov.ec)  
        U.WW.der.b.omega.b.n <- t(matrix(apply(dat$omega.b,1,"%*%",t(U.WW.der.b)),nrow=n.cov.ec))
        
        

             
        
    
        #browser()
        ksi.n <- (weight.corr*G1.n.i)/G0.n.beta.corr-
            kronecker(weight.corr*G0.n.beta.i,t(G1.n.corr))/G0.n.beta.corr^2
        #ita.n <- (-G2.n.ext.i/G0.n.beta+
        #        t(apply(G1.n.ext.beta.i,1,kronecker,t(G1.n)))/G0.n.beta^2)   # each row is the matrix for subject i stored by rows
        #ita.n <- t(matrix(t(ita.n),nrow=1+n.cov))     # each n.cov.ec rows are the n.cov.ec*(1+n.cov) matrix for the ith subject
        #ita.n <- t(matrix(apply(ita.n,1,'%*%',G2.WM.n.ext.inv),ncol=n)) # n*(n.cov.ec*(1+n.cov))
        G1.n.ext.beta <- apply(weight.corr*G1.n.ext.beta.i,2,mean)
        G2.n.ext <- t(matrix(apply(weight.corr*G2.n.ext.i,2,mean),nrow=1+n.cov))
        ita.n <- (-G2.n.ext/G0.n.beta.corr+G1.n.corr%*%t(G1.n.ext.beta)/G0.n.beta.corr^2)%*%(G2.WM.n.ext.inv)
        #tmp.ita <- ita.n     
        #ita.n <- c(t(ita.n))        
        #ita.n <- ita.n*t(rep(omega.n,n.cov.ec))    # n*(n.cov.ec*(1+n.cov))
        #ita.n <- matrix(ita.n,nrow=1+n.cov)           # (1+n.cov)*(n*n.cov.ec)
        #ita.n <- apply(ita.n,2,sum)
        #ita.n <- t(matrix(ita.n,nrow=n.cov.ec))     # n*(1+n.cov.ec)          
        ita.n <- t(matrix(apply(omega.n,1,"%*%",t(ita.n)),nrow=n.cov.ec))
             
        
        
        G1.b.p.2 <- apply(WS,2,"*",weight.corr.der*dat$pi.est*(1-dat$pi.est)*G0.p.i)
        G1.b.p.2 <- matrix(apply(G1.b.p.2,2,mean),nrow=n.cov.ec)
        G1.b.p.1 <- apply(weight.corr.der*dat$pi.est*(1-dat$pi.est)*apply(dat$S,2,"*",G0.p.i),2,mean)
        delta.p <- t(matrix(apply(dat$omega.b,1,"%*%",t(G1.b.p.2/G0.p.beta.corr-G1.p.corr%*%t(G1.b.p.1)/G0.p.beta.corr^2)),nrow=n.cov.ec))
        
        
        U.WM.der.b <- apply(MS.ext,2,"*",dat$r*(1-1/dat$pi.est)*G0.p.i)
        U.WM.der.b <- matrix(apply(U.WM.der.b,2,mean),nrow=1+n.cov)        
        U.WM.der.b.omega.b.p <- t(matrix(apply(dat$omega.b,1,"%*%",t(U.WM.der.b)),nrow=1+n.cov))
        omega.p <- apply(dat$M.ext,2,"*",G0.p.i)*dat$r/dat$pi.est-U.WM.der.b.omega.b.p 
        
        #browser() 
        
        tmp.p.cr <- apply(t(c.p),2,'*',G0.p.beta.i)
        U.WW.der.b  <- apply(WS,2,"*",weight.obj.der*dat$pi.est*(1-dat$pi.est)*(G0.p.i-tmp.p.cr))
        U.WW.der.b <- matrix(apply(U.WW.der.b,2,mean),nrow=n.cov.ec)  
        U.WW.der.b.omega.b.p <- t(matrix(apply(dat$omega.b,1,"%*%",t(U.WW.der.b)),nrow=n.cov.ec))
        
        
    
        #browser()
        ksi.p <- (weight.corr*G1.p.i)/G0.p.beta.corr-
            kronecker(weight.corr*G0.p.beta.i,t(G1.p.corr))/G0.p.beta.corr^2
        #ita.p <- (-G2.p.ext.i/G0.p.beta+
        #        t(apply(G1.p.ext.beta.i,1,kronecker,t(G1.p)))/G0.p.beta^2)   # each row is the matrix for subject i stored by rows
        #ita.p <- t(matrix(t(ita.p),nrow=1+n.cov))     # each n.cov.ec rows are the n.cov.ec*(1+n.cov) matrix for the ith subject
        #ita.p <- t(matrix(apply(ita.p,1,'%*%',G2.WM.p.ext.inv),ncol=n)) # n*(n.cov.ec*(1+n.cov))
        G1.p.ext.beta <- apply(weight.corr*G1.p.ext.beta.i,2,mean)
        G2.p.ext <- t(matrix(apply(weight.corr*G2.p.ext.i,2,mean),nrow=1+n.cov))
        ita.p <- (-G2.p.ext/G0.p.beta.corr+G1.p.corr%*%t(G1.p.ext.beta)/G0.p.beta.corr^2)%*%(G2.WM.p.ext.inv)
        #ita.p <- c(t(ita.p))        
        #ita.p <- ita.p*t(rep(omega.p,n.cov.ec))    # n*(n.cov.ec*(1+n.cov))
        #ita.p <- matrix(ita.p,nrow=1+n.cov)           # (1+n.cov)*(n*n.cov.ec)
        #ita.p <- apply(ita.p,2,sum)
        #ita.p <- t(matrix(ita.p,nrow=n.cov.ec))     # n*(1+n.cov.ec)          
        ita.p <- t(matrix(apply(omega.p,1,"%*%",t(ita.p)),nrow=n.cov.ec))            
          
              
        phi <- rho+
            cbind(matrix(0,n,2*(1+n.cov)),G0.der.c.n*(ksi.n+ita.n+delta.n),matrix(0,n,n.cov.ec))+
            cbind(matrix(0,n,2*(1+n.cov)),matrix(0,n,n.cov.ec),G0.der.c.p*(ksi.p+ita.p+delta.p))+
            cbind(matrix(0,n,2*(1+n.cov)),U.WW.der.b.omega.b.n,matrix(0,n,n.cov.ec))+
            cbind(matrix(0,n,2*(1+n.cov)),matrix(0,n,n.cov.ec),U.WW.der.b.omega.b.p)+
            cbind(matrix(0,n,2*(1+n.cov.ef)),U.WM.der.b.omega.b.n[,-(1:(1+n.cov.ef))],matrix(0,n,3*n.cov.ec))+
            cbind(matrix(0,n,2*(1+n.cov.ef)+n.cov.ec),U.WM.der.b.omega.b.p[,-(1:(1+n.cov.ef))],matrix(0,n,2*n.cov.ec))
            
        Sigma <- apply(matrix(apply(phi,2,"*",phi),nrow=length(dat$Y)),2,mean)
        Sigma <- matrix(Sigma,nrow=2*(1+n.cov+n.cov.ec))
        }
    }

    return(Sigma)
}



GetBootStrapSE <- function(dat)
{
    beta.bs.ideal <- matrix(0,1+n.cov,n.boot)
    beta.bs.nr <- matrix(0,1+n.cov,n.boot)
    beta.bs.sc <- matrix(0,2+n.cov,n.boot)
    beta.bs.sc.n <- matrix(0,1+n.cov,n.boot)
    beta.bs.sc.p <- matrix(0,1+n.cov,n.boot)
    beta.bs.ic <- array(0,c(2+n.cov,length(obj.dat),n.boot))   
    
#    browser()

    for(ns in 1:n.boot)
    {

        # Generate BootStrap Data
        index<-sample(1:n,n,replace=T)
        #dat.bs <- data.frame(dat)[index,]
        dat.bs <- list()
        dat.bs$X <- dat$X[index]
        dat.bs$W <- as.matrix(dat$W[index,],ncol=n.cov.ec)
        dat.bs$M <- as.matrix(dat$M[index,],ncol=n.cov.ec)
        dat.bs$Y <- dat$Y[index]
        dat.bs$r <- dat$r[index]
        dat.bs$W.ext <- dat$W.ext[index,]
        dat.bs$M.ext <- dat$M.ext[index,]
        dat.bs$S <- dat$S[index,]
        dat.bs$pi.true <- dat$pi.true[index]               
                
        fit <- GetEstimate.ideal(dat.bs)  
        beta.bs.ideal[,ns] <- fit$est
        fit <- GetEstimate.naive(dat.bs)
        beta.bs.nr[,ns] <- fit$est
        glm.b <- glm(r ~ Y + W, family=binomial,dat=dat.bs)
        b <- glm.b$coef
        pi.est <- 1/(1+exp(-b[1]-b[2]*dat.bs$Y-b[3]*dat.bs$W))
        dat.bs$pi.est <- c(pi.est)
        # Compute omega.b
        SST <- matrix(apply(dat.bs$S,2,"*",dat.bs$S),nrow=n)
        SST <- apply(SST,2,"*",dat.bs$pi.est*(1-dat.bs$pi.est))  
        S.b.der <- matrix(apply(SST,2,mean),nrow=2+n.cov) # S.b, score function for b
        S.b.der.inv <- solve(S.b.der)
        omega.b <- apply((dat.bs$r-dat.bs$pi.est)*dat.bs$S,1,"%*%",S.b.der.inv) # for the ith subject, t(S) X S.b.der stored in the ith column
        dat.bs$omega.b <- t(omega.b)
        fit <- GetEstimate.sc(beta.true,dat.bs)
        beta.bs.sc.n[,ns] <- fit$beta.n
        beta.bs.sc.p[,ns] <- fit$beta.p
        beta.bs.sc[,ns] <- fit$est
        if(!is.na(beta.bs.sc[1,ns]))
            beta.bs.start <- t(beta.bs.sc[,ns])
        else    
            beta.bs.start <- c(beta.true[1],beta.true)
        for(k in 1:length(obj.dat))
        {
            fit <- GetEstimate.ic(beta.bs.start,beta.bs.sc.n[,ns],beta.bs.sc.p[,ns],dat.bs,obj.dat[k],corr.dat[k])
            beta.bs.ic[,k,ns] <- fit$est
        }
    }
    
    # Compute SE
    
    #browser()
   
    #se.ideal <- apply(beta.bs.ideal,1,sd,na.rm=TRUE)
    #se.nr <- apply(beta.bs.nr,1,sd,na.rm=TRUE)
    #se.sc.n <- apply(beta.bs.sc.n,1,sd,na.rm=TRUE)
    #se.sc.p <- apply(beta.bs.sc.p,1,sd,na.rm=TRUE)
    #se.sc <- apply(beta.bs.sc,1,sd,na.rm=TRUE)
    #se.ic <- apply(beta.bs.ic,c(1,2),sd,na.rm=TRUE)
    
    #mad.ideal <- apply(beta.bs.ideal,1,sd,na.rm=TRUE)
    #mad.nr <- apply(beta.bs.nr,1,sd,na.rm=TRUE)
    #mad.sc.n <- apply(beta.bs.sc.n,1,sd,na.rm=TRUE)
    #mad.sc.p <- apply(beta.bs.sc.p,1,sd,na.rm=TRUE)
    #mad.sc <- apply(beta.bs.sc,1,sd,na.rm=TRUE)
    #mad.ic <- apply(beta.bs.ic,c(1,2),sd,na.rm=TRUE)
    
    ci.ideal <- apply(beta.bs.ideal,1,quantile,probs=c(0.025,0.975),na.rm=TRUE)
    ci.nr <- apply(beta.bs.nr,1,quantile,probs=c(0.025,0.975),na.rm=TRUE)
    ci.sc.n <- apply(beta.bs.sc.n,1,quantile,probs=c(0.025,0.975),na.rm=TRUE)
    ci.sc.p <- apply(beta.bs.sc.p,1,quantile,probs=c(0.025,0.975),na.rm=TRUE)
    ci.sc <- apply(beta.bs.sc,1,quantile,probs=c(0.025,0.975),na.rm=TRUE)
    ci.ic <- apply(beta.bs.ic,c(1,2),quantile,probs=c(0.025,0.975),na.rm=TRUE)
    
    
    
    
    #return(list(se.ideal=se.ideal,se.nr=se.nr,se.sc.n=se.sc.n,se.sc.p=se.sc.p,se.sc=se.sc,se.ic=se.ic,
    #    mad.ideal=mad.ideal,mad.nr=mad.nr,mad.sc.n=mad.sc.n,mad.sc.p=mad.sc.p,mad.sc=mad.sc,mad.ic=mad.ic)) 
    return(list(ci.ideal=ci.ideal,ci.nr=ci.nr,ci.sc.n=ci.sc.n,ci.sc.p=ci.sc.p,ci.sc=ci.sc,ci.ic=ci.ic))
}




# ideal estimator
beta.ideal <- matrix(NA,1+n.cov,n.simul)
beta.se.ideal <- matrix(NA,1+n.cov,n.simul)
beta.se.ideal.bs <- matrix(NA,1+n.cov,n.simul)
beta.mad.ideal.bs <- matrix(NA,1+n.cov,n.simul)  
beta.ci.ideal.bs <- array(NA,c(2,1+n.cov,n.simul))

# naive estimator    
beta.nr <- matrix(NA,1+n.cov,n.simul)
beta.se.nr <- matrix(NA,1+n.cov,n.simul)
beta.se.nr.bs <- matrix(NA,1+n.cov,n.simul)
beta.mad.nr.bs <- matrix(NA,1+n.cov,n.simul)   
beta.ci.nr.bs <- array(NA,c(2,1+n.cov,n.simul))

# simple correction estimator
beta.sc <- matrix(NA,2+n.cov,n.simul)
beta.se.sc <- matrix(NA,2+n.cov,n.simul)
beta.se.sc.bs <- matrix(NA,2+n.cov,n.simul)
beta.mad.sc.bs <- matrix(NA,2+n.cov,n.simul)  
beta.ci.sc.bs <- array(NA,c(2,2+n.cov,n.simul))   

beta.sc.n <- matrix(NA,1+n.cov,n.simul)
beta.se.sc.n <- matrix(NA,1+n.cov,n.simul)
beta.se.sc.n.bs <- matrix(NA,1+n.cov,n.simul)
beta.mad.sc.n.bs <- matrix(NA,1+n.cov,n.simul)   
beta.ci.sc.n.bs <- array(NA,c(2,1+n.cov,n.simul))   

beta.sc.p <- matrix(NA,1+n.cov,n.simul)
beta.se.sc.p <- matrix(NA,1+n.cov,n.simul)
beta.se.sc.p.bs <- matrix(NA,1+n.cov,n.simul)
beta.mad.sc.p.bs <- matrix(NA,1+n.cov,n.simul) 
beta.ci.sc.p.bs <- array(NA,c(2,1+n.cov,n.simul))   


# improved correction estimator
beta.ic <- array(NA,c(2+n.cov,length(obj.dat),n.simul))
beta.se.ic <- array(NA,c(2+n.cov,length(obj.dat),n.simul))
beta.se.ic.bs <- array(NA,c(2+n.cov,length(obj.dat),n.simul))
beta.mad.ic.bs <- array(NA,c(2+n.cov,length(obj.dat),n.simul))      
beta.ci.ic.bs <- array(NA,c(2,2+n.cov,length(obj.dat),n.simul))   

valid.rate <- 0


sink()
for(ns in 1:n.simul)
{
#browser()
    cat("+++++++++++++++++++++ simulation #", ns, "++++++++++++++++++++++++++++\n")
    dat <- GenerateData(n,beta.true)
    valid.rate <- valid.rate+sum(dat$r)
            
    #beta.ideal[ns] <- aft.fun(data$X,log(data$V),data$delta)$beta[2,]
    #beta.nr[ns] <- aft.fun(data$W,log(data$V),data$delta)$beta[2,]
    fit <- GetEstimate.ideal(dat)  
    beta.ideal[,ns] <- fit$est
    beta.se.ideal[,ns]  <- fit$se
    fit <- GetEstimate.naive(dat)
    beta.nr[,ns] <- fit$est
    beta.se.nr[,ns]  <- fit$se
    glm.b <- glm(r ~ S-1, family=binomial,data=dat)   
    b <- glm.b$coef
    pi.est <- 1/(1+exp(-dat$S%*%b))
    dat$pi.est <- c(pi.est)
    # Compute omega.b
    SST <- matrix(apply(dat$S,2,"*",dat$S),nrow=n)
    SST <- apply(SST,2,"*",dat$pi.est*(1-dat$pi.est))  
    S.b.der <- matrix(apply(SST,2,mean),nrow=2+n.cov) # S.b, score function for b
    S.b.der.inv <- solve(S.b.der)
    omega.b <- apply((dat$r-dat$pi.est)*dat$S,1,"%*%",S.b.der.inv) # for the ith subject, t(S) X S.b.der stored in the ith column
    dat$omega.b <- t(omega.b)
    fit <- GetEstimate.sc(beta.true,dat)
    beta.sc.n[,ns] <- fit$beta.n
    beta.sc.p[,ns] <- fit$beta.p
    beta.sc[,ns] <- fit$est
    beta.se.sc.n[,ns]  <- fit$se.n
    beta.se.sc.p[,ns]  <- fit$se.p
    beta.se.sc[,ns]  <- fit$se
    if(!is.na(beta.sc[1,ns]))
    {
        beta.start <- t(beta.sc[,ns])
    #else    
    #    beta.start <- c(beta.true[1],beta.true)
        for(k in 1:length(obj.dat))
        {
            fit <- GetEstimate.ic(beta.start,beta.sc.n[,ns],beta.sc.p[,ns],dat,obj.dat[k],corr.dat[k])
            beta.ic[,k,ns] <- fit$est
            beta.se.ic[,k,ns]  <- fit$se
        }
    }
    
    #fit <- GetBootStrapSE(dat) 
    #beta.ci.ideal.bs[,,ns] <- fit$ci.ideal
    #beta.ci.nr.bs[,,ns] <- fit$ci.nr
    #beta.ci.sc.n.bs[,,ns] <- fit$ci.sc.n
    #beta.ci.sc.p.bs[,,ns] <- fit$ci.sc.p
    #beta.ci.sc.bs[,,ns] <- fit$ci.sc
    #beta.ci.ic.bs[,,,ns] <- fit$ci.ic
    
    #beta.se.ideal.bs[,ns] <- fit$se.ideal
    #beta.se.nr.bs[,ns] <- fit$se.nr
    #beta.se.sc.n.bs[,ns] <- fit$se.sc.n
    #beta.se.sc.p.bs[,ns] <- fit$se.sc.p
    #beta.se.sc.bs[,ns] <- fit$se.sc
    #beta.se.ic.bs[,,ns] <- fit$se.ic
    
    #beta.mad.ideal.bs[,ns] <- fit$mad.ideal
    #beta.mad.nr.bs[,ns] <- fit$mad.nr
    #beta.mad.sc.n.bs[,ns] <- fit$mad.sc.n
    #beta.mad.sc.p.bs[,ns] <- fit$mad.sc.p
    #beta.mad.sc.bs[,ns] <- fit$mad.sc   
    #beta.mad.ic.bs[,,ns] <- fit$mad.ic  
    
}

sink(result.file,append=TRUE)


sd.ideal <- apply(beta.ideal,1,sd,na.rm=TRUE)
se.ideal <-  apply(beta.se.ideal,1,mean,na.rm=TRUE)

index <- 1:n.simul

est.ideal <- apply(beta.ideal[,index],1,mean,na.rm=TRUE)
sd.ideal <- apply(beta.ideal[,index],1,sd,na.rm=TRUE)
se.ideal <-  apply(beta.se.ideal[,index],1,mean,na.rm=TRUE)
cp.ideal <- apply(abs(beta.ideal[,index]-beta.true)<=qnorm(0.975)*beta.se.ideal[,index],1,mean,na.rm=TRUE)
med.ideal <- apply(beta.ideal[,index],1,median,na.rm=TRUE)
mad.ideal <- apply(beta.ideal[,index],1,mad,na.rm=TRUE)
med.se.ideal <-  apply(beta.se.ideal[,index],1,median,na.rm=TRUE)
se.ideal.bs <-  apply(beta.se.ideal.bs,1,mean,na.rm=TRUE)
mad.ideal.bs <-  apply(beta.mad.ideal.bs,1,mean,na.rm=TRUE)
#cp.ideal.bs <- apply(abs(beta.ideal[,index]-beta.true)<=qnorm(0.975)*beta.se.ideal.bs[,index],1,mean,na.rm=TRUE)
cp.ideal.bs <- apply(beta.ci.ideal.bs[1,,index]-beta.true<=0 & beta.ci.ideal.bs[2,,index]-beta.true>=0,1,mean,na.rm=TRUE)
converge.ideal <- n.simul-sum(is.na(beta.ideal[1,]))
#converge.ideal <- sum(index,na.rm=TRUE)
bias.med.ideal <- med.ideal-beta.true

sd.nr <- apply(beta.nr,1,sd,na.rm=TRUE)
est.nr <- apply(beta.nr[,index],1,mean,na.rm=TRUE)
sd.nr <- apply(beta.nr[,index],1,sd,na.rm=TRUE)
se.nr <-  apply(beta.se.nr[,index],1,mean,na.rm=TRUE)
cp.nr <- apply(abs(beta.nr[,index]-beta.true)<=qnorm(0.975)*beta.se.nr[,index],1,mean,na.rm=TRUE)
med.nr <- apply(beta.nr[,index],1,median,na.rm=TRUE)
mad.nr <- apply(beta.nr[,index],1,mad,na.rm=TRUE)
med.se.nr <-  apply(beta.se.nr[,index],1,median,na.rm=TRUE)
se.nr.bs <-  apply(beta.se.nr.bs,1,mean,na.rm=TRUE)
mad.nr.bs <-  apply(beta.mad.nr.bs,1,mean,na.rm=TRUE)
#cp.nr.bs <- apply(abs(beta.nr[,index]-beta.true)<=qnorm(0.975)*beta.se.nr.bs[,index],1,mean,na.rm=TRUE)
cp.nr.bs <- apply(beta.ci.nr.bs[1,,index]-beta.true<=0 & beta.ci.nr.bs[2,,index]-beta.true>=0,1,mean,na.rm=TRUE)
converge.nr <- n.simul-sum(is.na(beta.nr[1,]))
#converge.nr <- sum(index,na.rm=TRUE)
bias.med.nr <- med.nr-beta.true

sd.sc <- apply(beta.sc,1,sd,na.rm=TRUE)
est.sc <- apply(beta.sc[,index],1,mean,na.rm=TRUE)
sd.sc <- apply(beta.sc[,index],1,sd,na.rm=TRUE)
se.sc <-  apply(beta.se.sc[,index],1,mean,na.rm=TRUE)
cp.sc <- apply(abs(beta.sc[,index]-c(beta.true[1],beta.true))<=qnorm(0.975)*beta.se.sc[,index],1,mean,na.rm=TRUE)
med.sc <- apply(beta.sc[,index],1,median,na.rm=TRUE)
mad.sc <- apply(beta.sc[,index],1,mad,na.rm=TRUE)
med.se.sc <-  apply(beta.se.sc[,index],1,median,na.rm=TRUE)
se.sc.bs <-  apply(beta.se.sc.bs,1,mean,na.rm=TRUE)
mad.sc.bs <-  apply(beta.mad.sc.bs,1,mean,na.rm=TRUE)
#cp.sc.bs <- apply(abs(beta.sc[,index]-c(beta.true[1],beta.true))<=qnorm(0.975)*beta.se.sc.bs[,index],1,mean,na.rm=TRUE)
cp.sc.bs <- apply(beta.ci.sc.bs[1,,index]-c(beta.true[1],beta.true)<=0 & beta.ci.sc.bs[2,,index]-c(beta.true[1],beta.true)>=0,1,mean,na.rm=TRUE)
converge.sc <- n.simul-sum(is.na(beta.sc[1,]))
#converge.sc <- sum(index,na.rm=TRUE)
bias.med.sc <- med.sc-c(beta.true[1],beta.true)

sd.sc.n <- apply(beta.sc.n,1,sd,na.rm=TRUE)
est.sc.n <- apply(beta.sc.n[,index],1,mean,na.rm=TRUE)
sd.sc.n <- apply(beta.sc.n[,index],1,sd,na.rm=TRUE)
se.sc.n <-  apply(beta.se.sc.n[,index],1,mean,na.rm=TRUE)
cp.sc.n <- apply(abs(beta.sc.n[,index]-beta.true)<=qnorm(0.975)*beta.se.sc.n[,index],1,mean,na.rm=TRUE)
med.sc.n <- apply(beta.sc.n[,index],1,median,na.rm=TRUE)
mad.sc.n <- apply(beta.sc.n[,index],1,mad,na.rm=TRUE)
med.se.sc.n <-  apply(beta.se.sc.n[,index],1,median,na.rm=TRUE)
se.sc.n.bs <-  apply(beta.se.sc.n.bs,1,mean,na.rm=TRUE)
mad.sc.n.bs <-  apply(beta.mad.sc.n.bs,1,mean,na.rm=TRUE)
#cp.sc.n.bs <- apply(abs(beta.sc.n[,index]-beta.true)<=qnorm(0.975)*beta.se.sc.n.bs[,index],1,mean,na.rm=TRUE)
cp.sc.n.bs <- apply(beta.ci.sc.n.bs[1,,index]-beta.true<=0 & beta.ci.sc.n.bs[2,,index]-beta.true>=0,1,mean,na.rm=TRUE)
converge.sc.n <- n.simul-sum(is.na(beta.sc.n[1,]))
#converge.sc.n <- sum(index,na.rm=TRUE)
bias.med.sc.n <- med.sc.n-beta.true

sd.sc.p <- apply(beta.sc.p,1,sd,na.rm=TRUE)

est.sc.p <- apply(beta.sc.p[,index],1,mean,na.rm=TRUE)
sd.sc.p <- apply(beta.sc.p[,index],1,sd,na.rm=TRUE)
se.sc.p <-  apply(beta.se.sc.p[,index],1,mean,na.rm=TRUE)
cp.sc.p <- apply(abs(beta.sc.p[,index]-beta.true)<=qnorm(0.975)*beta.se.sc.p[,index],1,mean,na.rm=TRUE)
med.sc.p <- apply(beta.sc.p[,index],1,median,na.rm=TRUE)
mad.sc.p <- apply(beta.sc.p[,index],1,mad,na.rm=TRUE)
med.se.sc.p <-  apply(beta.se.sc.p[,index],1,median,na.rm=TRUE)
se.sc.p.bs <-  apply(beta.se.sc.p.bs,1,mean,na.rm=TRUE)
mad.sc.p.bs <-  apply(beta.mad.sc.p.bs,1,mean,na.rm=TRUE)
#cp.sc.p.bs <- apply(abs(beta.sc.p[,index]-beta.true)<=qnorm(0.975)*beta.se.sc.p.bs[,index],1,mean,na.rm=TRUE)
cp.sc.p.bs <- apply(beta.ci.sc.p.bs[1,,index]-beta.true<=0 & beta.ci.sc.p.bs[2,,index]-beta.true>=0,1,mean,na.rm=TRUE)
converge.sc.p <- n.simul-sum(is.na(beta.sc.p[1,]))
#converge.sc.p <- sum(index,na.rm=TRUE)
bias.med.sc.p <- med.sc.p-beta.true
#for(k in 1:length(obj.dat))
#{
#    est.ic <- t(apply(beta.ic,c(1,2),mean,na.rm=TRUE)) # by each dataset variation
#    sd.ic <- t(apply(beta.ic,c(1,2),sd,na.rm=TRUE))
#    se.ic <- t(apply(beta.se.ic,c(1,2),mean,na.rm=TRUE))
#    cp.ic <- t(apply(abs(beta.ic-c(beta.true[1],beta.true))<=qnorm(0.975)*beta.se.ic,c(1,2),mean,na.rm=TRUE))
#    converge.ic <- n.simul-apply(is.na(beta.ic[1,,]),1,sum)
#}



print(cbind(converge.ideal,beta.true,est.ideal,sd.ideal,se.ideal,med.ideal,bias.med.ideal,mad.ideal,med.se.ideal,cp.ideal))
print(cbind(converge.nr,beta.true,est.nr,sd.nr,se.nr,med.nr,bias.med.nr,mad.nr,med.se.nr,cp.nr))
print(cbind(converge.sc.n,beta.true,est.sc.n,sd.sc.n,se.sc.n,med.sc.n,bias.med.sc.n,mad.sc.n,med.se.sc.n,cp.sc.n))
print(cbind(converge.sc.p,beta.true,est.sc.p,sd.sc.p,se.sc.p,med.sc.p,bias.med.sc.p,mad.sc.p,med.se.sc.p,cp.sc.p))
print(cbind(converge.sc,c(beta.true[1],beta.true),est.sc,sd.sc,se.sc,med.sc,bias.med.sc,mad.sc,med.se.sc,cp.sc))


#print(cbind(converge.ideal,beta.true,est.ideal,sd.ideal,se.ideal,se.ideal.bs,med.ideal,mad.ideal,cp.ideal,cp.ideal.bs))
#print(cbind(converge.nr,beta.true,est.nr,sd.nr,se.nr,se.nr.bs,med.nr,mad.nr,cp.nr,cp.nr.bs))
#print(cbind(converge.sc.n,beta.true,est.sc.n,sd.sc.n,se.sc.n,se.sc.n.bs,med.sc.n,mad.sc.n,cp.sc.n,cp.sc.n.bs))
#print(cbind(converge.sc.p,beta.true,est.sc.p,sd.sc.p,se.sc.p,se.sc.p.bs,med.sc.p,mad.sc.p,cp.sc.p,cp.sc.p.bs))
#print(cbind(converge.sc,c(beta.true[1],beta.true),est.sc,sd.sc,se.sc,se.sc.bs,med.sc,mad.sc,cp.sc,cp.sc.bs))



for(i in 1:length(obj.dat))
{
    sd.ic <- apply(beta.ic[,i,],1,sd,na.rm=TRUE)
    est.ic <- apply(beta.ic[,i,index],1,mean,na.rm=TRUE) # by each dataset variation
    sd.ic <- apply(beta.ic[,i,index],1,sd,na.rm=TRUE)
    se.ic <- apply(beta.se.ic[,i,index],1,mean,na.rm=TRUE)
    cp.ic <- apply(abs(beta.ic[,i,index]-c(beta.true[1],beta.true))<=qnorm(0.975)*beta.se.ic[,i,index],1,mean,na.rm=TRUE)
    med.ic <- apply(beta.ic[,i,index],1,median,na.rm=TRUE)
    mad.ic <- apply(beta.ic[,i,index],1,mad,na.rm=TRUE)
    med.se.ic <-  apply(beta.se.ic[,i,index],1,median,na.rm=TRUE)
    se.ic.bs <-  apply(beta.se.ic.bs[,i,index],1,mean,na.rm=TRUE)
    mad.ic.bs <-  apply(beta.mad.ic.bs[,i,index],1,mean,na.rm=TRUE)
    #cp.ic.bs <- apply(abs(beta.ic[,i,index]-c(beta.true[1],beta.true))<=qnorm(0.975)*beta.se.ic.bs[,i,index],1,mean,na.rm=TRUE)
    cp.ic.bs <- apply(beta.ci.ic.bs[1,,k,index]-c(beta.true[1],beta.true)<=0 & beta.ci.ic.bs[2,,k,index]-c(beta.true[1],beta.true)>=0,1,mean,na.rm=TRUE)
    converge.ic <- n.simul-sum(is.na(beta.ic[1,i,index]))
    #converge.ic <- sum(index,na.rm=TRUE)
    bias.med.ic <- med.ic-c(beta.true[1],beta.true)
    obj.dat.i <- obj.dat[i]
    corr.dat.i <- corr.dat[i]
    print(cbind(obj.dat.i,corr.dat.i,converge.ic,c(beta.true[1],beta.true),
        est.ic,sd.ic,se.ic,med.ic,bias.med.ic,mad.ic,med.se.ic,cp.ic))
}

valid.rate <- valid.rate/(n*n.simul)
cat("valid.rate=",valid.rate,"\n")

tend <- proc.time()
cat("start time: ", tstart, "\n")
cat("end time:   ", tend, "\n")
cat("processing time: ")
print(tend-tstart)

sink()
}
}
}
}
}
#quit("no")

##################################################
# This script use simulations to investigate the approximation of
# Cov(Theta_A,Theta_B) as rp^2*SE_A*SE_B with rp the phenotypic 
# correlation and Theta_A and B the estimates of genetic variance 
# for traits A and B.
#
# Authors : Loic Yengo and Valentin Hivert
##################################################

library(MASS)
library(parallel)

# l(x) = l(x0) + l'(x0)(x-x0) + l"(x-x0)^2 + o(x-x0), x~N(x0,S)
# So near x0, MLE is approximately equal to
# l(x) = l(x0) + l'(x0)x - l'(x0)x0 + l"(x0)x^2 - 2l"(x0)*x0*x + l"(x0)x0^2
# l(x) = [l"(x0)]x^2 + (l'(x0)-2l"(x0)*x0)*x + l(x0)- l'(x0)x0 + l"(x0)x0^2
# x    = -0.5*(l'(x0)-2l"(x0)*x0)/[l"(x0)] = -0.5*[l'(x0)/l"(x0)] + x0
# l is a linear function of [y*]^2
# cov[(x1-x01)*(x2-x02)] = 0.25 * E[l1'(x01)l2'(x02)] / E[l1"(x01)l2"(x02)] ~ rp^2?

N   <- 5000
M   <- 5000
X   <- matrix(rbinom(N*M,2,0.5),nrow=N,ncol=M)
Z   <- sqrt(2)*(X-1)
G   <- tcrossprod(Z) / M
EG  <- eigen(G)
h21 <- 0.5 # heritability - Trait 1
h22 <- 0.3 # heritability - Trait 2
rg  <- 0.0 # genetic correlation
re  <- 0.5 # environmental correlation
rp  <- rg*sqrt(h21*h22) + re*sqrt((1-h21)*(1-h22))
rm(X);rm(G)
gc()
sim <- function(h21,h22,rg,re){
  B   <- mvrnorm(n=M,mu=c(0,0),Sigma=matrix(c(1,rg,rg,1),2,2))
  E   <- mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,re,re,1),2,2))
  b1  <- sqrt(h21/M)*B[,1]
  b2  <- sqrt(h22/M)*B[,2]
  e1  <- sqrt(1-h21)*E[,1]
  e2  <- sqrt(1-h22)*E[,2]
  g1  <- c(Z%*%b1)
  g2  <- c(Z%*%b2)
  y1  <- g1 + e1
  y2  <- g2 + e2
  Y   <- cbind(y1,y2)
  s2  <- EG$values
  Yt  <- crossprod(EG$vectors,Y)
  Yt2 <- Yt^2
  
  ## GREML likelihood
  ll_greml <- function(theta,x){
    Vg <- theta[1]
    Ve <- theta[2]
    Vy <- Ve + Vg*s2
    sum(-log(Vy)-x/Vy)
  }
  
  mdML1 <- optim(par=c(Vg=h21,Ve=1-h21),fn=function(theta) ll_greml(theta,Yt2[,1]),
                 control=list(fnscale=-1,maxit=1e4))
  mdML2 <- optim(par=c(Vg=h22,Ve=1-h22),fn=function(theta) ll_greml(theta,Yt2[,2]),
                 control=list(fnscale=-1,maxit=1e4))
  ml1   <- mdML1$par
  ml2   <- mdML2$par
  rp    <- rg*sqrt(h21*h22) + re*sqrt((1-h21)*(1-h22))
  return(c(ml1,ml2,Erp=rp,Orp=cor(y1,y2)))
}

B  <- 1000 # number of replicates
f  <- function(r) sim(h21=0.5,h22=0.5,rg=r/2,re=r/2)
rs <- seq(0,1,len=11)
rx <- rep(rs,each=B)
fx <- do.call("rbind",mclapply(rx,f,mc.cores=7))
mx <- aggregate(fx,by=list(rx),FUN=mean)
sx <- aggregate(fx,by=list(rx),FUN=sd)
rt <- sapply(rs,function(r){
  l <- which(rx==r)
  cor(fx[l,1],fx[l,3])
})


plot(mx[,"Erp"],mx[,"Orp"],pch=19,axes=FALSE,
     xlab="Expected phenotypic correlation",
     ylab="Observed phenotypic correlation",
     xlim=c(0.2,0.7),ylim=c(0.2,0.7))
axis(1);axis(2)
abline(a=0,b=1,col="coral1")
for(i in 1:length(rs)){
  segments(mx[i,"Erp"],mx[i,"Orp"]-sx[i,"Orp"],
           mx[i,"Erp"],mx[i,"Orp"]+sx[i,"Orp"],col="dodgerblue")
}

Erp <- mx[,"Erp"]^2
ymin <- -0.1
ymax <- 0.35

par(mfrow=c(1,1),mar=c(5.1, 4.1, 4.1, 2.1))
plot(Erp,rt,pch=19,axes=FALSE,xlim=c(ymin,ymax),ylim=c(ymin,ymax),
     xlab=bquote("Expected correlation between ML estimator of"~h[SNP]^2),
     ylab=bquote("Observed correlation between ML estimator of"~h[SNP]^2))
axis(1);axis(2)
abline(a=0,b=1,col="coral1")
for(i in 1:length(rs)){
  segments(Erp[i],rt[i]-1.96/sqrt(B),Erp[i],rt[i]+1.96/sqrt(B),col="dodgerblue")
}
legend(ymin,ymax,legend="95%CI",col="dodgerblue",lty=1,box.lty=0)

###################################################
#
# Same for Additive-by-Additive
#
N   <- 5000
M   <- 5000
X   <- matrix(rbinom(N*M,2,0.5),nrow=N,ncol=M)
Z   <- sqrt(2)*(X-1)
G   <- tcrossprod(Z) / M
EG  <- eigen(G**2)
h21 <- 0.5 # heritability - Trait 1
h22 <- 0.3 # heritability - Trait 2
rg  <- 0.0 # genetic correlation
re  <- 0.5 # environmental correlation
rp  <- rg*sqrt(h21*h22) + re*sqrt((1-h21)*(1-h22))
sim <- function(h21,h22,rg,re){
  B   <- mvrnorm(n=M,mu=c(0,0),Sigma=matrix(c(1,rg,rg,1),2,2))
  E   <- mvrnorm(n=N,mu=c(0,0),Sigma=matrix(c(1,re,re,1),2,2))
  b1  <- sqrt(h21/M)*B[,1]
  b2  <- sqrt(h22/M)*B[,2]
  e1  <- sqrt(1-h21)*E[,1]
  e2  <- sqrt(1-h22)*E[,2]
  g1  <- c(Z%*%b1)
  g2  <- c(Z%*%b2)
  y1  <- g1 + e1
  y2  <- g2 + e2
  Y   <- cbind(y1,y2)
  s2  <- EG$values
  Yt  <- crossprod(EG$vectors,Y)
  Yt2 <- Yt^2
  
  ## GREML likelihood
  ll_greml <- function(theta,x){
    Vg <- theta[1]
    Ve <- theta[2]
    Vy <- Ve + Vg*s2
    sum(-log(Vy)-x/Vy)
  }
  
  mdML1 <- optim(par=c(Vg=h21,Ve=1-h21),fn=function(theta) ll_greml(theta,Yt2[,1]),
                 control=list(fnscale=-1,maxit=1e4))
  mdML2 <- optim(par=c(Vg=h22,Ve=1-h22),fn=function(theta) ll_greml(theta,Yt2[,2]),
                 control=list(fnscale=-1,maxit=1e4))
  ml1   <- mdML1$par
  ml2   <- mdML2$par
  rp    <- rg*sqrt(h21*h22) + re*sqrt((1-h21)*(1-h22))
  return(c(ml1,ml2,Erp=rp,Orp=cor(y1,y2)))
}

B  <- 1000 # number of replicates
f  <- function(r) sim(h21=0.5,h22=0.5,rg=r/2,re=r/2)
rs <- seq(0,1,len=11)
rx <- rep(rs,each=B)
fx <- do.call("rbind",mclapply(rx,f,mc.cores=4))
mx <- aggregate(fx,by=list(rx),FUN=mean)
sx <- aggregate(fx,by=list(rx),FUN=sd)
rt <- sapply(rs,function(r){
  l <- which(rx==r)
  cor(fx[l,1],fx[l,3])
})


plot(mx[,"Erp"],mx[,"Orp"],pch=19,axes=FALSE,
     xlab="Expected phenotypic correlation",
     ylab="Observed phenotypic correlation",
     xlim=c(0,0.7),ylim=c(0,0.7))
axis(1);axis(2)
abline(a=0,b=1,col="coral1")
for(i in 1:length(rs)){
  segments(mx[i,"Erp"],mx[i,"Orp"]-sx[i,"Orp"],
           mx[i,"Erp"],mx[i,"Orp"]+sx[i,"Orp"],col="dodgerblue")
}

Erp <- mx[,"Erp"]^2
ymin <- -0.1
ymax <- 0.35

plot(Erp,rt,pch=19,axes=FALSE,xlim=c(ymin,ymax),ylim=c(ymin,ymax),
     xlab="Expected correlation between ML estimator of Vg",
     ylab="Observed correlation between ML estimator of Vg")
axis(1);axis(2)
abline(a=0,b=1,col="coral1")
for(i in 1:length(rs)){
  segments(Erp[i],rt[i]-1.96/sqrt(B),Erp[i],rt[i]+1.96/sqrt(B),col="dodgerblue")
}
legend(ymin,ymax,legend="95%CI",col="dodgerblue",lty=1,box.lty=0)
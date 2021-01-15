#######################################################################
#
# Supplemental Figure 11 : Single-replicate estimates from unlinked 
# markers simulations using true or estimated allele frequencies to 
# standardize the genotypes.
#
#######################################################################

library(parallel)
library(plyr)
####################################################
p  <- 0.5
N  <- 1000
M  <- 1000
h2 <- 0.5
beta=rnorm(M,mean=0,sd=sqrt(h2/M))

sim <- function(){
  x     <- matrix(rbinom(N*M,size=2,prob=p),nrow=N,ncol=M)
  px <- colMeans(x)/2
  xd = sapply(1:ncol(x), function(j){
    suppressMessages(plyr::mapvalues(x[,j],c(1,2),c(2*px[j],4*px[j]-2)))
  })
  
  sx    <- 2*px*(1-px)
  z     <- do.call("cbind",lapply(1:M,function(j) (x[,j]-2*px[j])/sqrt(sx[j])))
  zd    <- do.call("cbind",lapply(1:M,function(j) (xd[,j]-2*px[j]**2)/sx[j]))
  G     <- tcrossprod(z)/M
  GD    <- tcrossprod(zd)/M
  G2    <- G^2/mean(diag(G^2))
  vecG <- G[upper.tri(G)]
  vecG2 <- G2[upper.tri(G2)]
  vecGD <- GD[upper.tri(GD)]
  vG    <- var(vecG)
  vG2   <- var(vecG2)
  rg    <- c(cor(vecG,vecGD), cor(vecG,vecG2), cor(vecGD,vecG2))
  
  g     <- c(z%*%beta)
  e     <- rnorm(N,sd=sqrt(1-h2))
  y     <- g + e
  
  Y     <- tcrossprod(y)
  vecY  <- Y[upper.tri(Y)]
  X     <- cbind(1,vecG,vecG2)
  XX    <- crossprod(X)
  XY    <- crossprod(X,vecY)
  ref   <- c(solve(XX)%*%XY)     ## HE: Y~G + G^2
  HE_AA <- cov(vecY,vecG2) / vG2 ## HE: Y~G^2
  
  px <- colMeans(x)/2
  xd = sapply(1:ncol(x), function(j){
    suppressMessages(plyr::mapvalues(x[,j],c(1,2),c(2*px[j],4*px[j]-2)))
  })
  
  #True All. freq
  px <- rep(p,M)
  xd = sapply(1:ncol(x), function(j){
    suppressMessages(plyr::mapvalues(x[,j],c(1,2),c(2*px[j],4*px[j]-2)))
  })
  sx    <- 2*px*(1-px)
  z     <- do.call("cbind",lapply(1:M,function(j) (x[,j]-2*px[j])/sqrt(sx[j])))
  zd    <- do.call("cbind",lapply(1:M,function(j) (xd[,j]-2*px[j]**2)/sx[j]))
  G     <- tcrossprod(z)/M
  GD    <- tcrossprod(zd)/M
  G2    <- G^2/mean(diag(G^2))
  vecG <- G[upper.tri(G)]
  vecG2 <- G2[upper.tri(G2)]
  vecGD <- GD[upper.tri(GD)]
  vG    <- var(vecG)
  vG2   <- var(vecG2)
  rg.truefreq    <- c(cor(vecG,vecGD), cor(vecG,vecG2), cor(vecGD,vecG2))
  
  g     <- c(z%*%beta)
  e     <- rnorm(N,sd=sqrt(1-h2))
  y     <- g + e
  
  Y     <- tcrossprod(y)
  vecY  <- Y[upper.tri(Y)]
  X     <- cbind(1,vecG,vecG2)
  XX    <- crossprod(X)
  XY    <- crossprod(X,vecY)
  ref.truefreq   <- c(solve(XX)%*%XY)     ## HE: Y~G + G^2
  HE_AA.truefreq <- cov(vecY,vecG2) / vG2 ## HE: Y~G^2
  
  
  return(c(rg,ref,HE_AA,rg.truefreq,ref.truefreq,HE_AA.truefreq))
}

B=1000
res <- do.call("rbind",mclapply(1:B, function(k) sim(),mc.cores=7))

par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(6.1, 6.1, 4.1, 3.1))
plot(x=res[,6],y=res[,7],pch=16,xlab=bquote(widehat(eta)[SNP]^2~(AAA)),ylab=bquote(widehat(eta)[SNP]^2~(AA)),cex.lab=1.5)
points(x=res[,13],y=res[,14],pch=16,col="red")
points(0,0,pch = 3,cex=2,lwd=3)
abline(0,1)
legend(x="bottomright",legend = c("True","Estimated"),title="Allele frequencies",pch = 16,col = c("red","black"),bty="n",cex = 2)
dev.off()

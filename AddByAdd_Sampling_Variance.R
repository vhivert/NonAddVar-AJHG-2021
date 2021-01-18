#####################################################################################
# This script reproduce the Figure S12 and S13.
#
# Using simulated data, this script investigate the accuracy of the estimates of the
# GRM_AA eigenvalues' variance using Equation 3 of the Supplemental Note 1.
# We also plot the estimates of sampling variance of the additive-by-additive variance 
# under REML using the equation 4 and 6 of the Supplemental Note 2.
#
# Author: Valentin Hivert
# Date: 18/01/21
#####################################################################################

library(RColorBrewer)
palette(brewer.pal(10,"Paired"))
library(Matrix)

#Function to compute the variance of AA (Hadamard of G matrix) eigenvalues and the corresponding SE from G (variance of diagonal and off-diag elements)
VarEigAA<-function(X){ #X is the additive GRM
  X=as.matrix(X)
  N=ncol(X)
  varGii=var(diag(X))
  diag(X)="NA"
  varGik=var(as.vector(X[lower.tri(X)]),na.rm=T)
  
  vareigAA=2*varGii^2 + 4*varGii + 3*(N-1)*varGik^2
  vareigAA_approx=14*varGik+(3*N)*varGik^2
  
  se=sqrt(2/(N*vareigAA))
  se_Approx=sqrt(2/(N*vareigAA_approx))
  se_OLS=sqrt(2/(N^2 * 2*varGik^2))
  return(list(VarEig=vareigAA,ApproxVarEig=vareigAA_approx,se=se,se_Approx=se_Approx,se_OLS=se_OLS))
}

#Simulations
nrep=100 #number of replicate
eig_Hadamard=matrix("NA",nrow = 1000,ncol = nrep)
res=matrix(0,nrow = nrep,ncol = 5)
M=10000 #Number of markers
N=500   #Sample size
for (r in 1:nrep) {
  print(r)
  freq_MAF=runif(M,min = 0.01,max = 0.5)#Uniform all. freq with MAF=0.01
  W       <- do.call("cbind",lapply(1:M,function(j,p,N) rbinom(N,2,prob=p[j]),p=freq_MAF,N=N))
  for(j in 1:M){
    while(length(unique(W[,j]))==1){W[,j]=rbinom(N,2,prob=freq_MAF[j])} #get biallelic SNPs
  }
  W = scale(W)
  GRM = W%*%t(W)/ (M)
  eig_Hadamard[,r]=as.numeric(eigen(GRM*GRM)$values)
  res[r,]=unlist(VarEigAA(GRM))
}

TrueVarEigenAA=apply(eig_Hadamard,2,var)#True eigenvalues of GRM_AA
varEigenAA=res[,1]
varEigenAA_Approx=res[,2]

#Figure S12
par(mfrow=c(1,1),mar = c(4, 4.9, 2, 0.75) + 0.1)
plot(x = TrueVarEigenAA,y=varEigenAA,xlab = bquote(Var(lambda[(G.G)])),ylab = bquote("Estimated" ~ Var(lambda[(G.G)])),pch=19,ylim=c(min(c(varEigenAA,varEigenAA_Approx)),max(c(varEigenAA,varEigenAA_Approx))),col=1,cex.lab=2,cex=1.5)
abline(a=0,b=1,lty=2,lwd=1)
abline(v = mean(TrueVarEigenAA),lty=2,lwd=2)
abline(h = mean(varEigenAA),lty=2,col=1,lwd=2)

#Figure S13
par(mfrow=c(1,1),mar = c(4, 4.5, 2, 0.75) + 0.1)
VarLambda=14/M +3*N/M^2#Approximation using Equation 5
VarApprox=2/(N*VarLambda) #Approximation using Equation 6
hist((res[,3])**2,breaks = 10,col=1,main="",xlab = bquote(Var(widehat(eta)[SNP[REML]]^2)),cex.lab=1.5)
abline(v=mean((res[,3])**2),lty=2,lwd=2)
abline(v=VarApprox,lty=2,lwd=2,col="red")


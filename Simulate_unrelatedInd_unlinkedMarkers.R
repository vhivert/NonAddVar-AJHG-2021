#######################################################################
#
# This script simulate 100,000 genotype data for 35000 unrelated individuals 
# using uniform distribution of allele frequencies with MAF>1%.
# Data output are in ped/fam format.
#
# Author: Valentin Hivert
#######################################################################

makeID <- function(i,base="IID"){
  if(i<10)           return(paste0(base,"000",i))
  if(i>=10  & i<100)   return(paste0(base,"00",i))
  if(i>=100 & i<1000)   return(paste0(base,"0",i))
  if(i>=1000) return(paste0(base,i))
}
makeSNP <- function(i){
  if(i<10)           return(paste0("rs000",i))
  if(i>=10   & i<100)   return(paste0("rs00",i))
  if(i>=100  & i<1000)   return(paste0("rs0",i))
  if(i>=1000 & i<10000) return(paste0("rs0",i))
  if(i>=10000) return(paste0("rs",i))
}  
simPlinkTextData <- function(N,  # Sample size 
                             M,  # Number of unliked markers
                             prefixOutput){
  ## ped/map file
  chr     <- rep(1,M)
  pos     <- sort(sample(1:1e6,M))
  a1a2    <- do.call("rbind",lapply(1:M,function(j) sample(c("A","C","G","T"),2)))
  snp     <- sapply(1:M,makeSNP)
  ## ped/geno
  mapData <- cbind.data.frame(chr,snp,0,pos)
  write.table(mapData,paste0(prefixOutput,".map"),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
  
  if(N*M>35000000){
    N=N/100
    freq=runif(M,min = 0.01,max = 0.99)
    
    for(i in 1:100){
      print(paste0("Part ",i,"/100"))
      
      X       <- do.call("cbind",lapply(1:M,function(j,p,N) rbinom(N,2,prob=p[j]),p=freq,N=N))
      
      for(j in 1:M){
        while(length(unique(X[,j]))==1){X[,j]=rbinom(N,2,prob=freq[j])}
      }
      
      refGeno <- t(sapply(1:M,function(j) c(paste0(a1a2[j,1],"\t",a1a2[j,1]),paste0(a1a2[j,1],"\t",a1a2[j,2]),paste0(a1a2[j,2],"\t",a1a2[j,2]))))
      ped     <- do.call("cbind",lapply(1:M,function(j) refGeno[j,1+X[,j]]))
      ## fam file
      iid    <- sapply(((i-1)*N+1):(i*N),makeID,"IID")
      fid    <- iid
      pid    <- rep(0,N)
      mid    <- rep(0,N)
      sex    <- sample(c(1,2),N,replace=TRUE)
      pheno  <- rep(-9,N)
      fam    <- cbind.data.frame(fid,iid,pid,mid,sex,pheno)
      if(i==1){
        write.table(cbind.data.frame(fam,ped),paste0(prefixOutput,".ped"),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
      }else{write.table(cbind.data.frame(fam,ped),paste0(prefixOutput,".ped"),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",append = T)}
    }
  }else{
    X       <- do.call("cbind",lapply(1:M,function(j) rbinom(N,2,prob=runif(1,min = 0.01,max = 0.99))))
    X       <- apply(X,2,function(j,N){
      while(length(unique(j))==1){j=rbinom(N,2,prob=runif(1))}
      return(j)
    },N=N)
    
    summary(colMeans(X)/2)
    refGeno <- t(sapply(1:M,function(j) c(paste0(a1a2[j,1],"\t",a1a2[j,1]),paste0(a1a2[j,1],"\t",a1a2[j,2]),paste0(a1a2[j,2],"\t",a1a2[j,2]))))
    ped     <- do.call("cbind",lapply(1:M,function(j) refGeno[j,1+X[,j]]))
    
    ## fam file
    iid    <- sapply(1:N,makeID,"IID")
    fid    <- iid
    pid    <- rep(0,N)
    mid    <- rep(0,N)
    sex    <- sample(c(1,2),N,replace=TRUE)
    pheno  <- rep(-9,N)
    fam    <- cbind.data.frame(fid,iid,pid,mid,sex,pheno)
    
    ## ped/geno
    mapData <- cbind.data.frame(chr,snp,0,pos)
    pedData <- cbind.data.frame(fam,ped)
    
    write.table(mapData,paste0(prefixOutput,".map"),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    write.table(pedData,paste0(prefixOutput,".ped"),quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
    
  }
  
}
simPlinkTextData(N = 35000,M = 100000,prefixOutput = "data_Simu_35K_ind_100K_unlinkedSNP")

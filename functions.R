### import all the functions needed ###
library(igraph)
library(network)
library(MASS)
library(xtable)
library(parallel)
library(foreach)
library(HHG)
#########################################
# disToRanks - compute a rank of a distance matrix
disToRanks <- function(dis) {
  # Transform from distance to ranking, order from 0,...,n-1.
  # For ties, the minimum ranking is used.
  n=nrow(dis);
  disRank=matrix(0,n,n);
  for (i in (1:n)){
    v=dis[,i];
    disRank[,i]=rank(v,ties.method="min")-1;
    #disRank[i,i]=0;
  }
  return(disRank);
}

# Local Graph Correlation
localGraphCorr <- function(X,Y,option, disRank){ # Calculate local graph correlation
  # Author: Cencheng Shen
  # Implements local graph correlation from Shen, Jovo, CEP 2016.
  # By specifying option=1, 2, or 3, it calculates the LGC statistics
  # by mcorr, dcorr, and Mantel
  if (missing(option)){
    option=1; # By default use mcorr
  }
  if (missing(disRank)){
    disRank=cbind(disToRanks(X), disToRanks(Y)); # Sort distances within columns, if the ranks are not provided
  }
  
  n=nrow(X);
  RX=disRank[,1:n]; # The ranks for X
  RY=disRank[,(n+1):(2*n)]; # The ranks for X
  corrXY=matrix(0,n,n);
  varX=rep(0,n);
  varY=rep(0,n);
  
  # Depending on the choice of global test, calculate the entries of A and B accordingly for late multiplication.
  if (option!=3){
    # Double centering for mcorr/dcorr
    H=diag(n)-(1/n)*matrix(1,n,n);
    A=H%*%X%*%H;
    B=H%*%Y%*%H;
    # For mcorr, further adjust the double centered matrices to remove high-dimensional bias
    if (option==1){
      A=A-X/n;
      B=B-Y/n;
      # The diagonals of mcorr are set to zero, instead of the original formulation of mcorr
      for (j in (1:n)){
        A[j,j]=0;
        B[j,j]=0;
      }
    }
  } else {
    # Single centering for Mantel
    EX=sum(X)/n/(n-1);
    EY=sum(Y)/n/(n-1);
    A=X-EX;
    B=Y-EY;
    # Mantel does not use diagonal entries, which is equivalent to set them zero
    for (j in (1:n)){
      A[j,j]=0;
      B[j,j]=0;
    }
  }
  
  # Summing up the entriwise product of A and B based on the ranks, which
  # yields the local family of covariance and variances
  for (j in (1:n)){
    for (i in (1:n)){
      a=A[i,j];
      b=B[i,j];
      # If there are ties, set all rank 0 entries to the diagonal entry
      if (RX[i,j]==0){
        a=A[j,j];
      }
      if (RY[i,j]==0){
        b=B[j,j];
      }
      tmp1=RX[i,j]+1;
      tmp2=RY[i,j]+1;
      corrXY[tmp1:n, tmp2:n]=corrXY[tmp1:n, tmp2:n]+a*b;
      varX[tmp1:n]=varX[tmp1:n]+a^2;
      varY[tmp2:n]=varY[tmp2:n]+b^2;
    }
  }
  
  # Normalizing the covariance by the variances yields the local family of correlation.
  options(warn=-1);
  corrXY=corrXY/Re(sqrt(varX%*%t(varY))); 
  options(warn=0);
  
  # Set any local correlation to 0 if any corresponding local variance is no larger than 0
  for (i in (1:n)){
    if (varX[i]<=0){
      corrXY[i,]=rep(0,n);
    }
    if (varY[i]<=0){
      corrXY[,i]=rep(0,n);
    }
  }
  result=list(corr=corrXY,varX=varX,varY=varY);
  return(result);
}

# Independence Test between two matrices?
IndependenceTest <-function(C,P,rep){
  # Author: Cencheng Shen
  # The independence Test for given data, an auxiliary function of the main
  # CorrPermDIstTest function
  if (missing(rep)){
    rep=1000; # By default use 1000 random permutations
  }
  
  n=nrow(C);
  alpha=0.05; # type 1 error level
  # Test statiscs under the null and the alternative
  dCor1N=array(0,dim=c(n,n,rep));dCor2N=array(0,dim=c(n,n,rep));dCor3N=array(0,dim=c(n,n,rep));
  dCor1A=array(0,dim=c(n,n,rep));dCor2A=array(0,dim=c(n,n,rep));dCor3A=array(0,dim=c(n,n,rep));
  # Powers for LGC by mcorr/dcorr/Mantel
  power1=matrix(0,n,n);power2=matrix(0,n,n);power3=matrix(0,n,n);
  
  for (r in (1:rep)){
    # Random sampling with replacement
    per=sample(n,replace=TRUE);
    Ca=C[per,per];
    Pa=P[per,per];
    disRank=cbind(disToRanks(Ca),disToRanks(Pa));
    # Calculate the test statistics under the alternative
    dCor1A[,,r]=localGraphCorr(Ca,Pa,1,disRank)$corr;
    dCor2A[,,r]=localGraphCorr(Ca,Pa,2,disRank)$corr;
    dCor3A[,,r]=localGraphCorr(Ca,Pa,3,disRank)$corr;
    
    # A different random sampling
    perN=sample(n,replace=TRUE);
    Pa=P[perN,perN];
    disRank=cbind(disToRanks(Ca),disToRanks(Pa));
    # Calculate the test statistics under the null
    dCor1N[,,r]=localGraphCorr(Ca,Pa,1,disRank)$corr;
    dCor2N[,,r]=localGraphCorr(Ca,Pa,2,disRank)$corr;
    dCor3N[,,r]=localGraphCorr(Ca,Pa,3,disRank)$corr;
  }
  
  # For each local test, estimate the critical value from the test statistics under the null,
  # then estimate the power from the test statistics under the alternative.
  for (i in (1:n)){
    for (j in (1:n)){
      dCorT=sort(dCor1N[i,j,],decreasing=TRUE);
      cut1=dCorT[ceiling(rep*alpha)];
      power1[i,j]=mean(dCor1A[i,j,]>cut1);
      
      dCorT=sort(dCor2N[i,j,],decreasing=TRUE);
      cut2=dCorT[ceiling(rep*alpha)];
      power2[i,j]=mean(dCor2A[i,j,]>cut2);
      
      dCorT=sort(dCor3N[i,j,],decreasing=TRUE);
      cut3=dCorT[ceiling(rep*alpha)];
      power3[i,j]=mean(dCor3A[i,j,]>cut3);
    }
  }
  
  # Set the powers of all local tests at rank 1 to 0 
  power1[1,]=rep(0,n);
  power2[1,]=rep(0,n);
  power3[1,]=rep(0,n);
  power1[,1]=rep(0,n);
  power2[,1]=rep(0,n);
  power3[,1]=rep(0,n);
  
  result=list(power1=power1,power2=power2,power3=power3);
  return(result);
}

# Permutation tests - for comparion between different distance correlation. 
PermutationTest <-function(C,P,rep,allP,option){
  # Author: Cencheng Shen
  # The permutation Test for given data, an auxiliary function of the main
  # CorrPermDIstTest function
  if (missing(rep)){
    rep=1000;
  }
  if (missing(allP)){
    allP=0; 
  }
  if (missing(option)){
    option=c(1,1,1,1); 
  }
  if (allP!=0){
    PALL=permn(n); 
    rep=length(PALL);
  }
  
  # P-values for LGC by mcorr, LGC by dcorr, LGC by Mantel, and HHG
  n=nrow(C);
  p1=matrix(0,n,n); p2=matrix(0,n,n);p3=matrix(0,n,n);p4=0;
  
  # Calculate the observed test statistics for the given data sets
  disRankC=disToRanks(C);
  disRankP=disToRanks(P);
  disRank=cbind(disRankC,disRankP);
  if (option[1]!=0){
    cut1=localGraphCorr(C,P,1,disRank)$corr;
  }
  if (option[2]!=0){
    cut2=localGraphCorr(C,P,2,disRank)$corr;
  }
  if (option[3]!=0){
    cut3=localGraphCorr(C,P,3,disRank)$corr;
  }
  if (option[4]!=0){
    cut4=hhg.test(C,P,nr.perm=0);
    cut4=unlist(cut4[1]);
  }
  
  # Now Permute the second dataset for rep times, and calculate the p-values
  for (r2 in (1:rep)){
    # Use random permutations; if allP is not 0, use all possible permutations
    per=sample(n);
    if (allP!=0){
      per=unlist(PAll[r2]);
    }
    Pa=P[per,per];
    disRank=cbind(disRankC,disRankP[per, per]);
    if (option[1]!=0){
      dCor1=localGraphCorr(C,Pa,1,disRank)$corr;
      p1=p1+(dCor1<cut1)*1/rep;
    }
    if (option[2]!=0){
      dCor2=localGraphCorr(C,Pa,2,disRank)$corr;
      p2=p2+(dCor2<cut2)*1/rep;
    }
    if (option[3]!=0){
      dCor3=localGraphCorr(C,Pa,3,disRank)$corr;
      p3=p3+(dCor3<cut3)*1/rep;
    }
    if (option[4]!=0){
      dCor4=hhg.test(C,Pa,nr.perm=0);
      dCor4=unlist(dCor4[1]);
      p4=p4+(dCor4<cut4)*1/rep;
    }
  }
  
  # Output the p-value
  p1=1-p1;
  p2=1-p2;
  p3=1-p3;
  p4=1-p4;
  
  result=list(LGCmcorr=p1,LGCdcorr=p2,LGCMantel=p3,HHG=p4,mcorr=p1[n,n],dcorr=p2[n,n],Mantel=p3[n,n]);
  return(result);
}

# verify neighbors
verifyNeighbors <-function(p,thres){
  # p is a matrix of pvalues for all family of tests?
  if (missing(thres)){
    thres=0; 
  }
  n=nrow(p);
  pmin=min(p);
  
  neighbor=which((p-pmin)<=thres);
  #neighbor=matrix(0,length(nind),2);
  #nind=nind[length(nind)];
  
  #   ind=((p-pmin)<thres);
  #   nind=which(ind==TRUE);
  #   nind=nind[length(nind)];
  #   
  #neighbor[,2]=ceiling(nind/n);
  #neighbor[,1]=nind-(neighbor[,2]-1)*(n);
  return(neighbor);
}

# Permutation tests for identifying dependency
CorrPermDistTest <- function(type,rep,cv,titlechar,allP,option){
  # Author: Cencheng Shen
  # Permutation Tests for identifying dependency.
  # The output are the p-values of LGC by mcorr/dcorr/Mantel, and HHG;
  # followed by the estimated optimal neighborhood for LGC by
  # mcorr/dcorr/Mantel.
  # Note that the local family include the global test at the last entry.
  
  # Parameters:
  # type should be a n*2n matrix, a concatenation of two distance matrices,
  # rep specifies the number of random permutations to use,
  # cv specifies the number of bootstrap samples to use for neighborhood validation,
  # set allP to non-zero will use all permutations instead,
  # option specifies whether each test statistic is calculated or not.
  if (missing(rep)){
    rep=1000; # Default random permutation numbers
  }
  if (missing(allP)){
    allP=0; # If set to other value, will override rep and use all permutations; unfeasible for n large
  }
  if (missing(titlechar)){
    titlechar="Data";
  }
  if (missing(cv)){
    cv=1000; # Default bootstrap replicates to estimate the optimal neighborhood 
  }
  if (missing(option)){
    option=c(1,1,1,1);  # Default option. Setting any to 0 to disable the calculation of mcorr/dcorr/Mantel/HHG.
  }
  
  n=nrow(type);
  C=type[, 1:n];
  P=type[, (n+1):(2*n)];
  
  # If cv is not 0, use resampling to estimate the optimal neighborhood by the testing powers
  ps1=matrix(0,n,n);ps2=matrix(0,n,n);
  if (cv!=0){
    testP=IndependenceTest(C,P,cv);
    neighbor1=verifyNeighbors(1-testP$power1);
    neighbor2=verifyNeighbors(1-testP$power2);
    neighbor3=verifyNeighbors(1-testP$power3);
  }
  # Return p-values from the permutation test
  testP=PermutationTest(C,P,rep,allP,option);
  if (cv==0){
    neighbor1=verifyNeighbors(testP$p1);
    neighbor2=verifyNeighbors(testP$p2);
    neighbor3=verifyNeighbors(testP$p3);
  }
  
  #   # Plot level plot
  #   max=0.2;p=testP$LGCmcorr;
  #   p[which(p>max)]=max;
  #   myAt=seq(0,max,0.02);
  #   interval=5;
  #   ckey=list(at=myAt,labels=list(cex=2));
  #   col.l <- colorRampPalette(c('red', 'orange', 'yellow', 'green', 'cyan', 'blue'))
  #   ckey=list(labels=list(cex=2));
  #   lplot=levelplot(p,zscaleLog="e",col.regions = terrain.colors(100),at=myAt,scales=list(x=list(at=seq(interval,n,interval), cex=2), y=list(at=seq(interval,n,interval), cex=2)),xlab=list(label="Neighborhood Choice of X",cex=2),ylab=list(label="Neighborhood Choice of Y",cex=2),main=list(label="Permutation Test P-Value",cex=2),colorkey=ckey)
  #   
  # Output
  #output=list(titlechar=titlechar,LGCmcorr=mean(testP$LGCmcorr[neighbor1]),LGCdcorr=mean(testP$LGCdcorr[neighbor2]),LGCMantel=mean(testP$LGCMantel[neighbor3]),HHG=testP$HHG, mcorr=testP$mcorr,dcorr=testP$dcorr,Mantel=testP$Mantel,n=n,rep=rep,allP=allP,option=option,levelPlot=lplot);
  output=list(titlechar=titlechar,LGCmcorr=testP$LGCmcorr,LGCdcorr=testP$LGCdcorr,LGCMantel=testP$LGCMantel,HHG=testP$HHG, mcorr=testP$mcorr,dcorr=testP$dcorr,Mantel=testP$Mantel,n=n,rep=rep,allP=allP,option=option,neighbormcorr=neighbor1,neighbordcorr=neighbor2,neighborMantel=neighbor3);
  
  filename = paste("CorrPermDistTestType", titlechar,".RData",sep = "");
  save(output,file = filename);
  return(output);
}

# univariate distance matrix
uni.dist <- function(vec){
  n <- length(vec)
  mat <- matrix(0, n, n)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      mat[i,j] <- sqrt((vec[i] - vec[j])^2)
      mat[j,i] <- mat[i,j]
    }
  }
  return(mat)
}

save.image("data/functions.RData")
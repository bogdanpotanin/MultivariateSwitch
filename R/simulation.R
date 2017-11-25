simulation <- function(n, x=NULL, nX=NULL, y.index, z.index, group, zo3, beta=NULL, alpha=NULL, sigmaX=NULL, sigma=NULL, betaRange=NULL, sigmaRange=c(0,1))
{
  if (is.null(nX)) nX=NCOL(x)
  #Number of selection equations
  ns=NCOL(zo3)
  #Number of groups
  ngroup=length(group)
  maxGroup=max(group)
  #Gegerating exogeneous variables
  if (is.null(sigmaX)) 
    {
      sigmaX=genPositiveDefMat(nX, rangeVar = sigmaRange, covMethod = "onion")$Sigma
      for (i in 1:nX) #Standartising
      {
        divider=sqrt(sigmaX[i,i])
        sigmaX[i,]=sigmaX[i,]/divider
        sigmaX[,i]=sigmaX[,i]/divider
      }
  }
  if (is.null(x)) x=mvrnorm(n,matrix(0,1,nX),sigmaX)
  #Generating distrurbances
  if (is.null(sigma))
  {
    sigma=genPositiveDefMat(ns+maxGroup, rangeVar = sigmaRange, covMethod = "onion")$Sigma
    for (i in 1:ns) #Standartising
    {
      divider=sqrt(sigma[i,i])
      sigma[i,]=sigma[i,]/divider
      sigma[,i]=sigma[,i]/divider
    }
  }
  disturbances=mvrnorm(n,matrix(0,1,ns+maxGroup),sigma)
  #Storing independend variables for z
  zh=matrix(list(), ns, 1)
  for (i in 1:ns)
  {
    zh[[i]]=cbind(1,x[,z.index[[i]]])
    colnames(zh[[i]])=c(0,z.index[[i]])
  }
  #Coefficients
  if (is.null(beta))
  {
    if (is.null(betaRange)) betaRange=c(-1.5,1.5)
    for (i in 1:maxGroup)
    {
      signs=runif(length(y.index[[i]])+1, min = -1, max = 1)
      signs[signs<0]=-1
      signs[signs>=0]=1
      beta[[i]]=signs*runif(length(y.index[[i]])+1, min = betaRange[1]*sqrt(sigma[ns+i-1,ns+i-1]), max = betaRange[2]*sqrt(sigma[ns+i-1,ns+i-1]))
    }
  }
  if (is.null(alpha))
  {
    for (i in 1:ns)
    {
      alpha[[i]]=runif(length(z.index[[i]])+1, min = -1, max = 1)
    }
  }
  #Gererating true values of z
  z0=matrix(ncol=ns, nrow=n)
  for (i in 1:ns)
  {
    z0[,i]=zh[[i]]%*%alpha[[i]]+disturbances[,i];
  }
  #Observable values of z0
  z=matrix(ncol=ns, nrow=n)
  for (i in 1:ns)
  {
    z[z0[,i]<0,i]=-1
    z[z0[,i]>=0,i]=1
  }
  #Correction for z selection
  groupMembership=matrix(FALSE,nrow=n,ncol=ngroup);
  for (i in 1:ngroup)
  {
    #Determining to which zo3 and outcome observation belongs
    #The problem is that if we have, for example 1 -1 0 0 and 1 -1 1 0 then double true could arise
    l=apply(z,1,function(x) all(x*abs(zo3[i,])==zo3[i,]))
    groupMembership[l,i]=TRUE
  }
  #Solving the problem of multiple TRUE lines
  for (i in 1:n)
  {
    trueInd=which(groupMembership[i,])
    if (length(trueInd)>1)
    {
    trueInd1=sample(trueInd,1)
    groupMembership[i,-trueInd1]=FALSE
    }
  }
  #Correcting z for ommitments
  for (i in 1:ngroup)
  {
    for (j in 1:ns)
    {
      z[groupMembership[,i] & zo3[i,j]==0,j]=0
    }
  }
  #Storing independend variables for y
  yh=matrix(list(), maxGroup, 1);
  for (i in 1:maxGroup)
  {
    yh[[i]]=cbind(1,x[,y.index[[i]]])
    colnames(yh[[i]])=c(0,y.index[[i]])
  }
  #Main variable
  y=matrix(NA,n,1);
  for (i in 1:ngroup)
  {
    if (group[i]!=0) y[groupMembership[,i],]=yh[[group[i]]][groupMembership[,i],]%*%beta[[group[i]]]+disturbances[groupMembership[,i],ns+group[i]]
  }
  return(list('y'=y, 'X'=x, 'yh'=yh,'z'=z,'zh'=zh,'group'=group,'zo3'=zo3,'beta'=beta,'alpha'=alpha,'sigma'=sigma,'sigmaX'=sigmaX))
}
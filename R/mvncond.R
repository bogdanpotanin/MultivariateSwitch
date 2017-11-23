mvncond<-function(mean, sigma, dependent.ind, given.ind, X.given)
{
  if (length(sigma)==1) {return(list(sigma,mean));}
  n=dim(sigma)[1];
  newOrder=c(dependent.ind,given.ind)
  mean=mean[,newOrder]
  sigma=covmatrixReplace(sigma,newOrder)
  ind=length(dependent.ind);
  s12s22=sigma[1:ind,(ind+1):n]%*%solve(sigma[(ind+1):n,(ind+1):n]);
  meanCond=mean[,1:ind]+(X.given-mean[,(ind+1):n])%*%t(s12s22);
  sigmaCond=sigma[1:ind,1:ind]-s12s22%*%sigma[(ind+1):n,1:ind];
  return(list("sigmaCond"=sigmaCond,"meanCond"=meanCond,"s12s22"=s12s22))
}
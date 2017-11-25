covmatrixReplace<-function(sigma,newOrder)
{
  if (length(sigma)==1) {return(sigma);}
  n=dim(sigma)[1];
  sigmaNew=matrix(0,n,n)
  if (length(unique(newOrder[newOrder>=1 & newOrder<=n]))!=n) {stop("newOrder vector length should match sigma dimension");}
  for (i in 1:n)
    {
      sigmaNew[i,]=sigma[newOrder[i],newOrder]
    }
  return(sigmaNew)
}

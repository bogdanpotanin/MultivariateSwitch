triangular<-function(triangularPart,diagonalPart)
{
  triangularPart=as.matrix(triangularPart)
  diagonalPart=as.matrix(diagonalPart)
  n=max(dim(diagonalPart));
  triang=matrix(ncol=n,nrow=n);
  a=1;
  for (i in 1:n)
  {
    a=a+(i>1)*(n-i+1);
    b=a+n-i-1;
    if (a<=b)
    {
      triang[i,]=matrix(c(rep(0,i),triangularPart[a:b]),nrow=1);
    }
    else
    {
      triang[i,]=0
    }
  }
  return(triang+as.matrix(Diagonal(n=n,x=diagonalPart))+t(triang));
}

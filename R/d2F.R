d2F<-function (x,Sigma,mu,rho=FALSE)
{
  #Do not use twice jk and kj
  n=dim(x)[1]
  ns=dim(Sigma)[1]
  d2FValue=array(0,c(n,ns,ns))
  if (ns<=1) return (d2FValue)
  if (mu==0) {mu=x*0;}
  x=x-mu;
  for (i in 1:ns)
  {
    listCond=mvncond(mean=matrix(0,n,ns), sigma=Sigma, dependent.ind=c((1:ns)[-i]), given.ind=i, X.given=x[,i])
    SC=listCond[[1]]
    meanCond=listCond[[2]]
    for (j in 1:(ns-1))
    {
      j1=j+(i<=j)
      listCond1=mvncond(mean=meanCond[,c((1:(ns-1))[-j],j),drop=FALSE], sigma=SC, dependent.ind=c((1:(ns-1))[-j]), given.ind=j, X.given=x[,j1])
      SC1=listCond1[[1]]
      meanCond1=listCond1[[2]]
      k1=dmvnorm(x[,c(i,j1)],c(0,0),Sigma[c(i,j1),c(i,j1)])
      if (ns>=3)
      {
        k2=pmnorm(x[,-c(i,j1)]-meanCond1, varcov = SC1)
      }
      else
      {
        k2=1;
      }
    d2FValue[,i,j1]=k1*k2;
    if (rho) {d2FValue[,i,j1]=d2FValue[,i,j1]*Sigma[i,j1]}
    }
  }
  return(d2FValue)
}
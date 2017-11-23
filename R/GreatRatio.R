dF<-function (x,Sigma,mu,denominator=FALSE)
{
#???????? ?? pmnorm
if (is.null(dim(x))) {n= length(x);}
else {n=dim(x)[1];}
ns=dim(Sigma)[1];
if (mu==0) {mu=x*0;}
GreatRatioValue=matrix(0,n,ns);
x=x-mu;
if (denominator) {Func=pmnorm(x, varcov=Sigma)}
else {Func=1;}
for (i in 1:ns)
{
listCond=mvncond(mean=matrix(0,n,ns), sigma=Sigma, dependent.ind=c((1:ns)[-i]), given.ind=i, X.given=x[,i]);
SC=listCond[[1]]
meanCond=listCond[[2]]
k1=dnorm(x[,i],0,sqrt(Sigma[i,i]));
if (ns>=2)
{
  k2=apply(X=x[,-i]-meanCond,MARGIN=1, FUN=function(q) pmvnorm(mean = 0, sigma = SC, lower=rep(-Inf,ns-1), upper=q));
}
else
{
  k2=1;
}
GreatRatioValue[,i]=k1*k2;
}
return(GreatRatioValue/Func)
}

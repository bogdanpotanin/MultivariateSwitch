#' Use this function to estimate Multivariate Switch model
#' @param data data frame containing the variables in the model
#' @param outcome main (continious) equation formula
#' @param selection1 the first selecton equation formula
#' @param selection2 the second selecton equation formula
#' @param selection3 the third selecton equation formula
#' @param selection4 the fourth selecton equation formula
#' @param selection5 the fifth selecton equation formula
#' @param group vector which determines outcome to each group
#' @param zo3 matrix which rows correspond to possible combination of selection equations values.
#' @param ShowInfo shows likelihood function optimization info
#' @param onlyTwostep if true then only two-step procedure used
#' @param opts options that are passed to nlopt
#' @param x0 optional initial values to be used while solving optimization task
#' @details This function estimates Multivariate Switch model via maximum-likelihood and two-step procedures.
#' This model was developed by Kossova E.V. and Potanin B.S.
#' Dependent variables in selection equations should have values -1 or 1.
#' Also there is a special value 0 which indicates that this observation is unobservable but it is necessary
#' to take into consideration other selection equation information while calculating likelihood function.
#' The i-th row of zo3 corresponds to i-th element of group. zo3 rows should contain information regarding
#' possible combinations of selection equations values. While group determines the outcome for each of this
#' possible combinations. Special 0 value for group responsible for sample selection.
gheckman<-function(data, outcome, selection1=NULL, selection2=NULL, selection3=NULL, selection4=NULL, selection5=NULL, group=NULL, zo3=NULL, ShowInfo=TRUE, onlyTwostep=FALSE, opts=list("algorithm" = "NLOPT_LD_TNEWTON", "xtol_rel" = 1e-16, "print_level" = 1, maxeval = 1000000), x0=NULL)
{
print("Version 1.0.0")
#PHASE 0: Extracting data from formulas
yh=model.frame(formula = outcome, data=data, na.action = NULL);#data for main equation
n=length(yh[,1]);#total number of observations
yh[rowSums(is.na(yh))>0,]=NA;#remove y that could not be calculated
y=yh[,1];#independend variable
yh[,1]=1;#intercept
colnames(yh)[1]="intercept";
#Data for selection equations
zh=matrix(list(),5,1)#selection equation independent variables
z=matrix(NA,n,5);#selection equation dedepndent variables
for (i in 1:5)#for each possible selection equation
{
  selection=paste("selection",toString(i),sep="");#name of current selection equation
  if (!is.null(get(selection)))#if this selection equation assigned
  {
    zh[[i]]=model.frame(formula = as.formula(get(selection)), data=data, na.action = NULL)
    zh[[i]][rowSums(is.na(zh[[i]]))>0,]=NA;#remove z that could not be calculated
    z[,i]=zh[[i]][,1];
    zh[[i]][,1]=rep(1,n);#intercept
    colnames(zh[[i]])[1]="intercept";
    zh[[i]][is.na(zh[[i]])]=0;
  }
  else
  {
    zh=zh[1:(i-1)];#preserve only existing selection equation
    z=z[,1:(i-1)];
    break
  }#if it is no such selection equation
}
#PHASE 1: Preparing data
sortList<-gheckmanSort(y, yh, z, zh, group, zo3);
y=sortList[[1]];
yh=sortList[[2]];
group=sortList[[5]];
zo3=sortList[[6]];
ngroup=sortList[[7]];
noutcome=sortList[[8]];
zo3Converter=sortList[[9]];
zo=sortList[[10]];
ns=sortList[[11]];
nSigma=sortList[[12]];
nsMax=sortList[[13]];
coef=sortList[[14]];
ndz=sortList[[15]];
groupsize=sortList[[16]];
nyh=sortList[[17]];
nzh=sortList[[18]];
nrhoY=sortList[[19]];
nrhoZ=sum(1:(nsMax-1));
rhoSigma=matrix(list(),noutcome);
parameters=vector(length = length(x0));#Store variables names
opts = opts;#setting optimization options max(maxeval/15,nsMax*50)
#PHASE 2: Two-step method
numberOfParameters=ndz[[nsMax]][dim(ndz[[nsMax]])[1]]
numberOfParametersY=coef[[noutcome]][dim(coef[[noutcome]])[1]]
x0=matrix(0,numberOfParameters);
#get initial values
if (nsMax==1) {z=matrix(z,ncol=1);}
for (i in 1:nsMax)
{
#-1 in order to exclude constant
x0[ndz[[i]]]=coef(myGlm<-glm(I((z[z[,i]!=0,i]+1)/2)~.-1,data=data.frame(zh[[i]][z[,i]!=0,]),family=binomial(link="probit")));
parameters[ndz[[i]]]=all.vars(formula(myGlm)[-2]);
}
zhColnames=matrix(list(),nsMax);#zhColnames
for (i in 1:nsMax) {zhColnames[[i]]=colnames(zh[[i]]);}
zh=sortList[[4]];#Only now we can load it from sort method result
z=sortList[[3]];#Only now we can load it from sort method result
#Set lower and upper bound constraints for parameters
lb=c(rep(-0.9999999,nrhoZ),rep(0,numberOfParametersY-nrhoZ),rep(-Inf,numberOfParameters-ndz[[1]][1]+1));
ub=c(rep(0.9999999,nrhoZ),rep(0,numberOfParametersY-nrhoZ),rep(Inf,numberOfParameters-ndz[[1]][1]+1));
#Estimate coefficients and store them to x0
f<-nloptr(x0=x0, eval_f=gheckmanLikelihood,opts=opts, lb=lb, ub=ub, y=y, zh=zh, yh=yh, zo=zo, ns=ns, ndz=ndz, nSigma=nSigma, coef=coef, group=group*0, ngroup=ngroup, nsMax=nsMax, zo3Converter=zo3Converter, noutcome=noutcome, zo3=zo3, groupsize=groupsize, nrhoY=nrhoY, ShowInfo=ShowInfo, maximization=FALSE);
x0=f$solution;
#x0=f$par;
#Storing covariance matrix
covmatrix=jacobian(func = gheckmanGradient, x = x0, y=y, zh=zh, yh=yh, zo=zo, ns=ns, ndz=ndz, nSigma=nSigma, coef=coef, group=group*0, ngroup=ngroup, nsMax=nsMax, zo3Converter=zo3Converter, noutcome=noutcome, zo3=zo3, groupsize=groupsize, nrhoY=nrhoY, ShowInfo=FALSE);
covmatrixGamma=solve(covmatrix[ndz[[1]][1]:length(x0),ndz[[1]][1]:length(x0)]);
#Calculating lambda
zTilde=matrix(list(), ngroup, 1);
lambda_z=matrix(list(), ngroup, 1);
lambda=matrix(list(), ngroup, 1);
lambda2=matrix(list(), ngroup, 1);
lambda_zTilde=matrix(list(), ngroup, 1);
LAMBDA=matrix(list(), ngroup, 1);
zG=matrix(list(), ngroup, 1);
zTildeG=matrix(list(), ngroup, 1);
Sigma=matrix(list(), ngroup, 1);
Sigma0=triangular(x0[1:nrhoZ],rep(1,nsMax));
for (i in 1:ngroup)
{
  Sigma[[i]]=Sigma0[zo3[i,]!=0,zo3[i,]!=0];
  zTilde[[i]]=z[[i]]*NA;
  if (group[i]!=0)
  {
    SigmaCond=as.matrix(Sigma[[i]]);
    for (j in 1:ns[i])
    {
      SigmaCond[,j]=SigmaCond[,j]*(-zo[[i]][j]);
      SigmaCond[j,]=SigmaCond[j,]*(-zo[[i]][j]);
      zij=zh[[i,zo3Converter[[i]][j]]]%*%as.matrix(x0[ndz[[zo3Converter[[i]][j]]]]);#%predicted values for zi
      zTilde[[i]][,j]=zij*z[[i]][,j];#%upper adjusted values
    }
FzTilde=pmnorm(as.matrix(zTilde[[i]]), varcov = SigmaCond);
lambdai=dF(zTilde[[i]],SigmaCond,0, denominator=FALSE)/FzTilde;#general lambda
lambdai_z=z[[i]]*lambdai;#lambda sign adjusted
lambdai2=d2F(zTilde[[i]],SigmaCond,0)/FzTilde;#great lambda
lambda_z[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = nsMax);#lambda*zi
lambda[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = nsMax);#lambda
lambda2[[i]]=array(0,dim=c(length(lambdai[,1]),nsMax, nsMax));#Big lambda
lambda_zTilde[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = nsMax);#lambda1*zTilde
LAMBDA[[i]]=array(0,dim=c(length(lambdai[,1]),nsMax, nsMax));#Big lambda
zG[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = nsMax);
zTildeG[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = nsMax);
#Adding zeros instead of ommited equations in lambda
counterz=0;
  for (j in 1:nsMax)
  {
    if (zo3[i,j]!=0)
    {
      counterz=counterz+1;
      counterz1=0;
      lambda_z[[i]][,j]=lambdai_z[,counterz];
      lambda[[i]][,j]=lambdai[,counterz];
      lambda_zTilde[[i]][,j]=lambdai[,counterz]*zTilde[[i]][,counterz];
      zG[[i]][,j]=z[[i]][,counterz];
      zTildeG[[i]][,j]=zTilde[[i]][,counterz];
      #Hessian of lambda
      for (k in 1:nsMax)
      {
        if (zo3[i,k]!=0)
        {
          counterz1=counterz1+1;
          LAMBDA[[i]][,j,k]=lambdai2[,counterz,counterz1];
          lambda2[[i]][,j,k]=lambdai2[,counterz,counterz1]*SigmaCond[counterz,counterz1];
        }
      }
    }
  }
 }
}
#Combining lambda by group into lambda by outcome
lambda_zG=lambda_z; lambda_z=matrix(list(), noutcome, 1);
lambdaG=lambda; lambda=matrix(list(), noutcome, 1);
lambdaG2=lambda2; lambda2=matrix(list(), noutcome, 1);
lambda_zTildeG=lambda_zTilde; lambda_zTilde=matrix(list(), noutcome, 1);
LAMBDAG=LAMBDA; LAMBDA=matrix(list(), noutcome, 1);
zTildeOutcome=matrix(list(), noutcome, 1);
zOutcome=matrix(list(), noutcome, 1);
for (i in 1:ngroup)
{
  if (group[i]!=0)
  {
    lambda_z[[group[i]]]=rbind(lambda_z[[group[i]]],lambda_zG[[i]]);
    lambda[[group[i]]]=rbind(lambda[[group[i]]],lambdaG[[i]]);
    lambda2[[group[i]]]=abind(lambda2[[group[i]]],lambdaG2[[i]],along=1);
    lambda_zTilde[[group[i]]]=rbind(lambda_zTilde[[group[i]]],lambda_zTildeG[[i]]);
    LAMBDA[[group[i]]]=abind(LAMBDA[[group[i]]],LAMBDAG[[i]],along=1);
    zOutcome[[group[i]]]=rbind(zOutcome[[group[i]]],zG[[i]]);
    zTildeOutcome[[group[i]]]=rbind(zTildeOutcome[[group[i]]],zTildeG[[i]]);
  }
}
#Combining y and yh by outcome
y1=matrix(list(), noutcome, 1);
yh1=matrix(list(), noutcome, 1);
for (i in 1:ngroup)
{
  if (group[i]!=0)
  {
    y1[[group[i]]]=rbind(y1[[group[i]]],y[[i]]);
    yh1[[group[i]]]=rbind(yh1[[group[i]]],yh[[i]]);
  }
}
twostep=matrix(list(), noutcome, 1);
twostepOLS=matrix(list(), noutcome, 1);
nCoef=matrix(0,1,noutcome);
X=matrix(list(),noutcome,1);
CovB=matrix(list(), noutcome, 1);
errors=matrix(list(),noutcome,1);
rhoSigma=matrix(list(),noutcome,1);
rhoY=matrix(list(),noutcome,1);
sigma=matrix(0,noutcome,1)
Ggamma=matrix(list(),1,noutcome);
yVariance=matrix(list(),1,noutcome);#variance of y
for (i in 1:noutcome)
{
X[[i]]=data.frame(yh1[[i]],lambda_z[[i]])
model=lm(y1[[i]]~.,data=X[[i]][,-1]);
twostep[[i]]=model;
nCoef[i]=length(x0[coef[[i]]])
x0[coef[[i]]]=coef(model)[1:nCoef[i]];#coefficients
coefLambda=coef(model)[(nCoef[i]+1):length(coef(model))];#coefficients
coefLambda[is.na(coefLambda)]=0;
parameters[coef[[i]]]=variable.names(model)[1:nCoef[i]];#Store coefficients names
Ggamma[[i]]=matrix(0,length(y1[[i]]),sum(nzh));
#sigma and ? calculation
startCoef1=1;
yVariance[[i]]=0;
yVariance[[i]]=lambda_zTilde[[i]]%*%coefLambda^2+
                            (lambda_z[[i]]%*%coefLambda)^2;
for (t in (1:ngroup)[group==i])#For each group corresponding to this outcome
{
  startCoef=1;
  endCoef1=startCoef1+groupsize[t]-1#get rows of this group
  counterz=0;
  for (k in 1:nsMax)#for each selection equation
  {
    l=0#part common for all coefficients of this selection equation
    endCoef=startCoef+nzh[k]-1;#number of coefficients
    if (zo3[t,k]!=0)
    {
      counterz=counterz+1;
      counterz1=0;
      for (j in 1:nsMax)#for each other equation
      {
        if (zo3[t,j]!=0)
        {
          counterz1=counterz1+1
          if (k!=j)
          {
            yVariance[[i]][startCoef1:endCoef1]=yVariance[[i]][startCoef1:endCoef1]-coefLambda[k]*z[[t]][,counterz]*z[[t]][,counterz1]*(coefLambda[j]-Sigma0[k,j]*coefLambda[k])*LAMBDAG[[t]][,k,j];
            l=l+coefLambda[j]*z[[t]][,counterz]*z[[t]][,counterz1]*(LAMBDAG[[t]][,k,j]-lambdaG[[t]][,k]*lambdaG[[t]][,j])-coefLambda[k]*LAMBDAG[[t]][,k,j]*z[[t]][,counterz]*z[[t]][,counterz1]*Sigma0[k,j];
          }
        }
      }
      Ggamma[[i]][startCoef1:endCoef1,startCoef:endCoef]=zh[[t,k]]*(l-coefLambda[k]*(lambda_zTildeG[[t]][,k]+lambdaG[[t]][,k]^2));
    }
    startCoef=endCoef+1;
  }
  startCoef1=endCoef1+1;
}
errors[[i]]=predict(model)-y1[[i]];
x0[nSigma-noutcome+i]=sqrt((sum(errors[[i]]^2)+sum(yVariance[[i]]))/length(y1[[i]]));
Ggamma[[i]]=Ggamma[[i]]/x0[nSigma-noutcome+i];#adjust for beta instead or rho
yVariance[[i]]=x0[nSigma-noutcome+i]^2-yVariance[[i]]
#Dealing with rho
rhoSigma[[i]]=coef(model)[(nCoef[i]+1):(nCoef[i]+nsMax)];
rhoSigma[[i]][is.na(rhoSigma[[i]])]=0;
sigma[i]=x0[nSigma-noutcome+i];
rhoY[[i]]=rhoSigma[[i]]/sigma[i];
}
#Delta method
  #G and triangle matrices
triangle=matrix(list(),noutcome,1);
for (i in 1:noutcome)
{
  triangle[[i]]=diag(1-yVariance[[i]]/sigma[i]^2);
}
  #Covariance matrix
for (i in 1:noutcome)
{
Q=as.matrix(X[[i]]);
Id=diag(rep(1,length(yh1[[i]][,1])));
CovB[[i]]=sigma[i]^2*(solve(t(Q)%*%Q)%*%t(Q)%*%((Id-triangle[[i]])+Ggamma[[i]]%*%covmatrixGamma%*%t(Ggamma[[i]]))%*%Q%*%solve(t(Q)%*%Q));
  #Change covmatrix for summary
newCoefs=coeftest(twostep[[i]], vcov. = CovB[[i]])
twostepOLS[[i]]=twostep[[i]];
twostep[[i]]=summary(twostep[[i]]);
twostep[[i]]$coefficients=newCoefs;
twostep[[i]]$sigma=sigma[i];
}
if (noutcome==1)
{
  twostep=twostep[[1]];
}
if (onlyTwostep) {return(list("model"=twostep, "covmatrix"=CovB,"sigma"=sigma, "twostepLS"=twostepOLS, "sortList"=sortList));}
#PAHSE 3: MLE with Two-step initial values
#Set lower and upper bound constraints for parameters
rhoSize=nSigma-noutcome;
sizeAboveRho=length(x0)-rhoSize;
lb=(c(rep(-1,rhoSize),rep(-Inf,sizeAboveRho)));
ub=(c(rep(1,rhoSize),rep(Inf,sizeAboveRho)));
x0[is.na(x0)]=0;
#Estimate coefficients and store them to MLE
f<-nloptr(x0=x0, eval_f=gheckmanLikelihood,opts=opts, lb=lb, ub=ub, y=y, zh=zh, yh=yh, zo=zo, ns=ns, ndz=ndz, nSigma=nSigma, coef=coef, group=group, ngroup=ngroup, nsMax=nsMax, zo3Converter=zo3Converter, noutcome=noutcome, zo3=zo3, groupsize=groupsize, nrhoY=nrhoY, ShowInfo=ShowInfo, maximization=FALSE);
stdev=jacobian(func = gheckmanGradient, x = f$solution, y=y, zh=zh, yh=yh, zo=zo, ns=ns, ndz=ndz, nSigma=nSigma, coef=coef, group=group, ngroup=ngroup, nsMax=nsMax, zo3Converter=zo3Converter, noutcome=noutcome, zo3=zo3, groupsize=groupsize, nrhoY=nrhoY, ShowInfo=FALSE);
CovM=solve(stdev)
stdev=sqrt(diag(CovM));
solution=f$solution;
x0=f$solution;
aSTD=pnorm(solution/stdev);
bSTD=pnorm(-solution/stdev);
pvalue=x0*0;
for (i in (1:length(x0)))
{
  pvalue[i]=1-max(aSTD[i],bSTD[i])+min(aSTD[i],bSTD[i]);
}
#Orginizing output
counter=0;
for (i in 1:nsMax)
{
  for (j in (i+1):nsMax)
  {
    if (j>i & j<=nsMax)
    {
    counter=counter+1;
    parameters[counter]=paste(c("rhoZ",i,j), collapse = " ");
    }
  }
}
for (i in 1:noutcome)
{
  for (j in 1:nsMax)
  {
    counter=counter+1;
    parameters[counter]=paste(c("rhoY",i,j), collapse = " ");
  }
}
for (i in 1:noutcome)
{
  counter=counter+1;
  parameters[counter]=paste(c("sigma",i), collapse = " ");
}
#Calculating lambda for MLE
#Calculating lambda
zTilde=matrix(list(), ngroup, 1);
lambda_z=matrix(list(), ngroup, 1);
lambda=matrix(list(), ngroup, 1);
lambda2=matrix(list(), ngroup, 1);
lambda_zTilde=matrix(list(), ngroup, 1);
LAMBDA=matrix(list(), ngroup, 1);
zG=matrix(list(), ngroup, 1);
zTildeG=matrix(list(), ngroup, 1);
Sigma=matrix(list(), ngroup, 1);
Sigma0=triangular(x0[1:nrhoZ],rep(1,nsMax));
for (i in 1:ngroup)
{
  Sigma[[i]]=Sigma0[zo3[i,]!=0,zo3[i,]!=0];
  zTilde[[i]]=z[[i]]*NA;
  if (group[i]!=0)
  {
    SigmaCond=as.matrix(Sigma[[i]]);
    for (j in 1:ns[i])
    {
      SigmaCond[,j]=SigmaCond[,j]*(-zo[[i]][j]);
      SigmaCond[j,]=SigmaCond[j,]*(-zo[[i]][j]);
      zij=zh[[i,zo3Converter[[i]][j]]]%*%as.matrix(x0[ndz[[zo3Converter[[i]][j]]]]);#%predicted values for zi
      zTilde[[i]][,j]=zij*z[[i]][,j];#%upper adjusted values
    }
    FzTilde=pmnorm(as.matrix(zTilde[[i]]), varcov = SigmaCond);
    lambdai=dF(zTilde[[i]],SigmaCond,0, denominator=FALSE)/FzTilde;#general lambda
    lambdai_z=z[[i]]*lambdai;#lambda sign adjusted
    lambdai2=d2F(zTilde[[i]],SigmaCond,0)/FzTilde;#great lambda
    lambda_z[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = nsMax);#lambda*zi
    lambda[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = nsMax);#lambda
    lambda2[[i]]=array(0,dim=c(length(lambdai[,1]),nsMax, nsMax));#Big lambda
    lambda_zTilde[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = nsMax);#lambda1*zTilde
    LAMBDA[[i]]=array(0,dim=c(length(lambdai[,1]),nsMax, nsMax));#Big lambda
    zG[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = nsMax);
    zTildeG[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = nsMax);
    #Adding zeros instead of ommited equations in lambda
    counterz=0;
    for (j in 1:nsMax)
    {
      if (zo3[i,j]!=0)
      {
        counterz=counterz+1;
        counterz1=0;
        lambda_z[[i]][,j]=lambdai_z[,counterz];
        lambda[[i]][,j]=lambdai[,counterz];
        lambda_zTilde[[i]][,j]=lambdai[,counterz]*zTilde[[i]][,counterz];
        zG[[i]][,j]=z[[i]][,counterz];
        zTildeG[[i]][,j]=zTilde[[i]][,counterz];
        #Hessian of lambda
        for (k in 1:nsMax)
        {
          if (zo3[i,k]!=0)
          {
            counterz1=counterz1+1;
            LAMBDA[[i]][,j,k]=lambdai2[,counterz,counterz1];
            lambda2[[i]][,j,k]=lambdai2[,counterz,counterz1]*SigmaCond[counterz,counterz1];
          }
        }
      }
    }
  }
}
#Combining lambda by group into lambda by outcome
lambda_zG=lambda_z; lambda_z=matrix(list(), noutcome, 1);
lambdaG=lambda; lambda=matrix(list(), noutcome, 1);
lambdaG2=lambda2; lambda2=matrix(list(), noutcome, 1);
lambda_zTildeG=lambda_zTilde; lambda_zTilde=matrix(list(), noutcome, 1);
LAMBDAG=LAMBDA; LAMBDA=matrix(list(), noutcome, 1);
zTildeOutcome=matrix(list(), noutcome, 1);
zOutcome=matrix(list(), noutcome, 1);
for (i in 1:ngroup)
{
  if (group[i]!=0)
  {
    lambda_z[[group[i]]]=rbind(lambda_z[[group[i]]],lambda_zG[[i]]);
    lambda[[group[i]]]=rbind(lambda[[group[i]]],lambdaG[[i]]);
    lambda2[[group[i]]]=abind(lambda2[[group[i]]],lambdaG2[[i]],along=1);
    lambda_zTilde[[group[i]]]=rbind(lambda_zTilde[[group[i]]],lambda_zTildeG[[i]]);
    LAMBDA[[group[i]]]=abind(LAMBDA[[group[i]]],LAMBDAG[[i]],along=1);
    zOutcome[[group[i]]]=rbind(zOutcome[[group[i]]],zG[[i]]);
    zTildeOutcome[[group[i]]]=rbind(zTildeOutcome[[group[i]]],zTildeG[[i]]);
  }
}
result=noquote(cbind(parameters,f$solution,stdev,pvalue));
colnames(result)=c("Parameter","value","stdev","p-value");
return(list("mle"=list("result" = result, "coefficients"=f$solution,"stdev"=stdev,"p-value"=pvalue,"names"=parameters), "twostep"=list("model"=twostep,"covmatrix"=CovB,"sigma"=sigma, "twostepLS"=twostepOLS), "logLikelihood"=-f$objective, "x0"=f$solution,"lambda"=lambdaG, "sortList"=sortList, "CovM"=CovM))
}

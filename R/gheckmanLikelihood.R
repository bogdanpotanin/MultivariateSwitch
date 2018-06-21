gheckmanLikelihood<-function(x, y, zh, yh, zo, ns, ndz, nSigma, coef, group, ngroup, nsMax, zo3Converter, noutcome, zo3, groupsize, nrhoY, ShowInfo=FALSE, maximization=FALSE)
{
  maxNumber=Inf#sqrt(.Machine$double.xma);
  #PHASE 1: Preparing data
  #Coefficients variables initialization
  f=0;
  beta=matrix(list(), noutcome, 1);#coefficients for y
  for (i in 1:noutcome)
  {
    beta[[i]]=as.matrix(x[coef[[i]]]);#where beta[[i]] is a column of i-th group coefficients
    #colnames(beta[[i]])=colnames(yh[[i]])
  }
  alpha=matrix(list(), nsMax, 1)#coefficients for each z[i]
  for (i in 1:nsMax )
  {
    alpha[[i]]=t(x[ndz[[i]]]);#where alpha(:,i) is column of zi coefficients
    #colnames(alpha[[i]])=colnames(zh[[i]])
  }
  #Covariation matrix
  #distinguishing different types of elements
  rhoZ=x[1:sum(1:(nsMax-1))];#correlations between different z disturbances with each other
  rhoY=matrix(list(), nrow=noutcome, ncol=1);#correlations between different z disturbances and y disturbance
  sigma=matrix(nrow=noutcome,ncol=1);#disturbances for different outcomes
  for (i in 1:noutcome)
  {
    rhoY[[i]]=as.matrix(x[nrhoY[[i]]]);#correlations between z and y
    sigma[i]=x[nSigma-noutcome+i];#sigma represents variance of y disturbances
  }
  #if (sigma==0) {return(Inf);}#For GenSa algorithm
  #Matrix of z
  Sigma0=triangular(rhoZ,rep(1,nsMax));#triangular(right corner elements, diagonal emelents)
  Sigma0=rbind(cbind(Sigma0,rep(0,nsMax)),rep(0,nsMax+1));
  if (nsMax==1) {Sigma0=matrix(c(1,0,0,0),ncol=2);}
  Sigma=matrix(list(), noutcome, 1);
  #Adding y disturbances to the matrix
  if (max(group!=0))
  {
    for (i in 1:noutcome)
    {
      Sigma[[i]]=Sigma0;
      Sigma[[i]][1:(nsMax+1),nsMax+1]=t(c(rhoY[[i]],sigma[i]))*sigma[i];#sigma])*sigma in order to include sigma^2 in matrix
      Sigma[[i]][nsMax+1,1:nsMax]=rhoY[[i]]*sigma[i];
    }
  }
  #else {Sigma0=Sigma0[1:nsMax,1:nsMax]}#for GenSa algorithm
  #Assigning vectors for gradient
  g=x*0; #final gradient vector
  rhoZGradient=rhoZ*0;
  rhoYGradient=matrix(list(), noutcome, 1);
  for (i in 1:noutcome)
  {
    rhoYGradient[[i]]=matrix(0,nrow=nsMax,ncol=1);
  }
  betaGradient=matrix(list(), noutcome, 1);
  for (i in 1:noutcome)
  {
    betaGradient[[i]]=beta[[i]]*0;
  }
  alphaGradient=matrix(list(), nsMax, 1);
  for (i in 1:nsMax)
  {
    alphaGradient[[i]]=matrix(0,nrow=length(ndz[[i]]),ncol=1);
  }
  sigmaGradient=matrix(rep(0,noutcome),ncol=1);
  #Disturbances and selection probabilities estimation
  SigmaPositively=0;
  if(!is.null(Sigma[[1]]))
  {
    for (i in 1:noutcome)
    {
      SigmaPositively=SigmaPositively+is.positive.definite(as.matrix(Sigma[[i]]), tol=1e-16);
    }
  }
  else {if (is.positive.semi.definite(Sigma0)){SigmaPositively=noutcome;}}
  if (SigmaPositively==noutcome) #check wheather Sigma is positively defined
  {
    epsilon=matrix(list(), noutcome, 1);#disturbances of y
    #y disturbances estimation
    for (i in 1:(ngroup))
    {
      if (group[i]!=0)
      {
        epsilon[[i]]=y[[i]]-as.matrix(yh[[i]])%*%beta[[group[i]]];
      }
    }
  }
  else
  {
    #if (!grad) {return(Inf)}#For GenSa
    if (!maximization) {return(list(objective=maxNumber, gradient=rep(maxNumber,length(x))));}
    else {return(list(objective=-maxNumber, gradient=rep(-maxNumber,length(x))));}
  }
  #z selection probabilities
  z=matrix(list(), noutcome, 1);
  for (i in 1:ngroup)
  {
    z[[i]]=matrix(0,ncol=ns[i],nrow=groupsize[i])
    for (j in 1:ns[i])
    {
        z[[i]][,j]=zh[[i,zo3Converter[[i]][j]]]%*%t(alpha[[zo3Converter[[i]][j]]]);
    }
  }
  #PHASE 2: Likelihood estimation
  #Variables to store observations likelihoods
  observed1=matrix(list(), ngroup, 1);#probability of z conditional on y disturbances
  observed2=matrix(list(), ngroup, 1);#unconditional probability of y disturbances
  #Conditional covariance matrices and mean
  SigmaCond=matrix(list(), ngroup, 1);#array of conditional matrices for each group
  k=matrix(list(), ngroup, 1);#argument of multinominal normal probability function
  #print(Sigma[[1]])
  #Likelihood calculation for each group
  for (i in (1:(ngroup)))
  {
    if (group[i]!=0) #if we have observations for y in this group
    {
      #Calculating conditional mean and covariance matrix for i-th group
      listCond=mvncond(mean=matrix(0,nrow=groupsize[i],ncol=ns[i]+1), sigma=Sigma[[group[i]]][zo3[i,]!=0,zo3[i,]!=0], dependent.ind = 1:ns[i], given.ind = ns[i]+1, X.given=epsilon[[i]]);
      SigmaCond[[i]]=listCond[[1]];
      meanCond=listCond[[2]];
      #meanCond1-meanCond
      #Estimating unconditional probabilities of y disturbances
      observed2[[i]]=dnorm(epsilon[[i]], mean=0, sd=sigma[group[i]]);
    }
    else #so if we have no observations for y
    {
      SigmaCond[[i]]=as.matrix(Sigma0[c(zo3[i,]!=0,FALSE),c(zo3[i,]!=0,FALSE)]);#Choose part of Sigma where selection equation presists
      meanCond=0;
    }
    #Adjusting signs of covariance matrix
    for (j in 1:(ns[i]))
    {
      SigmaCond[[i]][,j]=SigmaCond[[i]][,j]*zo[[i]][j];
      SigmaCond[[i]][j,]=SigmaCond[[i]][j,]*zo[[i]][j];
    }
    #Calculating argument for multinominal normal probability
    k[[i]]<-sweep((z[[i]]+meanCond),MARGIN=2,as.vector(zo[[i]]),'*');
    if (ns[i]==2)
      {
        rho=SigmaCond[[i]][1,2]/sqrt(SigmaCond[[i]][1,1]*SigmaCond[[i]][2,2]);#correlation for standtadtised distribution
        if (abs(rho)<=0.99) {observed1[[i]]=pbivnorm(x = cbind(k[[i]][,1]/sqrt(SigmaCond[[i]][1,1]), k[[i]][,2]/sqrt(SigmaCond[[i]][2,2])), rho = rho);}
        else {return(list(objective=maxNumber, gradient=rep(maxNumber,length(x))));}#if standartisation close to singular
      }
    else {observed1[[i]]=pmnorm(as.matrix(k[[i]]), varcov = SigmaCond[[i]]);}
    #Summarizing the results
    if (group[i]!=0)
    {
      f=f-sum(log(observed1[[i]])+log(observed2[[i]]));
    }
    else
    {
      f=f-sum(log(observed1[[i]]));
    }
  }
  #if (!grad) {return(f);}#For GenSa
  #PHASE 3: Gradient calculation
  for (i in 1:ngroup)#%for each group
  {
    rhoZGradientCounter=1;
    #Assigning variables that change per group
    SigmaCond2=array(0,dim=c(ns[i]-1,ns[i]-1,ns[i]));#array recalculated for each group
    k2=array(0,dim=c(groupsize[i],ns[i]-1,ns[i]));#(observation, conditioned arguments)
    SigmaConds12s22=matrix(0,nrow=ns[i]-1,ncol=ns[i]);#parts by which mean multiplied under condition
    dFdk=matrix(0,nrow=groupsize[i],ncol=ns[i])
    FkCondition2=matrix(0,nrow=groupsize[i],ncol=ns[i]);
    fkCondition2=matrix(0,nrow=groupsize[i],ncol=ns[i]);
    for (j in 1:ns[i])#for each argument of mvncdf
    {
      #Lets denote dydx partial derivative of y respect to x
      #Seperate cycle in order to have all parameters while dealing
      #with second derevatives
      #Estimating conditional parameters
      if (ns[i]>1)#because else dF(x1|x2)/dx2 includes only pdf part
      {
        listCond=mvncond(mean=matrix(0,nrow=groupsize[i],ncol=ns[i]), sigma=SigmaCond[[i]], dependent.ind=c((1:ns[i])[-j]), given.ind = j, X.given=k[[i]][,j]);
        SigmaCond2[,,j]=listCond[[1]];
        meanCond2=listCond[[2]];
        SigmaConds12s22[,j]=listCond[[3]];
        #Calculating conditional probabilities
        k2[,,j]=k[[i]][,-j]-meanCond2;#argument with j-th element conditioned
        FkCondition2[,j]=pmnorm(as.matrix(k2[,,j]),varcov=SigmaCond2[,,j]);
      }
      else
      {
        FkCondition2[,j]=1;
      }
      fkCondition2[,j]=dnorm(k[[i]][,j], mean=0, sd=sqrt(SigmaCond[[i]][j,j]));
    }
      for (j in 1:ns[i])#for each argument of mvncdf
      {
        jConverted=zo3Converter[[i]][j];#the number of j in zo3
        #Partial derivatives estimation
        dFdk[,j]=FkCondition2[,j]*fkCondition2[,j]/observed1[[i]];
        dkda=zh[[i,jConverted]]*zo[[i]][j];
        #Adding derivative in this group for alpha gradient
        alphaGradient[[jConverted]]=alphaGradient[[jConverted]]-t(dkda)%*%dFdk[,j];
        #Calculating second derivatives for rhoZ and rhoY
      if (j<=(ns[i]-1))
      {
          for (l in j:(ns[i]-1))#for each unconditioned argument of mvncdf greater then j
        {
          #Calculating conditional probabilities
          j2=l+1;#which of ns arguments is conditioned second
          j2Converted=zo3Converter[[i]][j2];
          if (ns[i]>2)#because else dF(x1|x2)/dx2 includes only pdf part
          {
            listCond=mvncond(mean=matrix(0,nrow=groupsize[i],ncol=ns[i]-1), sigma=SigmaCond2[,,j], dependent.ind=c((1:(ns[i]-1))[-l]), given.ind=l, X.given=k2[,l,j]);
            SigmaCond3=listCond[[1]];
            meanCond3=listCond[[2]];
            k3=k2[,-l,j]-meanCond3;#argument with j-th and j2-th elements conditioned
            FkCondition3=pmnorm(x=as.matrix(k3), varcov = SigmaCond3);
            #FkCondition3=apply(X=k3,MARGIN=1, FUN=function(x) pmvnorm(mean = 0, sigma = SigmaCond3, lower=rep(-Inf,ns[i]-2), upper=x, algorithm = algorithm));
          }
          else
          {
            FkCondition3=1;
          }
          fkCondition3=dmvnorm(cbind(k[[i]][,j],k[[i]][,j2]),c(0,0),SigmaCond[[i]][c(j,j2),c(j,j2)])
          #Calculating derivative of F respect to j and j2
          dFdk2=t(FkCondition3*fkCondition3/observed1[[i]]);#transpose to simplify furher multiplication
          #Adding this derivative to the rhoZ gradient
          shareblePart=sum(dFdk2)*zo[[i]][j]*zo[[i]][j2];#part common both for rhoZ and convariance part of rhoY
          rhoZGradient[zo3Converter[[i]][rhoZGradientCounter]]=rhoZGradient[zo3Converter[[i]][rhoZGradientCounter]]-shareblePart;
          rhoZGradientCounter=rhoZGradientCounter+1;
          #Adding adjusted derivatives to the rhoY gradient (plus
          #because multiplied by -rhoY and -s12s22)
          if (group[i]!=0)
          {
            #covariance part
            rhoYGradient[[group[i]]][jConverted]=rhoYGradient[[group[i]]][jConverted]+shareblePart*rhoY[[group[i]]][j2Converted];
            rhoYGradient[[group[i]]][j2Converted]=rhoYGradient[[group[i]]][j2Converted]+shareblePart*rhoY[[group[i]]][jConverted];
            #variance second part
            rhoYGradient[[group[i]]][jConverted]=rhoYGradient[[group[i]]][jConverted]-sum(dFdk2*SigmaConds12s22[l,j])*(rhoY[[group[i]]][jConverted]);#do not multiplied by 2 because devided by 2 due to sigma derivative formula
            rhoYGradient[[group[i]]][j2Converted]=rhoYGradient[[group[i]]][j2Converted]-sum(dFdk2*SigmaConds12s22[j,j2])*(rhoY[[group[i]]][j2Converted]);
          }
        }
      }
        if (group[i]!=0)
        {
          #Gradient for beta: mvncdf part
          dkdb=-yh[[i]]*rhoY[[group[i]]][jConverted]*zo[[i]][j]/sigma[group[i],1];
          dFdb=t(dkdb)%*%dFdk[,j];
          betaGradient[[group[i]]]=betaGradient[[group[i]]]-dFdb;
          #Gradient for rhoY variance first part
          dfdkj=-fkCondition2[,j]*k[[i]][,j]/SigmaCond[[i]][j,j];
          rhoYGradient[[group[i]]][j]=rhoYGradient[[group[i]]][j]+sum(FkCondition2[,j]*dfdkj/observed1[[i]])*(rhoY[[group[i]]][jConverted]);
          #Gradient for rhoY mvncdf argument part
          dFdrhoYj=t(dFdk[,j])%*%epsilon[[i]]*(zo[[i]][j]/sigma[group[i],1]);
          rhoYGradient[[group[i]]][j]=rhoYGradient[[group[i]]][j]-sum(dFdrhoYj);
          #Gradient for sigma mvncdf argument part
          sigmaGradient[group[i],1]=sigmaGradient[group[i],1]+(dFdrhoYj*rhoY[[group[i]]][j])/sigma[group[i],1];
        }
      }
      #Gradient for beta: mvnpdf part
      if (group[i]!=0)
      {
        dfdb=t(yh[[i]])%*%epsilon[[i]]/(sigma[group[i],1]^2);
        betaGradient[[group[i]]]=betaGradient[[group[i]]]-dfdb;
        #Gradient for sigma mvnpdf part
        sigmaGradient[group[i],1]=sigmaGradient[group[i],1]-sum((epsilon[[i]]^2-sigma[group[i],1]^2)/(sigma[group[i],1]^3));
      }
  }
    #Gradient vector construction
    g=rhoZGradient;
    for (i in 1:noutcome)
    {
      g=c(g,rhoYGradient[[i]]);
    }
    for (i in 1:noutcome)
    {
      g=c(g,t(sigmaGradient[i]));
    }
    for (i in 1:noutcome)
    {
      g=c(g,t(betaGradient[[i]]));
    }
    for (i in 1:nsMax)
    {
      g=c(g,t(alphaGradient[[i]]));
    }
    #Show info
    if(ShowInfo)
    {
      print('Likelihood')
      print(f)
      if (max(group!=0))
      {
      for (i in 1:noutcome)
      {
      print(cat('Covariance matrix ',toString(i)))
      print(Sigma[[i]])
      }
      }
      else
      {
        print(Sigma0)
      }
      if (max(group!=0))
      {
      for (i in 1:max(group))
      {
      print(cat('y',toString(i)))
      print(t(beta[[i]]))
      }
      }
      for (i in 1:nsMax)
      {
      print(cat('z',toString(i)))
      print((alpha[[i]]))
      }
    }
    if (nsMax==1) {g=g[-1];}
    if(is.infinite(f)) {f=maxNumber/10;}
    g[is.infinite(g) | is.na(g)]=maxNumber/10;
    if (!maximization) {return(list(objective=f, gradient=t(g)));}
    else {return(list(objective=-f, gradient=-t(g)));}
}
gheckmanGradient<-function(x, y, zh, yh, zo, ns, ndz, nSigma, coef, group, ngroup, nsMax, zo3Converter, noutcome, zo3, groupsize, nrhoY, ShowInfo=FALSE, maximization=FALSE)
{
  return(gheckmanLikelihood(x, y, zh, yh, zo, ns, ndz, nSigma, coef, group, ngroup, nsMax, zo3Converter, noutcome, zo3, groupsize, nrhoY, ShowInfo=FALSE, maximization=FALSE)[[2]])
}

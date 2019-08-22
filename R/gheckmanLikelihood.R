gheckmanLikelihood<-function(x, y, z_variables, y_variables, rules_no_ommit, n_selection_equations, coef_z, sigma_last_index, coef_y, groups, n_groups, n_selection_equations_max, rules_converter, n_outcome, rules, groups_observations, rho_y_indices, show_info=FALSE, maximization=FALSE)
{
  maxNumber=Inf#sqrt(.Machine$double.xma);
  #PHASE 1: Preparing data
  #Coefficients variables initialization
  f=0;
  beta=matrix(list(), n_outcome, 1);#coefficients for y
  for (i in 1:n_outcome)
  {
    beta[[i]]=as.matrix(x[coef_y[[i]]]);#where beta[[i]] is a column of i-th groups coefficients
    #colnames(beta[[i]])=colnames(y_variables[[i]])
  }
  alpha=matrix(list(), n_selection_equations_max, 1)#coefficients for each z[i]
  for (i in 1:n_selection_equations_max )
  {
    alpha[[i]]=t(x[coef_z[[i]]]);#where alpha(:,i) is column of zi coefficients
    #colnames(alpha[[i]])=colnames(z_variables[[i]])
  }
  #Covariation matrix
  #distinguishing different types of elements
  rhoZ=x[1:sum(1:(n_selection_equations_max-1))];#correlations between different z disturbances with each other
  rhoY=matrix(list(), nrow=n_outcome, ncol=1);#correlations between different z disturbances and y disturbance
  sigma=matrix(nrow=n_outcome,ncol=1);#disturbances for different outcomes
  for (i in 1:n_outcome)
  {
    rhoY[[i]]=as.matrix(x[rho_y_indices[[i]]]);#correlations between z and y disturbances
    sigma[i]=x[sigma_last_index-n_outcome+i];#sigma represents variance of y disturbances
  }
  #if (sigma==0) {return(Inf);}#For GenSa algorithm
  #Matrix of z
  Sigma0=triangular(rhoZ,rep(1,n_selection_equations_max));#triangular(right corner elements, diagonal emelents)
  Sigma0=rbind(cbind(Sigma0,rep(0,n_selection_equations_max)),rep(0,n_selection_equations_max+1));
  if (n_selection_equations_max==1) {Sigma0=matrix(c(1,0,0,0),ncol=2);}
  Sigma=matrix(list(), n_outcome, 1);
  #Adding y disturbances to the matrix
  if (max(groups!=0))
  {
    for (i in 1:n_outcome)
    {
      Sigma[[i]]=Sigma0;
      Sigma[[i]][1:(n_selection_equations_max+1),n_selection_equations_max+1]=t(c(rhoY[[i]],sigma[i]))*sigma[i];#sigma])*sigma in order to include sigma^2 in matrix
      Sigma[[i]][n_selection_equations_max+1,1:n_selection_equations_max]=rhoY[[i]]*sigma[i];
    }
  }
  #else {Sigma0=Sigma0[1:n_selection_equations_max,1:n_selection_equations_max]}#for GenSa algorithm
  #Assigning vectors for gradient
  g=x*0; #final gradient vector
  rhoZGradient=rhoZ*0;
  rhoYGradient=matrix(list(), n_outcome, 1);
  for (i in 1:n_outcome)
  {
    rhoYGradient[[i]]=matrix(0,nrow=n_selection_equations_max,ncol=1);
  }
  betaGradient=matrix(list(), n_outcome, 1);
  for (i in 1:n_outcome)
  {
    betaGradient[[i]]=beta[[i]]*0;
  }
  alphaGradient=matrix(list(), n_selection_equations_max, 1);
  for (i in 1:n_selection_equations_max)
  {
    alphaGradient[[i]]=matrix(0,nrow=length(coef_z[[i]]),ncol=1);
  }
  sigmaGradient=matrix(rep(0,n_outcome),ncol=1);
  #Disturbances and selection probabilities estimation
  SigmaPositively=0;
  if(!is.null(Sigma[[1]]))
  {
    for (i in 1:n_outcome)
    {
      tryCatch({
      SigmaPositively=SigmaPositively+is.positive.definite(as.matrix(Sigma[[i]]), tol=1e-16);
      }, error = function(cond) {SigmaPositively = 0;} )
    }
  }
  else {if (is.positive.semi.definite(Sigma0)){SigmaPositively=n_outcome;}}
  if (SigmaPositively==n_outcome) #check wheather Sigma is positively defined
  {
    epsilon=matrix(list(), n_outcome, 1);#disturbances of y
    #y disturbances estimation
    for (i in 1:(n_groups))
    {
      if (groups[i]!=0)
      {
        epsilon[[i]]=y[[i]]-as.matrix(y_variables[[i]])%*%beta[[groups[i]]];
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
  z=matrix(list(), n_outcome, 1);
  for (i in 1:n_groups)
  {
    z[[i]]=matrix(0,ncol=n_selection_equations[i],nrow=groups_observations[i])
    for (j in 1:n_selection_equations[i])
    {
        z[[i]][,j]=z_variables[[i,rules_converter[[i]][j]]]%*%t(alpha[[rules_converter[[i]][j]]]);
    }
  }
  #PHASE 2: Likelihood estimation
  #Variables to store observations likelihoods
  observed1=matrix(list(), n_groups, 1);#probability of z conditional on y disturbances
  observed2=matrix(list(), n_groups, 1);#unconditional probability of y disturbances
  #Conditional covariance matrices and mean
  SigmaCond=matrix(list(), n_groups, 1);#array of conditional matrices for each groups
  k=matrix(list(), n_groups, 1);#argument of multinominal normal probability function
  #Likelihood calculation for each groups
  for (i in (1:(n_groups)))
  {
    if (groups[i]!=0) #if we have observations for y in this groups
    {
      #Calculating conditional mean and covariance matrix for i-th groups
      listCond=mvncond(mean=matrix(0,nrow=groups_observations[i],ncol=n_selection_equations[i]+1),
                       sigma=Sigma[[groups[i]]][c(rules[i,]!=0,TRUE),c(rules[i,]!=0,TRUE)],
                       dependent.ind = 1:n_selection_equations[i], given.ind = n_selection_equations[i]+1,
                       X.given=epsilon[[i]]);
      SigmaCond[[i]]=listCond[[1]];
      meanCond=listCond[[2]];
      #meanCond1-meanCond
      #Estimating unconditional probabilities of y disturbances
      observed2[[i]]=dnorm(epsilon[[i]], mean=0, sd=sigma[groups[i]]);
    }
    else #so if we have no observations for y
    {
      SigmaCond[[i]]=as.matrix(Sigma0[c(rules[i,]!=0,FALSE),c(rules[i,]!=0,FALSE)]);#Choose part of Sigma where selection equation presists
      meanCond=0;
    }
    #Adjusting signs of covariance matrix
    for (j in 1:(n_selection_equations[i]))
    {
      SigmaCond[[i]][,j]=SigmaCond[[i]][,j]*rules_no_ommit[[i]][j];
      SigmaCond[[i]][j,]=SigmaCond[[i]][j,]*rules_no_ommit[[i]][j];
    }
    #Calculating argument for multinominal normal probability
    k[[i]]<-sweep((z[[i]]+meanCond),MARGIN=2,as.vector(rules_no_ommit[[i]]),'*');
    if (n_selection_equations[i]==2)
      {
        tryCatch(
          {
        rho=SigmaCond[[i]][1,2]/sqrt(SigmaCond[[i]][1,1]*SigmaCond[[i]][2,2]);#correlation for standtadtised distribution
        if (abs(rho)<=0.99) {observed1[[i]]=pbivnorm(x = cbind(k[[i]][,1]/sqrt(SigmaCond[[i]][1,1]), k[[i]][,2]/sqrt(SigmaCond[[i]][2,2])), rho = rho);}
        else {return(list(objective=maxNumber, gradient=rep(maxNumber,length(x))));}#if standartisation close to singular
          },
          error=function(cond)
          {
            observed1[[i]]=pmnorm(as.matrix(k[[i]]), varcov = SigmaCond[[i]]);
            print(123)
          }
        )
      }
    else {observed1[[i]]=pmnorm(as.matrix(k[[i]]), varcov = SigmaCond[[i]]);}
    #Summarizing the results
    if (groups[i]!=0)
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
  for (i in 1:n_groups)#%for each groups
  {
    rhoZGradientCounter=1;
    #Assigning variables that change per groups
    SigmaCond2=array(0,dim=c(n_selection_equations[i]-1,n_selection_equations[i]-1,n_selection_equations[i]));#array recalculated for each groups
    k2=array(0,dim=c(groups_observations[i],n_selection_equations[i]-1,n_selection_equations[i]));#(observation, conditioned arguments)
    SigmaConds12s22=matrix(0,nrow=n_selection_equations[i]-1,ncol=n_selection_equations[i]);#parts by which mean multiplied under condition
    dFdk=matrix(0,nrow=groups_observations[i],ncol=n_selection_equations[i])
    FkCondition2=matrix(0,nrow=groups_observations[i],ncol=n_selection_equations[i]);
    fkCondition2=matrix(0,nrow=groups_observations[i],ncol=n_selection_equations[i]);
    for (j in 1:n_selection_equations[i])#for each argument of mvncdf
    {
      #Lets denote dydx partial derivative of y respect to x
      #Seperate cycle in order to have all parameters while dealing
      #with second derevatives
      #Estimating conditional parameters
      if (n_selection_equations[i]>1)#because else dF(x1|x2)/dx2 includes only pdf part
      {
        listCond=mvncond(mean=matrix(0,nrow=groups_observations[i],ncol=n_selection_equations[i]), sigma=SigmaCond[[i]], dependent.ind=c((1:n_selection_equations[i])[-j]), given.ind = j, X.given=k[[i]][,j]);
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
      for (j in 1:n_selection_equations[i])#for each argument of mvncdf
      {
        jConverted=rules_converter[[i]][j];#the number of j in rules
        #Partial derivatives estimation
        dFdk[,j]=FkCondition2[,j]*fkCondition2[,j]/observed1[[i]];
        dkda=z_variables[[i,jConverted]]*rules_no_ommit[[i]][j];
        #Adding derivative in this groups for alpha gradient
        alphaGradient[[jConverted]]=alphaGradient[[jConverted]]-t(dkda)%*%dFdk[,j];
        #Calculating second derivatives for rhoZ and rhoY
      if (j<=(n_selection_equations[i]-1))
      {
          for (l in j:(n_selection_equations[i]-1))#for each unconditioned argument of mvncdf greater then j
        {
          #Calculating conditional probabilities
          j2=l+1;#which of n_selection_equations arguments is conditioned second
          j2Converted=rules_converter[[i]][j2];
          if (n_selection_equations[i]>2)#because else dF(x1|x2)/dx2 includes only pdf part
          {
            listCond=mvncond(mean=matrix(0,nrow=groups_observations[i],ncol=n_selection_equations[i]-1), sigma=SigmaCond2[,,j], dependent.ind=c((1:(n_selection_equations[i]-1))[-l]), given.ind=l, X.given=k2[,l,j]);
            SigmaCond3=listCond[[1]];
            meanCond3=listCond[[2]];
            k3=k2[,-l,j]-meanCond3;#argument with j-th and j2-th elements conditioned
            FkCondition3=pmnorm(x=as.matrix(k3), varcov = SigmaCond3);
            #FkCondition3=apply(X=k3,MARGIN=1, FUN=function(x) pmvnorm(mean = 0, sigma = SigmaCond3, lower=rep(-Inf,n_selection_equations[i]-2), upper=x, algorithm = algorithm));
          }
          else
          {
            FkCondition3=1;
          }
          fkCondition3=dmvnorm(cbind(k[[i]][,j],k[[i]][,j2]),c(0,0),SigmaCond[[i]][c(j,j2),c(j,j2)])
          #Calculating derivative of F respect to j and j2
          dFdk2=t(FkCondition3*fkCondition3/observed1[[i]]);#transpose to simplify furher multiplication
          #Adding this derivative to the rhoZ gradient
          shareblePart=sum(dFdk2)*rules_no_ommit[[i]][j]*rules_no_ommit[[i]][j2];#part common both for rhoZ and convariance part of rhoY
          rhoZGradient[rules_converter[[i]][rhoZGradientCounter]]=rhoZGradient[rules_converter[[i]][rhoZGradientCounter]]-shareblePart;
          rhoZGradientCounter=rhoZGradientCounter+1;
          #Adding adjusted derivatives to the rhoY gradient (plus
          #because multiplied by -rhoY and -s12s22)
          if (groups[i]!=0)
          {
            #covariance part
            rhoYGradient[[groups[i]]][jConverted]=rhoYGradient[[groups[i]]][jConverted]+shareblePart*rhoY[[groups[i]]][j2Converted];
            rhoYGradient[[groups[i]]][j2Converted]=rhoYGradient[[groups[i]]][j2Converted]+shareblePart*rhoY[[groups[i]]][jConverted];
            #variance second part
            rhoYGradient[[groups[i]]][jConverted]=rhoYGradient[[groups[i]]][jConverted]-sum(dFdk2*SigmaConds12s22[l,j])*(rhoY[[groups[i]]][jConverted]);#do not multiplied by 2 because devided by 2 due to sigma derivative formula
            rhoYGradient[[groups[i]]][j2Converted]=rhoYGradient[[groups[i]]][j2Converted]-sum(dFdk2*SigmaConds12s22[j,j2])*(rhoY[[groups[i]]][j2Converted]);
          }
        }
      }
        if (groups[i]!=0)
        {
          #Gradient for beta: mvncdf part
          dkdb=-y_variables[[i]]*rhoY[[groups[i]]][jConverted]*rules_no_ommit[[i]][j]/sigma[groups[i],1];
          dFdb=t(dkdb)%*%dFdk[,j];
          betaGradient[[groups[i]]]=betaGradient[[groups[i]]]-dFdb;
          #Gradient for rhoY variance first part
          dfdkj=-fkCondition2[,j]*k[[i]][,j]/SigmaCond[[i]][j,j];
          rhoYGradient[[groups[i]]][j]=rhoYGradient[[groups[i]]][j]+sum(FkCondition2[,j]*dfdkj/observed1[[i]])*(rhoY[[groups[i]]][jConverted]);
          #Gradient for rhoY mvncdf argument part
          dFdrhoYj=t(dFdk[,j])%*%epsilon[[i]]*(rules_no_ommit[[i]][j]/sigma[groups[i],1]);
          rhoYGradient[[groups[i]]][j]=rhoYGradient[[groups[i]]][j]-sum(dFdrhoYj);
          #Gradient for sigma mvncdf argument part
          sigmaGradient[groups[i],1]=sigmaGradient[groups[i],1]+(dFdrhoYj*rhoY[[groups[i]]][j])/sigma[groups[i],1];
        }
      }
      #Gradient for beta: mvnpdf part
      if (groups[i]!=0)
      {
        dfdb=t(y_variables[[i]])%*%epsilon[[i]]/(sigma[groups[i],1]^2);
        betaGradient[[groups[i]]]=betaGradient[[groups[i]]]-dfdb;
        #Gradient for sigma mvnpdf part
        sigmaGradient[groups[i],1]=sigmaGradient[groups[i],1]-sum((epsilon[[i]]^2-sigma[groups[i],1]^2)/(sigma[groups[i],1]^3));
      }
  }
    #Gradient vector construction
    g=rhoZGradient;
    for (i in 1:n_outcome)
    {
      g=c(g,rhoYGradient[[i]]);
    }
    for (i in 1:n_outcome)
    {
      g=c(g,t(sigmaGradient[i]));
    }
    for (i in 1:n_outcome)
    {
      g=c(g,t(betaGradient[[i]]));
    }
    for (i in 1:n_selection_equations_max)
    {
      g=c(g,t(alphaGradient[[i]]));
    }
    #Show info
    if(show_info)
    {
      print('Likelihood')
      print(f)
      if (max(groups!=0))
      {
      for (i in 1:n_outcome)
      {
      print(cat('Covariance matrix ',toString(i)))
      print(Sigma[[i]])
      }
      }
      else
      {
        print(Sigma0)
      }
      if (max(groups!=0))
      {
      for (i in 1:max(groups))
      {
      print(cat('y',toString(i)))
      print(t(beta[[i]]))
      }
      }
      for (i in 1:n_selection_equations_max)
      {
      print(cat('z',toString(i)))
      print((alpha[[i]]))
      }
    }
    if (n_selection_equations_max==1) {g=g[-1];}
    if(is.infinite(f)) {f=maxNumber/10;}
    g[is.infinite(g) | is.na(g)]=maxNumber/10;
    if (!maximization) {return(list(objective=f, gradient=t(g)));}
    else {return(list(objective=-f, gradient=-t(g)));}
}
gheckmanGradient<-function(x, y, z_variables, y_variables, rules_no_ommit, n_selection_equations, coef_z, sigma_last_index, coef_y, groups, n_groups, n_selection_equations_max, rules_converter, n_outcome, rules, groups_observations, rho_y_indices, show_info=FALSE, maximization=FALSE)
{
  return(gheckmanLikelihood(x, y, z_variables, y_variables, rules_no_ommit, n_selection_equations, coef_z, sigma_last_index, coef_y, groups, n_groups, n_selection_equations_max, rules_converter, n_outcome, rules, groups_observations, rho_y_indices, show_info=FALSE, maximization=FALSE)[[2]])
}

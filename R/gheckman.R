#' Use this function to estimate Multivariate Switch model
#' @param data data frame containing the variables in the model
#' @param outcome main (continious) equation formula
#' @param selection1 the first selecton equation formula
#' @param selection2 the second selecton equation formula
#' @param selection3 the third selecton equation formula
#' @param selection4 the fourth selecton equation formula
#' @param selection5 the fifth selecton equation formula
#' @param groups vector which determines each groups outcome
#' @param rules matrix which rows correspond to possible combination of selection equations values.
#' @param ShowInfo shows likelihood function optimization info
#' @param only_twostep if true then only two-step procedure used
#' @param opts options that are passed to nlopt
#' @param x0 optional initial values to be used while solving optimization task
#' @param remove_zero_columns set true for switch regression (more then 1 groups) if your dependent variable is a part of dummy variable
#' @details This function estimates Multivariate Switch-Selection model via maximum-likelihood and two-step procedures.
#' This model was developed by E.V. Kossova and B.S. Potanin
#' Please cite as follows:
#' "Kossova, Elena & Potanin, Bogdan, 2018. 'Heckman method and switching regression model multivariate
#' generalization,' Applied Econometrics, Publishing House 'SINERGIA PRESS', vol. 50, pages 114-143."
#' Dependent variables in selection equations should have values -1 or 1.
#' Also there is a special value 0 which indicates that this observation is unobservable but it is necessary
#' to take into consideration other selection equation information while calculating likelihood function.
#' The i-th row of rules corresponds to i-th element of groups. rules rows should contain information regarding
#' possible combinations of selection equations values. While groups determines the outcome for each of this
#' possible combinations. Special 0 value for groups responsible for sample selection.
gheckman<-function(data, outcome, selection1=NULL, selection2=NULL, selection3=NULL,
                   selection4=NULL, selection5=NULL, groups=NULL, rules=NULL,
                   ShowInfo=FALSE, only_twostep=FALSE,
                   opts=list("algorithm" = "NLOPT_LD_TNEWTON", "xtol_rel" = 1e-16, "print_level" = 1, maxeval = 1000000),
                   x1=NULL, remove_zero_columns=FALSE)
{
  #PHASE 0: Extracting data from formulas
  y_variables=model.frame(formula = outcome, data=data, na.action = NULL);#data for main equation
  n_observations_total=length(y_variables[,1]);#total number of observations
  y_variables[rowSums(is.na(y_variables))>0,]=NA;#remove y that could not be calculated
  y=y_variables[,1];#independend variable
  y_variables[,1]=1;#intercept
  colnames(y_variables)[1]="intercept";
  #Data for selection equations
  z_variables=matrix(list(),5,1)#selection equation independent variables
  z=matrix(NA,n_observations_total,5);#selection equation dedepndent variables
  for (i in 1:5)#for each possible selection equation
  {
    selection=paste("selection",toString(i),sep="");#name of current selection equation
    if (!is.null(get(selection)))#if this selection equation assigned
    {
      z_variables[[i]]=model.frame(formula = as.formula(get(selection)), data=data, na.action = NULL)
      z_variables[[i]][rowSums(is.na(z_variables[[i]]))>0,]=NA;#remove z that could not be calculated
      z[,i]=z_variables[[i]][,1];
      z_variables[[i]][,1]=rep(1,n_observations_total);#intercept
      colnames(z_variables[[i]])[1]="intercept";
      z_variables[[i]][is.na(z_variables[[i]])]=0;
    }
    else#if it is no such selection equation
    {
      z_variables=z_variables[1:(i-1)];#preserve only existing selection equation
      z=z[,1:(i-1)];
      break
    }
  }
  #PHASE 1: Preparing data
  sortList<-gheckmanSort(y, y_variables, z, z_variables, groups, rules, remove_zero_columns);
  y=sortList$y;
  y_variables=sortList$y_variables;
  groups=sortList$groups;
  rules=sortList$rules;
  n_groups=sortList$n_groups;
  n_outcome=sortList$n_outcome;
  rules_converter=sortList$rules_converter;
  rules_no_ommit=sortList$rules_no_ommit;
  n_selection_equations=sortList$n_selection_equations;
  sigma_last_index=sortList$sigma_last_index;
  n_selection_equations_max=sortList$n_selection_equations_max;
  coef_y=sortList$coef_y;
  coef_z=sortList$coef_z;
  groupsize=sortList$groups_observations;
  nyh=sortList$n_y_variables;
  nzh=sortList$n_z_variables;
  rho_y_indices=sortList$rho_y_indices;
  rho_z_n=sum(1:(n_selection_equations_max-1));
  rhoSigma=matrix(list(),n_outcome);
  opts = opts;#setting optimization options max(maxeval/15,n_selection_equations_max*50)
  #PHASE 2: Two-step method
  numberOfParameters=coef_z[[n_selection_equations_max]][dim(coef_z[[n_selection_equations_max]])[1]]
  numberOfParametersY=coef_y[[n_outcome]][dim(coef_y[[n_outcome]])[1]]
  x0=matrix(0,numberOfParameters);
  parameters=vector(length = length(x0));#Store variables names
  #get initial values
  if (n_selection_equations_max==1) {z=matrix(z,ncol=1);}
  for (i in 1:n_selection_equations_max)
  {
    #-1 in order to exclude constant
    x0[coef_z[[i]]]=coef_y(myGlm<-glm(I((z[z[,i]!=0,i]+1)/2)~.-1,data=data.frame(z_variables[[i]][z[,i]!=0,]),family=binomial(link="probit")));
    parameters[coef_z[[i]]]=all.vars(formula(myGlm)[-2]);
  }
  zhColnames=matrix(list(),n_selection_equations_max);#zhColnames
  for (i in 1:n_selection_equations_max) {zhColnames[[i]]=colnames(z_variables[[i]]);}
  z_variables=sortList[[4]];#Only now we can load it from sort method result
  z=sortList[[3]];#Only now we can load it from sort method result
  #Set lower and upper bound constraints for parameters
  lb=c(rep(-0.9999999,rho_z_n),rep(0,numberOfParametersY-rho_z_n),rep(-Inf,numberOfParameters-coef_z[[1]][1]+1));
  ub=c(rep(0.9999999,rho_z_n),rep(0,numberOfParametersY-rho_z_n),rep(Inf,numberOfParameters-coef_z[[1]][1]+1));
  #Estimate coefficients and store them to x0
  f<-nloptr(x0=x0, eval_f=gheckmanLikelihood,opts=opts, lb=lb, ub=ub, y=y, z_variables=z_variables, y_variables=y_variables, rules_no_ommit=rules_no_ommit, n_selection_equations=n_selection_equations, coef_z=coef_z, sigma_last_index=sigma_last_index, coef_y=coef_y, groups=groups*0, n_groups=n_groups, n_selection_equations_max=n_selection_equations_max, rules_converter=rules_converter, n_outcome=n_outcome, rules=rules, groupsize=groupsize, rho_y_indices=rho_y_indices, ShowInfo=ShowInfo, maximization=FALSE);
  x0=f$solution;
  #x0=f$par;
  #Storing covariance matrix
  covmatrix=jacobian(func = gheckmanGradient, x = x0, y=y, z_variables=z_variables, y_variables=y_variables, rules_no_ommit=rules_no_ommit, n_selection_equations=n_selection_equations, coef_z=coef_z, sigma_last_index=sigma_last_index, coef_y=coef_y, groups=groups*0, n_groups=n_groups, n_selection_equations_max=n_selection_equations_max, rules_converter=rules_converter, n_outcome=n_outcome, rules=rules, groupsize=groupsize, rho_y_indices=rho_y_indices, ShowInfo=FALSE);
  covmatrixGamma=solve(covmatrix[coef_z[[1]][1]:length(x0),coef_z[[1]][1]:length(x0)]);
  #Calculating lambda
  zTilde=matrix(list(), n_groups, 1);
  lambda_z=matrix(list(), n_groups, 1);
  lambda=matrix(list(), n_groups, 1);
  lambda2=matrix(list(), n_groups, 1);
  lambda_zTilde=matrix(list(), n_groups, 1);
  LAMBDA=matrix(list(), n_groups, 1);
  zG=matrix(list(), n_groups, 1);
  zTildeG=matrix(list(), n_groups, 1);
  Sigma=matrix(list(), n_groups, 1);
  Sigma0=triangular(x0[1:rho_z_n],rep(1,n_selection_equations_max));
  for (i in 1:n_groups)
  {
    Sigma[[i]]=Sigma0[rules[i,]!=0,rules[i,]!=0];
    zTilde[[i]]=z[[i]]*NA;
    if (groups[i]!=0)
    {
      SigmaCond=as.matrix(Sigma[[i]]);
      for (j in 1:n_selection_equations[i])
      {
        SigmaCond[,j]=SigmaCond[,j]*(-rules_no_ommit[[i]][j]);
        SigmaCond[j,]=SigmaCond[j,]*(-rules_no_ommit[[i]][j]);
        zij=z_variables[[i,rules_converter[[i]][j]]]%*%as.matrix(x0[coef_z[[rules_converter[[i]][j]]]]);#%predicted values for zi
        zTilde[[i]][,j]=zij*z[[i]][,j];#%upper adjusted values
      }
      FzTilde=pmnorm(as.matrix(zTilde[[i]]), varcov = SigmaCond);
      lambdai=dF(zTilde[[i]],SigmaCond,0, denominator=FALSE)/FzTilde;#general lambda
      lambdai_z=z[[i]]*lambdai;#lambda sign adjusted
      lambdai2=d2F(zTilde[[i]],SigmaCond,0)/FzTilde;#great lambda
      lambda_z[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = n_selection_equations_max);#lambda*zi
      lambda[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = n_selection_equations_max);#lambda
      lambda2[[i]]=array(0,dim=c(length(lambdai[,1]),n_selection_equations_max, n_selection_equations_max));#Big lambda
      lambda_zTilde[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = n_selection_equations_max);#lambda1*zTilde
      LAMBDA[[i]]=array(0,dim=c(length(lambdai[,1]),n_selection_equations_max, n_selection_equations_max));#Big lambda
      zG[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = n_selection_equations_max);
      zTildeG[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = n_selection_equations_max);
      #Adding zeros instead of ommited equations in lambda
      counterz=0;
      for (j in 1:n_selection_equations_max)
      {
        if (rules[i,j]!=0)
        {
          counterz=counterz+1;
          counterz1=0;
          lambda_z[[i]][,j]=lambdai_z[,counterz];
          lambda[[i]][,j]=lambdai[,counterz];
          lambda_zTilde[[i]][,j]=lambdai[,counterz]*zTilde[[i]][,counterz];
          zG[[i]][,j]=z[[i]][,counterz];
          zTildeG[[i]][,j]=zTilde[[i]][,counterz];
          #Hessian of lambda
          for (k in 1:n_selection_equations_max)
          {
            if (rules[i,k]!=0)
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
  #Combining lambda by groups into lambda by outcome
  lambda_zG=lambda_z; lambda_z=matrix(list(), n_outcome, 1);
  lambdaG=lambda; lambda=matrix(list(), n_outcome, 1);
  lambdaG2=lambda2; lambda2=matrix(list(), n_outcome, 1);
  lambda_zTildeG=lambda_zTilde; lambda_zTilde=matrix(list(), n_outcome, 1);
  LAMBDAG=LAMBDA; LAMBDA=matrix(list(), n_outcome, 1);
  zTildeOutcome=matrix(list(), n_outcome, 1);
  zOutcome=matrix(list(), n_outcome, 1);
  for (i in 1:n_groups)
  {
    if (groups[i]!=0)
    {
      lambda_z[[groups[i]]]=rbind(lambda_z[[groups[i]]],lambda_zG[[i]]);
      lambda[[groups[i]]]=rbind(lambda[[groups[i]]],lambdaG[[i]]);
      lambda2[[groups[i]]]=abind(lambda2[[groups[i]]],lambdaG2[[i]],along=1);
      lambda_zTilde[[groups[i]]]=rbind(lambda_zTilde[[groups[i]]],lambda_zTildeG[[i]]);
      LAMBDA[[groups[i]]]=abind(LAMBDA[[groups[i]]],LAMBDAG[[i]],along=1);
      zOutcome[[groups[i]]]=rbind(zOutcome[[groups[i]]],zG[[i]]);
      zTildeOutcome[[groups[i]]]=rbind(zTildeOutcome[[groups[i]]],zTildeG[[i]]);
    }
  }
  #Combining y and y_variables by outcome
  y1=matrix(list(), n_outcome, 1);
  yh1=matrix(list(), n_outcome, 1);
  for (i in 1:n_groups)
  {
    if (groups[i]!=0)
    {
      y1[[groups[i]]]=rbind(y1[[groups[i]]],y[[i]]);
      yh1[[groups[i]]]=rbind(yh1[[groups[i]]],y_variables[[i]]);
    }
  }
  twostep=matrix(list(), n_outcome, 1);
  twostepOLS=matrix(list(), n_outcome, 1);
  nCoef=matrix(0,1,n_outcome);
  X=matrix(list(),n_outcome,1);
  CovB=matrix(list(), n_outcome, 1);
  errors=matrix(list(),n_outcome,1);
  rhoSigma=matrix(list(),n_outcome,1);
  rhoY=matrix(list(),n_outcome,1);
  sigma=matrix(0,n_outcome,1)
  Ggamma=matrix(list(),1,n_outcome);
  yVariance=matrix(list(),1,n_outcome);#variance of y
  for (i in 1:n_outcome)
  {
    X[[i]]=data.frame(yh1[[i]],lambda_z[[i]])
    model=lm(y1[[i]]~.,data=X[[i]][,-1]);
    twostep[[i]]=model;
    nCoef[i]=length(x0[coef_y[[i]]])
    x0[coef_y[[i]]]=coef_y(model)[1:nCoef[i]];#coefficients
    coefLambda=coef_y(model)[(nCoef[i]+1):length(coef_y(model))];#coefficients
    coefLambda[is.na(coefLambda)]=0;
    parameters[coef_y[[i]]]=variable.names(model)[1:nCoef[i]];#Store coefficients names
    Ggamma[[i]]=matrix(0,length(y1[[i]]),sum(nzh));
    #sigma and ? calculation
    startCoef1=1;
    yVariance[[i]]=0;
    yVariance[[i]]=lambda_zTilde[[i]]%*%coefLambda^2+
      (lambda_z[[i]]%*%coefLambda)^2;
    for (t in (1:n_groups)[groups==i])#For each groups corresponding to this outcome
    {
      startCoef=1;
      endCoef1=startCoef1+groupsize[t]-1#get rows of this groups
      counterz=0;
      for (k in 1:n_selection_equations_max)#for each selection equation
      {
        l=0#part common for all coefficients of this selection equation
        endCoef=startCoef+nzh[k]-1;#number of coefficients
        if (rules[t,k]!=0)
        {
          counterz=counterz+1;
          counterz1=0;
          for (j in 1:n_selection_equations_max)#for each other equation
          {
            if (rules[t,j]!=0)
            {
              counterz1=counterz1+1
              if (k!=j)
              {
                yVariance[[i]][startCoef1:endCoef1]=yVariance[[i]][startCoef1:endCoef1]-coefLambda[k]*z[[t]][,counterz]*z[[t]][,counterz1]*(coefLambda[j]-Sigma0[k,j]*coefLambda[k])*LAMBDAG[[t]][,k,j];
                l=l+coefLambda[j]*z[[t]][,counterz]*z[[t]][,counterz1]*(LAMBDAG[[t]][,k,j]-lambdaG[[t]][,k]*lambdaG[[t]][,j])-coefLambda[k]*LAMBDAG[[t]][,k,j]*z[[t]][,counterz]*z[[t]][,counterz1]*Sigma0[k,j];
              }
            }
          }
          Ggamma[[i]][startCoef1:endCoef1,startCoef:endCoef]=z_variables[[t,k]]*(l-coefLambda[k]*(lambda_zTildeG[[t]][,k]+lambdaG[[t]][,k]^2));
        }
        startCoef=endCoef+1;
      }
      startCoef1=endCoef1+1;
    }
    errors[[i]]=predict(model)-y1[[i]];
    x0[sigma_last_index-n_outcome+i]=sqrt((sum(errors[[i]]^2)+sum(yVariance[[i]]))/length(y1[[i]]));
    Ggamma[[i]]=Ggamma[[i]]/x0[sigma_last_index-n_outcome+i];#adjust for beta instead or rho
    yVariance[[i]]=x0[sigma_last_index-n_outcome+i]^2-yVariance[[i]]
    #Dealing with rho
    rhoSigma[[i]]=coef_y(model)[(nCoef[i]+1):(nCoef[i]+n_selection_equations_max)];
    rhoSigma[[i]][is.na(rhoSigma[[i]])]=0;
    sigma[i]=x0[sigma_last_index-n_outcome+i];
    rhoY[[i]]=rhoSigma[[i]]/sigma[i];
  }
  #Delta method
  #G and triangle matrices
  triangle=matrix(list(),n_outcome,1);
  for (i in 1:n_outcome)
  {
    triangle[[i]]=diag(1-yVariance[[i]]/sigma[i]^2);
  }
  #Covariance matrix
  for (i in 1:n_outcome)
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
  if (n_outcome==1)
  {
    twostep=twostep[[1]];
  }
  if (only_twostep) {return(list("model"=twostep, "covmatrix"=CovB,"sigma"=sigma, "twostepLS"=twostepOLS, "sortList"=sortList));}
  #PAHSE 3: MLE with Two-step initial values
  #Set lower and upper bound constraints for parameters
  rhoSize=sigma_last_index-n_outcome;
  sizeAboveRho=length(x0)-rhoSize;
  lb=(c(rep(-1,rhoSize),rep(-Inf,sizeAboveRho)));
  ub=(c(rep(1,rhoSize),rep(Inf,sizeAboveRho)));
  x0[is.na(x0)]=0;
  x0_twostep=x0
  #Estimate coefficients and store them to MLE
  if (!is.null(x1)) {x0<-x1;};#substitute some initial value
  f<-nloptr(x0=x0, eval_f=gheckmanLikelihood,opts=opts, lb=lb, ub=ub, y=y, z_variables=z_variables, y_variables=y_variables, rules_no_ommit=rules_no_ommit, n_selection_equations=n_selection_equations, coef_z=coef_z, sigma_last_index=sigma_last_index, coef_y=coef_y, groups=groups, n_groups=n_groups, n_selection_equations_max=n_selection_equations_max, rules_converter=rules_converter, n_outcome=n_outcome, rules=rules, groupsize=groupsize, rho_y_indices=rho_y_indices, ShowInfo=ShowInfo, maximization=FALSE);
  stdev=jacobian(func = gheckmanGradient, x = f$solution, y=y, z_variables=z_variables, y_variables=y_variables, rules_no_ommit=rules_no_ommit, n_selection_equations=n_selection_equations, coef_z=coef_z, sigma_last_index=sigma_last_index, coef_y=coef_y, groups=groups, n_groups=n_groups, n_selection_equations_max=n_selection_equations_max, rules_converter=rules_converter, n_outcome=n_outcome, rules=rules, groupsize=groupsize, rho_y_indices=rho_y_indices, ShowInfo=FALSE);
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
  for (i in 1:n_selection_equations_max)
  {
    for (j in (i+1):n_selection_equations_max)
    {
      if (j>i & j<=n_selection_equations_max)
      {
        counter=counter+1;
        parameters[counter]=paste(c("rhoZ",i,j), collapse = " ");
      }
    }
  }
  for (i in 1:n_outcome)
  {
    for (j in 1:n_selection_equations_max)
    {
      counter=counter+1;
      parameters[counter]=paste(c("rhoY",i,j), collapse = " ");
    }
  }
  for (i in 1:n_outcome)
  {
    counter=counter+1;
    parameters[counter]=paste(c("sigma",i), collapse = " ");
  }
  #Calculating lambda for MLE
  #Calculating lambda
  zTilde=matrix(list(), n_groups, 1);
  lambda_z=matrix(list(), n_groups, 1);
  lambda=matrix(list(), n_groups, 1);
  lambda2=matrix(list(), n_groups, 1);
  lambda_zTilde=matrix(list(), n_groups, 1);
  LAMBDA=matrix(list(), n_groups, 1);
  zG=matrix(list(), n_groups, 1);
  zTildeG=matrix(list(), n_groups, 1);
  Sigma=matrix(list(), n_groups, 1);
  Sigma0=triangular(x0[1:rho_z_n],rep(1,n_selection_equations_max));
  for (i in 1:n_groups)
  {
    Sigma[[i]]=Sigma0[rules[i,]!=0,rules[i,]!=0];
    zTilde[[i]]=z[[i]]*NA;
    if (groups[i]!=0)
    {
      SigmaCond=as.matrix(Sigma[[i]]);
      for (j in 1:n_selection_equations[i])
      {
        SigmaCond[,j]=SigmaCond[,j]*(-rules_no_ommit[[i]][j]);
        SigmaCond[j,]=SigmaCond[j,]*(-rules_no_ommit[[i]][j]);
        zij=z_variables[[i,rules_converter[[i]][j]]]%*%as.matrix(x0[coef_z[[rules_converter[[i]][j]]]]);#%predicted values for zi
        zTilde[[i]][,j]=zij*z[[i]][,j];#%upper adjusted values
      }
      FzTilde=pmnorm(as.matrix(zTilde[[i]]), varcov = SigmaCond);
      lambdai=dF(zTilde[[i]],SigmaCond,0, denominator=FALSE)/FzTilde;#general lambda
      lambdai_z=z[[i]]*lambdai;#lambda sign adjusted
      lambdai2=d2F(zTilde[[i]],SigmaCond,0)/FzTilde;#great lambda
      lambda_z[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = n_selection_equations_max);#lambda*zi
      lambda[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = n_selection_equations_max);#lambda
      lambda2[[i]]=array(0,dim=c(length(lambdai[,1]),n_selection_equations_max, n_selection_equations_max));#Big lambda
      lambda_zTilde[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = n_selection_equations_max);#lambda1*zTilde
      LAMBDA[[i]]=array(0,dim=c(length(lambdai[,1]),n_selection_equations_max, n_selection_equations_max));#Big lambda
      zG[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = n_selection_equations_max);
      zTildeG[[i]]=matrix(0,nrow = length(lambdai[,1]),ncol = n_selection_equations_max);
      #Adding zeros instead of ommited equations in lambda
      counterz=0;
      for (j in 1:n_selection_equations_max)
      {
        if (rules[i,j]!=0)
        {
          counterz=counterz+1;
          counterz1=0;
          lambda_z[[i]][,j]=lambdai_z[,counterz];
          lambda[[i]][,j]=lambdai[,counterz];
          lambda_zTilde[[i]][,j]=lambdai[,counterz]*zTilde[[i]][,counterz];
          zG[[i]][,j]=z[[i]][,counterz];
          zTildeG[[i]][,j]=zTilde[[i]][,counterz];
          #Hessian of lambda
          for (k in 1:n_selection_equations_max)
          {
            if (rules[i,k]!=0)
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
  #Combining lambda by groups into lambda by outcome
  lambda_zG=lambda_z; lambda_z=matrix(list(), n_outcome, 1);
  lambdaG=lambda; lambda=matrix(list(), n_outcome, 1);
  lambdaG2=lambda2; lambda2=matrix(list(), n_outcome, 1);
  lambda_zTildeG=lambda_zTilde; lambda_zTilde=matrix(list(), n_outcome, 1);
  LAMBDAG=LAMBDA; LAMBDA=matrix(list(), n_outcome, 1);
  zTildeOutcome=matrix(list(), n_outcome, 1);
  zOutcome=matrix(list(), n_outcome, 1);
  for (i in 1:n_groups)
  {
    if (groups[i]!=0)
    {
      lambda_z[[groups[i]]]=rbind(lambda_z[[groups[i]]],lambda_zG[[i]]);
      lambda[[groups[i]]]=rbind(lambda[[groups[i]]],lambdaG[[i]]);
      lambda2[[groups[i]]]=abind(lambda2[[groups[i]]],lambdaG2[[i]],along=1);
      lambda_zTilde[[groups[i]]]=rbind(lambda_zTilde[[groups[i]]],lambda_zTildeG[[i]]);
      LAMBDA[[groups[i]]]=abind(LAMBDA[[groups[i]]],LAMBDAG[[i]],along=1);
      zOutcome[[groups[i]]]=rbind(zOutcome[[groups[i]]],zG[[i]]);
      zTildeOutcome[[groups[i]]]=rbind(zTildeOutcome[[groups[i]]],zTildeG[[i]]);
    }
  }
  result=noquote(cbind(parameters,f$solution,stdev,pvalue));
  colnames(result)=c("Parameter","value","stdev","p-value");
  citation="Kossova, Elena & Potanin, Bogdan, 2018. 'Heckman method and switching regression model multivariate generalization,' Applied Econometrics, Publishing House 'SINERGIA PRESS', vol. 50, pages 114-143."
  print("Please cite as :");
  print(citation);
  return(list("mle"=list("result" = result, "coefficients"=f$solution,"stdev"=stdev,"p-value"=pvalue,"names"=parameters), "twostep"=list("model"=twostep,"covmatrix"=CovB,"sigma"=sigma, "twostepLS"=twostepOLS, "x0"=x0_twostep), "logLikelihood"=-f$objective, "x0"=f$solution,"lambda"=lambdaG, "sortList"=sortList, "CovM"=CovM), "citation"=citation)
}

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
#' @param show_info shows likelihood function optimization info
#' @param only_twostep if true then only two-step procedure used
#' @param opts options that are passed to nlopt
#' @param x0 optional initial values to be used while solving optimization task
#' @param remove_zero_columns set true for switch regression (more then 1 groups) if your dependent variable is a part of dummy variable
#' @details This function estimates Multivariate Switch-Selection model via maximum-likelihood and two-step procedures.
#' This model was developed by E.V. Kossova and B.S. Potanin
#' Please cite as follows:
#' Kossova, Elena & Potanin, Bogdan, 2018. 'Heckman method and switching regression model multivariate
#' generalization,' Applied Econometrics, Publishing House 'SINERGIA PRESS', vol. 50, pages 114-143.
#' Dependent variables in selection equations should have values -1 or 1.
#' Also there is a special value 0 which indicates that this observation is unobservable but it is necessary
#' to take into consideration other selection equation information while calculating likelihood function.
#' The i-th row of rules corresponds to i-th element of groups. rules rows should contain information regarding
#' possible combinations of selection equations values. While groups determines the outcome for each of this
#' possible combinations. Special 0 value for groups responsible for sample selection.
gheckman<-function(data, outcome, selection1=NULL, selection2=NULL, selection3=NULL,
                   selection4=NULL, selection5=NULL, groups=NULL, rules=NULL,
                   show_info=FALSE, only_twostep=FALSE,
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
  sort_list<-gheckmanSort(y, y_variables, z, z_variables, groups, rules, remove_zero_columns);
  #Make sort_list elements variables
  for (i in 1:length(sort_list))
  {
    if(names(sort_list)[[i]]!="z" & names(sort_list)[[i]]!="z_variables")
    {
    do.call("<-",list(names(sort_list)[[i]], sort_list[[i]]));
    }
  }
  opts = opts;#setting optimization options max(maxeval/15,n_selection_equations_max*50)
  #PHASE 2: Two-step method
  n_parameters=coef_z[[n_selection_equations_max]][dim(coef_z[[n_selection_equations_max]])[1]]
  n_parameters_y=coef_y[[n_outcome]][dim(coef_y[[n_outcome]])[1]]
  x0=matrix(0,n_parameters);
  x0_names=vector(length = length(x0));#Store variables names
  #Get initial values
  if (n_selection_equations_max==1) 
    {
      z=matrix(z,ncol=1);
    }
  for (i in 1:n_selection_equations_max)
  {
    #-1 in order to exclude constant
    x0[coef_z[[i]]]=coef(probit_simple<-glm(I((z[z[,i]!=0,i]+1)/2)~.-1,
      data=data.frame(z_variables[[i]][z[,i]!=0,]),family=binomial(link="probit")));
    x0_names[coef_z[[i]]]=all.vars(formula(probit_simple)[-2]);#get z_variables variables names
  }
  z_variables_names=matrix(list(),n_selection_equations_max);#z_variables_names
  for (i in 1:n_selection_equations_max) {z_variables_names[[i]]=colnames(z_variables[[i]]);}
  #Load z data from sort_list
  z_variables=sort_list$z_variables;#Only now we can load it from sort method result
  z=sort_list$z;#Only now we can load it from sort method result
  #Set lower and upper bound constraints for x0_names
  lb=c(rep(-0.9999999,rho_z_n),rep(0,n_parameters_y-rho_z_n),rep(-Inf,n_parameters-coef_z[[1]][1]+1));
  ub=c(rep(0.9999999,rho_z_n),rep(0,n_parameters_y-rho_z_n),rep(Inf,n_parameters-coef_z[[1]][1]+1));
  #Estimate coefficients and store them to x0
  f<-nloptr(x0=x0, eval_f=gheckmanLikelihood,opts=opts, lb=lb, ub=ub, y=y, z_variables=z_variables, 
            y_variables=y_variables, rules_no_ommit=rules_no_ommit, n_selection_equations=n_selection_equations, 
            coef_z=coef_z, sigma_last_index=sigma_last_index, coef_y=coef_y, groups=groups*0, n_groups=n_groups, 
            n_selection_equations_max=n_selection_equations_max, rules_converter=rules_converter, 
            n_outcome=n_outcome, rules=rules, groups_observations=groups_observations, 
            rho_y_indices=rho_y_indices, show_info=show_info, maximization=FALSE);
  x0=f$solution;
  #Storing covariance matrix
  covmatrix=jacobian(func = gheckmanGradient, x = x0, y=y, z_variables=z_variables, y_variables=y_variables, rules_no_ommit=rules_no_ommit, n_selection_equations=n_selection_equations, coef_z=coef_z, sigma_last_index=sigma_last_index, coef_y=coef_y, groups=groups*0, n_groups=n_groups, n_selection_equations_max=n_selection_equations_max, rules_converter=rules_converter, n_outcome=n_outcome, rules=rules, groups_observations=groups_observations, rho_y_indices=rho_y_indices, show_info=FALSE);
  covmatrix_gamma=solve(covmatrix[coef_z[[1]][1]:length(x0),coef_z[[1]][1]:length(x0)]);
  sort_list$covmatrix_gamma=covmatrix_gamma;
  #Calculating lambda
    #Storing variables names to sort_list
  sort_list$x0_names=x0_names
    #One dimensional
  sort_list=lambdaCalculate(x0=x0, sort_list=sort_list, is_2D = TRUE);
  #Twostep LS
  sort_list=gheckmanLS(x0, sort_list);
  sort_list=gheckmanLSAdjustCovariance(sort_list);
  print(summary(sort_list$twostep_LS_unadjusted[[1]]))
  print(sort_list$twostep_LS[[1]])
  print(sort_list$sigma_twostep[1]);
  #Update variables
  for (i in 1:length(sort_list))
  {
    do.call("<-",list(names(sort_list)[[i]], sort_list[[i]]));
  }
  x0[(sigma_last_index-n_outcome+1):sigma_last_index]=sigma_twostep;
  #Return the result for twostep
  if (only_twostep) 
  {
    return(list('model'=twostep_LS, 'covmatrix'=Cov_B,'sigma'=sigma_twostep, 
                'model_unadjusted'=twostep_LS_unadjusted, 'sort_list'=sort_list));
  }

  #PAHSE 3: MLE with Two-step initial values
  #Set lower and upper bound constraints for x0_names
  rho_n=sigma_last_index-n_outcome;
  no_rho_n=length(x0)-rho_n;
  lb=(c(rep(-1,rho_n),rep(-Inf,no_rho_n)));
  ub=(c(rep(1,rho_n),rep(Inf,no_rho_n)));
  x0[is.na(x0)]=0;
  sort_list$x0=x0
  #Estimate coefficients and store them to MLE
  if (!is.null(x1)) {x0<-x1;};#substitute some initial value
  f<-nloptr(x0=x0, eval_f=gheckmanLikelihood,opts=opts, lb=lb, ub=ub, y=y, z_variables=z_variables, y_variables=y_variables, rules_no_ommit=rules_no_ommit, n_selection_equations=n_selection_equations, coef_z=coef_z, sigma_last_index=sigma_last_index, coef_y=coef_y, groups=groups, n_groups=n_groups, n_selection_equations_max=n_selection_equations_max, rules_converter=rules_converter, n_outcome=n_outcome, rules=rules, groups_observations=groups_observations, rho_y_indices=rho_y_indices, show_info=show_info, maximization=FALSE);
  #Standard deviations estimation
  st_dev=jacobian(func = gheckmanGradient, x = f$solution, y=y, z_variables=z_variables, y_variables=y_variables, rules_no_ommit=rules_no_ommit, n_selection_equations=n_selection_equations, coef_z=coef_z, sigma_last_index=sigma_last_index, coef_y=coef_y, groups=groups, n_groups=n_groups, n_selection_equations_max=n_selection_equations_max, rules_converter=rules_converter, n_outcome=n_outcome, rules=rules, groups_observations=groups_observations, rho_y_indices=rho_y_indices, show_info=FALSE);
  Cov_M=solve(st_dev)
  st_dev=sqrt(diag(Cov_M));
  x0=f$solution;
  aSTD=pnorm(x0/st_dev);
  bSTD=pnorm(-x0/st_dev);
  p_value=x0*0;
  for (i in (1:length(x0)))
  {
    p_value[i]=1-max(aSTD[i],bSTD[i])+min(aSTD[i],bSTD[i]);
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
        x0_names[counter]=paste(c('rho Z',i,j), collapse = " ");
      }
    }
  }
  for (i in 1:n_outcome)
  {
    for (j in 1:n_selection_equations_max)
    {
      counter=counter+1;
      x0_names[counter]=paste(c('rho Y',i,j), collapse = " ");
    }
  }
  for (i in 1:n_outcome)
  {
    counter=counter+1;
    x0_names[counter]=paste(c('sigma',i), collapse = " ");
  }
    #MLE result table
  result=noquote(cbind(x0_names,x0,st_dev,p_value));
  colnames(result)=c("Coefficients:","Estimate","Std. Error","p-value");
  result=as.table(result)
  rownames(result)<-rep("",length(x0_names));
    #Citation
  citation="Kossova, Elena & Potanin, Bogdan, 2018. 'Heckman method and switching regression model multivariate generalization,' Applied Econometrics, Publishing House 'SINERGIA PRESS', vol. 50, pages 114-143."
  print(noquote(paste("Please cite as :",citation)));
  return(list("mle"=list("result" = result, "coefficients"=x0,"st_dev"=st_dev,"p-value"=p_value), "twostep"=list("model"=twostep_LS,"covmatrix"=Cov_B,"sigma"=sigma, "twostep_LS"=twostep_LS_unadjusted, "x0"=sort_list$x0), "logLikelihood"=-f$objective, "x0"=x0, "sort_list"=sort_list, "Cov_M"=Cov_M, "citation"=citation))
}

#' Use this function to estimate lambda
#' @param x0 values obtained from optimization routine
#' @param sort_list list of variables that gheckman method returns
#' @param groups_indices_ascending groups for which you want to estimate labmdas. Insure that they go ascending order
#' @param new_rules custom rules which correspond to different groups
gheckmanLS<-function(x0, sort_list)
{
  #Make sort_list elements variables
  for (i in 1:length(sort_list))
  {
    do.call("<-",list(names(sort_list)[[i]], sort_list[[i]]));
  }
  #Initialazing variables
  twostep_LS=matrix(list(), n_outcome, 1);#lm with adjusted covariance matrix
  n_coef_outcome=matrix(0,1,n_outcome);#number of coefficients for given outcome
  X=matrix(list(),n_outcome,1);#regressors matrix for given outcome
  rho_multiply_sigma=matrix(list(),n_outcome,1);#lambda coefficients
  twostep_errors=matrix(list(), n_outcome, 1);
  #Regression for each outcome
  for (i in 1:n_outcome)
  {
    X[[i]]=data.frame(y_variables_outcome[[i]], lambda_outcome[[i]]*z_outcome[[i]])
    twostep_LS[[i]]=lm(y_outcome[[i]]~.,data=X[[i]][,-1]);#-1 to exclude constant
    n_coef_outcome[i]=length(x0[coef_y[[i]]]);
    x0[coef_y[[i]]]=coef(twostep_LS[[i]])[1:n_coef_outcome[i]];#coefficients
    rho_multiply_sigma[[i]]=coef(twostep_LS[[i]])[(n_coef_outcome[i]+1):length(coef(twostep_LS[[i]]))];#coefficients
    #rho_multiply_sigma[is.na(coefLambda)]=0;
    x0_names[coef_y[[i]]]=variable.names(twostep_LS[[i]])[1:n_coef_outcome[i]];#Store coefficients names
    twostep_errors[[i]]=residuals(twostep_LS[[i]]);
  }
  list_return=list('twostep_LS'=twostep_LS, 'rho_multiply_sigma'=rho_multiply_sigma, 
            'x0_names'=x0_names, 'twostep_errors'=twostep_errors, 'X'=X);
  sort_list$x0=x0;
  sort_list=append(sort_list,list_return);
return(sort_list)
}
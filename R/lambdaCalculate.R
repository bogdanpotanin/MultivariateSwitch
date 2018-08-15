#' Use this function to estimate lambda
#' @param x0 values obtained from optimization routine
#' @param sort_list list of variables that gheckman method returns
#' @param groups_indices_ascending groups for which you want to estimate labmdas. Insure that they go ascending order
#' @param new_rules custom rules which correspond to different groups
lambdaCalculate<-function(x0, sort_list, new_rules=NULL, groups_indices_ascending=NULL, is_1D=TRUE, is_2D=FALSE)
{
#Make sort_list elements variables
for (i in 1:length(sort_list))
{
  do.call("<-",list(names(sort_list)[[i]], sort_list[[i]]));
}
#Assisting variables
z_tilde=matrix(list(), n_groups, 1);
Sigma=matrix(list(), n_groups, 1);
Sigma_Cond=matrix(list(), n_groups, n_groups);
Sigma_no_y=triangular(x0[1:rho_z_n],rep(1,n_selection_equations_max));
F_z_tilde=matrix(list(), n_groups, n_groups);
rho_z_indices=matrix(list(), n_outcome, 1);
#Check for custom input
if (is.null(groups_indices_ascending))
{
  groups_indices_ascending=1:n_groups;
}
groups_iterated=groups_indices_ascending;
if (!is.null(new_rules))
{
  new_rules_adjusted=matrix(c(NA),nrow=n_groups,ncol = n_selection_equations_max);
  new_rules_adjusted[groups_indices_ascending,]=new_rules;
  new_rules=new_rules_adjusted;
  new_rules_no_ommit=matrix(list(), n_groups, 1);
  new_rules_converter=matrix(list(), n_groups, 1);#rules_converter[[i]] gives rules_no_ommit i-th element index in rules
  new_rules_index=1:n_selection_equations_max;#index columns of rules
  for (i in 1:groups_iterated)
  {
    new_rules_no_ommit[[i]]=matrix(new_rules[i,new_rules[i,]!=0],nrow=1);#removing 0 from rules_no_ommit
    new_rules_converter[[i]]=new_rules_index[which(new_rules[i,] != 0)];#determining converter mapping
  }
  rules=new_rules;
  rules_no_ommit=new_rules_no_ommit;
  rules_converter=new_rules_converter;
}
#Calculating variables values
for (i in groups_iterated)
{
  Sigma[[i]]=Sigma_no_y[rules[i,]!=0,rules[i,]!=0];#remove unobservable z
  z_tilde[[i]]=z[[i]]*NA;
  if (groups[i]!=0)#if it is not unobservable group
  {
    Sigma_Cond[[i]]=as.matrix(Sigma[[i]]);
    for (j in 1:n_selection_equations[i])
    {
      Sigma_Cond[[i]][,j]=Sigma_Cond[[i]][,j]*(-rules_no_ommit[[i]][j]);
      Sigma_Cond[[i]][j,]=Sigma_Cond[[i]][j,]*(-rules_no_ommit[[i]][j]);
      z_i_j=z_variables[[i,rules_converter[[i]][j]]]%*%as.matrix(x0[coef_z[[rules_converter[[i]][j]]]]);#%predicted values for zi
      z[[i]][,j]=z[[i]][,j]*sign(z[[i]][1,j]*rules_no_ommit[[i]][j]);
      z_tilde[[i]][,j]=z_i_j*z[[i]][,j];#%upper adjusted values
    }
    F_z_tilde[[i]]=pmnorm(as.matrix(z_tilde[[i]]), varcov = Sigma_Cond[[i]]);
  }
}
#Lambdas
  #One dimensional
lambda=matrix(list(), n_groups, 1);#for groups
lambda_outcome=matrix(list(), n_outcome, 1);#for outcomes
lambda_no_ommit=matrix(list(), n_groups, 1);#for groups without 0
lambda_outcome_no_ommit=matrix(list(), n_outcome, 1);#for outcomes without 0
z_tilde_outcome=matrix(list(), n_outcome, 1);#z_tilde for outcomes
z_tilde_outcome_no_ommit=matrix(list(), n_outcome, 1);#z_tilde for outcomes
if(is_1D)
{
  for (i in (groups_iterated)[groups[groups_iterated]!=0])
  {
    lambda[[i]]=matrix(0,nrow = groups_observations[i],ncol = n_selection_equations_max);
    lambda[[i]][,rules[i,]!=0]=dF(z_tilde[[i]],Sigma_Cond[[i]],0, denominator=FALSE)/F_z_tilde[[i]];
    lambda_outcome[[groups[i]]]=rbind(lambda_outcome[[groups[i]]],lambda[[i]]);
    lambda_no_ommit[[i]]=lambda[[i]][,rules[i,]!=0];
    lambda_outcome_no_ommit[[groups[i]]]=lambda_outcome[[groups[i]]][,rules[i,]!=0];
    z_tilde_outcome[[groups[i]]]=rbind(z_tilde_outcome[[groups[i]]],z_tilde[[i]]);
    z_tilde_outcome_no_ommit[[groups[i]]]=z_tilde_outcome[[groups[i]]][,rules[i,]!=0];
  }
}
  #Two dimensional
lambda2=matrix(list(), n_groups, 1);#for groups
lambda2_outcome=matrix(list(), n_outcome, 1);#for outcomes
lambda2_no_ommit=matrix(list(), n_groups, 1);#for groups without 0
lambda2_outcome_no_ommit=matrix(list(), n_outcome, 1);#for outcomes without 0
if(is_2D)
{
  for (i in (groups_iterated)[groups[groups_iterated]!=0])
  {
  lambda2[[i]]=array(0,dim=c(groups_observations[i],n_selection_equations_max, n_selection_equations_max));
  lambda2[[i]][,rules[i,]!=0,rules[i,]!=0]=d2F(z_tilde[[i]],Sigma_Cond[[i]],0)/F_z_tilde[[i]];
  lambda2_outcome[[groups[i]]]=abind(lambda2_outcome[[groups[i]]],lambda2[[i]], along = 1);
  lambda2_no_ommit[[i]]=lambda2[[i]][,rules[i,]!=0,rules[i,]!=0];
  lambda2_outcome_no_ommit[[groups[i]]]=lambda2_outcome[[groups[i]]][,rules[i,]!=0,rules[i,]!=0];
  }
}
list_return=list('lambda'=lambda, 'lambda_outcome'=lambda_outcome, 'lambda2'=lambda2, 'lambda2_outcome'=lambda2_outcome,
            'lambda_no_ommit'=lambda_no_ommit, 'lambda_outcome_no_ommit'=lambda_outcome_no_ommit,
            'z_tilde_outcome'=z_tilde_outcome, 'z_tilde_outcome_no_ommit'=z_tilde_outcome_no_ommit,
            'lambda2_no_ommit'=lambda2_no_ommit, 'lambda2_outcome_no_ommit'=lambda2_outcome_no_ommit,
            'Sigma'=Sigma, 'Sigma_no_y'=Sigma_no_y, 'Sigma_Cond'=Sigma_Cond, 'z_tilde'=z_tilde, 'F_z_tilde'=F_z_tilde);
sort_list=append(sort_list,list_return);
return(sort_list)
}
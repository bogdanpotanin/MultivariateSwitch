#' Use this function to estimate lambda
#Calculating lambda
lambdaCalculate<-function(x0, sort_list=NULL)
{
#Assisting variables
z_tilde=matrix(list(), sort_list$n_groups, 1);
Sigma=matrix(list(), sort_list$n_groups, 1);
Sigma_Cond=matrix(list(), sort_list$n_groups, sort_list$n_groups);
Sigma_no_y=triangular(x0[1:sort_list$rho_z_n],rep(1,sort_list$n_selection_equations_max));
F_z_tilde=matrix(list(), sort_list$n_groups, sort_list$n_groups);
for (i in 1:n_groups)
{
  Sigma[[i]]=Sigma_no_y[sort_list$rules[i,]!=0,sort_list$rules[i,]!=0];#remove unobservable z
  z_tilde[[i]]=z[[i]]*NA;
  if (groups[i]!=0)#if it is not unobservable group
  {
    Sigma_Cond[[i]]=as.matrix(Sigma[[i]]);
    for (j in 1:sort_list$n_selection_equations[i])
    {
      Sigma_Cond[[i]][,j]=Sigma_Cond[[i]][,j]*(-sort_list$rules_no_ommit[[i]][j]);
      Sigma_Cond[[i]][j,]=Sigma_Cond[[i]][j,]*(-sort_list$rules_no_ommit[[i]][j]);
      z_i_j=sort_list$z_variables[[i,sort_list$rules_converter[[i]][j]]]%*%as.matrix(x0[sort_list$coef_z[[rules_converter[[i]][j]]]]);#%predicted values for zi
      z_tilde[[i]][,j]=z_i_j*sort_list$z[[i]][,j];#%upper adjusted values
    }
    F_z_tilde[[i]]=pmnorm(as.matrix(z_tilde[[i]]), varcov = Sigma_Cond);
  }
}
#Lambdas
lambda=matrix(list(), sort_list$n_groups, 1);#for groups
lambda_outcome=matrix(list(), sort_list$n_outcome, 1);#for outcomes
  for (i in 1:n_groups)
  {
    lambda[[i]]=matrix(0,nrow = sort_list$groups_observations[i],ncol = sort_list$n_selection_equations_max);#lambda
    lambda[[i]][,rules!=0]=dF(z_tilde[[i]],Sigma_Cond,0, denominator=FALSE)/F_z_tilde[[i]];#general lambda
    lambda_outcome[[groups[i]]]=rbind(lambda_outcome[[groups[i]]],lambda[[i]]);
  }
return(list(lambda, lambda_outcome, Sigma, Sigma_no_y, Sigma_Cond, z_tilde, F_z_tilde))
}
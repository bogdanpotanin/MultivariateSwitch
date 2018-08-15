#' Use this function to estimate lambda
#' @param x0 values obtained from optimization routine
#' @param sort_list list of variables that gheckman method returns
#' @param groups_indices_ascending groups for which you want to estimate labmdas. Insure that they go ascending order
#' @param new_rules custom rules which correspond to different groups
gheckmanLSAdjustCovariance<-function(sort_list)
{
  #Make sort_list elements variables
  for (i in 1:length(sort_list))
  {
    do.call("<-",list(names(sort_list)[[i]], sort_list[[i]]));
  }
  #Calculating sigma
  sigma_twostep=matrix(0,n_outcome,1);
  for (i in (1:n_groups)[groups!=0])#for each outcome
  {
    for (k in 1:n_selection_equations[i])#for each selection equation
    {
      sigma_twostep[groups[i]]=sigma_twostep[groups[i]]+
        sum(rho_multiply_sigma[[groups[i]]][k]^2*z_tilde[[i]][,k]*lambda_no_ommit[[i]][,k]);
      sigma_twostep[groups[i]]=sigma_twostep[groups[i]]+
        sum((z[[i]][,k]*rho_multiply_sigma[[groups[i]]][k]*lambda_no_ommit[[i]][,k])^2);
      for (j in 1:n_selection_equations[i])
      {
        if(k!=j)
        {
          rho_z_k_j=0;
          if(n_selection_equations[i]>1) { rho_z_k_j=Sigma[[i]][k,j]; }
          sigma_twostep[groups[i]]=sigma_twostep[groups[i]]-
            sum(rho_multiply_sigma[[groups[i]]][k]*z[[i]][,k]*
             (z[[i]][,j]*(rho_multiply_sigma[[groups[i]]][j]-rho_multiply_sigma[[groups[i]]][k]*rho_z_k_j)*
                lambda2_no_ommit[[i]][,k,j]));
        }
      }
    }
  }
  #Calculating elements by group
  rho_y=matrix(list(),n_outcome,1);
  y_variance=matrix(list(),n_outcome,1);
  triangle_matrix=matrix(list(),n_outcome,1);
  for (i in 1:n_outcome)
  {
    y_variance[[i]]=twostep_errors[[i]]^2;
    sigma_twostep[i]=sqrt((sigma_twostep[i]+sum(y_variance[[i]]))/outcome_observations[i]);
    rho_y[[i]]=rho_multiply_sigma[[i]]/sigma_twostep[i];
    triangle_matrix[[i]]=diag(1-y_variance[[i]]/sigma_twostep[i]^2);
  }
#Calculating Gamma
  G_gamma=matrix(list(),1,n_groups);
  G_gamma_outcome=matrix(list(),1,n_outcome);
  for (i in (1:n_groups)[groups!=0])#for each outcome
  {
    G_gamma[[i]]=matrix(0,groups_observations[i],sum(n_z_variables[rules[i,]!=0]));
    n_z_variables_start=1;
    n_z_variables_end=n_z_variables[rules[i,]!=0][1];
    for (k in 1:n_selection_equations[i])#for each selection equation
    {
      G_gamma_inner_part=-rho_y[[group[i]]][k]*
        (z_tilde[[i]][,k]*lambda_no_ommit[[i]][,k]+lambda_no_ommit[[i]][,k]^2)
      for (j in 1:n_selection_equations[i])
      {
        if(k!=j)
        {
          rho_z_k_j=0;
          if(n_selection_equations[i]>1) { rho_z_k_j=Sigma[[i]][k,j]; }
          G_gamma_inner_part=G_gamma_inner_part-rho_y[[group[i]]][k]*
            rho_z_k_j*lambda2_no_ommit[[i]][,k,j]*z[[i]][,k]*z[[i]][,j];
          G_gamma_inner_part=G_gamma_inner_part+rho_y[[group[i]]][j]*
            z[[i]][,k]*z[[i]][,j]*
            (lambda2_no_ommit[[i]][,k,j]-lambda_no_ommit[[i]][,k]*lambda_no_ommit[[i]][,j]);
        }
      }
      G_gamma[[i]][,n_z_variables_start:n_z_variables_end]=
        as.vector(G_gamma_inner_part)*z_variables[[i,rules_converter[[i]][k]]];
      if(k<n_selection_equations[i])#to prevent out of bounds for k+1
      {
      n_z_variables_start=n_z_variables_end+1;
      n_z_variables_end=n_z_variables_start+n_z_variables[rules[i,]!=0][k+1]-1;
      }
    }
    G_gamma_outcome[[groups[i]]]=rbind(G_gamma_outcome[[groups[i]]],G_gamma[[i]]);
  }
  #Calculating covariance matrix
  Cov_B=matrix(list(),1,n_outcome);
  twostep_LS_unadjusted=twostep_LS;
  for (i in 1:n_outcome)
  {
    #Calculating covariance matrix
  Cov_B_2=diag(outcome_observations[i])-triangle_matrix[[i]];
  Cov_B_1=G_gamma_outcome[[groups[i]]]%*%covmatrix_gamma%*%t(G_gamma_outcome[[groups[i]]]);
  Q=as.matrix(X[[i]]);
  sigma_twostep[i]=summary(twostep_LS_unadjusted[[i]])$sigma
  Cov_B[[i]]=(sigma_twostep[i]^2*(solve(t(Q)%*%Q)%*%t(Q)%*%
                                    (Cov_B_1+Cov_B_2)
                                  %*%Q%*%solve(t(Q)%*%Q)));
    #Adjust twostep least squares covariance matrix
  coef_adjusted=coeftest(twostep_LS_unadjusted[[i]], vcov. = Cov_B[[i]])
  twostep_LS[[i]]=summary(twostep_LS_unadjusted[[i]]);
  twostep_LS[[i]]$coefficients=coef_adjusted;
  twostep_LS[[i]]$sigma=sigma_twostep[i];
  }
  list_return=list('sigma_twostep'=sigma_twostep, 'twostep_LS_unadjusted'=twostep_LS_unadjusted,
                   'Cov_B'=Cov_B);
  sort_list$twostep_LS=twostep_LS;
  sort_list=append(sort_list,list_return);
return(sort_list)
}
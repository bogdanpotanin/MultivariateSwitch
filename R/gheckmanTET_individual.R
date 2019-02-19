gheckmanTET_individual<-function(sort_list, group_treatment, group_control,
                       rule_treatment=NULL, rule_control=NULL,
                       outcome_treatment=NULL, outcome_control=NULL, at_point=NULL,
                       standard_errors=TRUE, endogeneous_treatment_variable_names=NULL)
{
  #Make sort_list elements variables
  for (i in 1:length(sort_list))
  {
    do.call("<-",list(names(sort_list)[[i]], sort_list[[i]]));
  }
  #Determining rule and outcome under treatment and control
  if (is.null(rule_treatment))
  {
    rule_treatment=rules[group_treatment,]
  }
  if (is.null(rule_control))
  {
    rule_control=rules[group_control,]
  }
  if (is.null(outcome_treatment))
  {
    outcome_treatment=groups[group_treatment]
  }
  if (is.null(outcome_control))
  {
    outcome_control=groups[group_control]
  }
  #Assigning coefficients
  y_variables_names=colnames(y_variables[[outcome_treatment]]);
  coef_y_treatment_ind=coef_y[[outcome_treatment]];
  coef_y_control_ind=coef_y[[outcome_control]];
  coef_y_treatment=x0[coef_y_treatment_ind];
  coef_y_control=x0[coef_y_control_ind];
  #Calculating outcome variables observable parts
  y_variables_treatment=y_variables[[group_treatment]];
  y_treatment_observed=as.matrix(y_variables_treatment)%*%
    matrix(coef_y_treatment,ncol=1);
  y_variables_control=y_variables[[group_treatment]];
  y_variables_control[,endogeneous_treatment_variable_names]=-1*   #-1 если нужно поменять знак, заменить потом на rule
    y_variables_control[,endogeneous_treatment_variable_names];
  y_control_observed=as.matrix(y_variables_control)%*%
    matrix(coef_y_control,ncol=1);
  #Individual TET_observed
  TET_observed=y_treatment_observed-y_control_observed;
  #Calculating lambdas
  lambda_treatment=lambdaCalculate(x0 = x0,sort_list = sort_list,
                                   new_rules = rule_treatment,
                                   groups_indices_ascending = group_treatment)$lambda_no_ommit[[1]];
  lambda_control=lambdaCalculate(x0 = x0,sort_list = sort_list,
                                 new_rules = rule_control,
                                 groups_indices_ascending = group_treatment)$lambda_no_ommit[[1]];
  #Calculating epsilons
    #Epsilon treatment for treatment
  epsilon_treatment_for_treatment=matrix(x0[sigma_last_index-n_outcome+outcome_treatment]*
    matrix(x0[rho_y_indices[[outcome_treatment]]],ncol=1)*
    rule_treatment, nrow = 1);
  epsilon_treatment_for_treatment=sweep(as.matrix(lambda_treatment),MARGIN=2,
                                        as.matrix(epsilon_treatment_for_treatment),'*');
    #Epsilon control for treatment
  epsilon_control_for_treatment=matrix(x0[sigma_last_index-n_outcome+outcome_control]*
    matrix(x0[rho_y_indices[[outcome_control]]],ncol=1)*
    rule_treatment, nrow = 1);
  epsilon_control_for_treatment=sweep(as.matrix(lambda_treatment),MARGIN=2,
                                      as.matrix(epsilon_control_for_treatment),'*');
    #Epsilon control for control
  epsilon_control_for_control=matrix(x0[sigma_last_index-n_outcome+outcome_control]*
    matrix(x0[rho_y_indices[[outcome_control]]],ncol=1)*
    rule_control, nrow = 1);
  epsilon_control_for_control=sweep(as.matrix(lambda_control),MARGIN=2,
                                    as.matrix(epsilon_control_for_control),'*');
  #Individual TET_unobserved
  TET_unobserved=epsilon_treatment_for_treatment-epsilon_control_for_treatment;
  #Selection bias
  selection_bias=epsilon_control_for_treatment-epsilon_control_for_control;
  #Output for TET
  Decopmposition=data.frame("total_effect"=rowSums(cbind(TET_observed,TET_unobserved,
      selection_bias)),
      "TET_observed"=TET_observed,
      "TET_unobserved"=rowSums(TET_unobserved), "selection_bias"=rowSums(selection_bias));
  #Dealing with standard errors
  if (!standard_errors)
  {
  return(Decopmposition)
  }
    #Standard errors for TET_observed
  TET_observed_std=(
    as.matrix(y_variables_treatment)%*%
    (Cov_M[coef_y_treatment_ind,
        coef_y_treatment_ind]+
    Cov_M[coef_y_control_ind,
        coef_y_control_ind]-
        2*Cov_M[coef_y_treatment_ind,
                coef_y_control_ind])%*%
  t(as.matrix(y_variables_treatment))
)
TET_observed_std=sqrt(diag(TET_observed_std));
    #Standard errors for TET_unobserved
    #Assume that lambdas are constants
      #Calculating gradient for function of random variables
delta_gradient_rho_treatment=diag(x0[sigma_last_index-n_outcome+outcome_treatment],
                                          n_selection_equations_max,n_selection_equations_max);
delta_gradient_rho_control=diag(-x0[sigma_last_index-n_outcome+outcome_control],
                                         n_selection_equations_max,n_selection_equations_max);
delta_gradient_sigma_treatment=x0[rho_y_indices[[outcome_treatment]]];
delta_gradient_sigma_control=-1*x0[rho_y_indices[[outcome_control]]];
delta_gradient=rbind(delta_gradient_rho_treatment,delta_gradient_rho_control,
                     delta_gradient_sigma_treatment,delta_gradient_sigma_control);
delta_indecies=c(rho_y_indices[[outcome_treatment]],rho_y_indices[[outcome_control]],
                 sigma_last_index-n_outcome+outcome_treatment,sigma_last_index-n_outcome+outcome_control);
delta_matrix=t(delta_gradient)%*%Cov_M[delta_indecies,delta_indecies]%*%delta_gradient;
      #Final part
lambda_treatment_under_rule=sweep(lambda_treatment,MARGIN=2,rule_treatment,'*');
TET_unobserved_std=as.matrix(lambda_treatment_under_rule)%*%delta_matrix%*%t(as.matrix(lambda_treatment_under_rule));
TET_unobserved_std=sqrt(diag(TET_unobserved_std));
      #Standard errors for selection_bias
      #Assume that lambdas are constants
        #Calculating gradient for function of random variables
delta_gradient_2=-rbind(delta_gradient_rho_control,delta_gradient_sigma_control);
delta_indecies_2=c(rho_y_indices[[outcome_control]],sigma_last_index-n_outcome+outcome_control);
delta_matrix_2=t(delta_gradient_2)%*%Cov_M[delta_indecies_2,delta_indecies_2]%*%delta_gradient_2;
lambda_control_under_rule=sweep(lambda_control,MARGIN=2,rule_control,'*');
selection_bias_std=as.matrix(lambda_treatment_under_rule-lambda_control_under_rule)%*%delta_matrix_2%*%t(as.matrix(lambda_treatment_under_rule-lambda_control_under_rule));
selection_bias_std=sqrt(diag(selection_bias_std));
      #Standard errors for total effect
      #It would be nice to vectorize
total_effect_std=matrix(0,groups_observations,groups_observations)
delta_indecies=c(coef_y_treatment_ind,coef_y_control_ind,
                 rho_y_indices[[outcome_treatment]],rho_y_indices[[outcome_control]],
                 sigma_last_index-n_outcome+outcome_treatment,sigma_last_index-n_outcome+outcome_control);
for(i in 1:groups_observations[group_treatment])
{
delta_gradient_beta_1=as.double(y_variables_treatment[i,]);
delta_gradient_beta_2=as.double(-y_variables_treatment[i,]);
delta_gradient_rho_treatment=as.double(x0[sigma_last_index-n_outcome+outcome_treatment]*
                                         lambda_treatment_under_rule[i,]);
delta_gradient_rho_control=as.double(-x0[sigma_last_index-n_outcome+outcome_control]*
                                       lambda_control_under_rule[i,]);
delta_gradient_sigma_treatment=as.double(sum(x0[rho_y_indices[[outcome_treatment]]]*
                                               lambda_treatment_under_rule[i,]));
delta_gradient_sigma_control=as.double(sum(-1*x0[rho_y_indices[[outcome_control]]]*
                                             lambda_control_under_rule[i,]));
delta_gradient=matrix(c(delta_gradient_beta_1,delta_gradient_beta_2,
                            delta_gradient_rho_treatment,
                            delta_gradient_rho_control,
                            delta_gradient_sigma_treatment,
                            delta_gradient_sigma_control),ncol=1);
total_effect_std[i,]=t(delta_gradient)%*%Cov_M[delta_indecies,delta_indecies]%*%delta_gradient;
}
total_effect_std=sqrt(diag(total_effect_std));
#Output for TET_std
Decopmposition_std=data.frame("TET_observed_std"=TET_observed_std,
                          "TET_unobserved_std"=TET_unobserved_std, "selection_bias_std"=selection_bias_std, "total_effect_std"=total_effect_std);
TET_output=data.frame(Decopmposition,Decopmposition_std);
return(TET_output)
}

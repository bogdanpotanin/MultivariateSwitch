gheckmanATET<-function(sort_list, group_treatment, group_control,
                       rule_treatment, rule_control,
                       outcome_treatment, outcome_control)
{
  #Make sort_list elements variables
  for (i in 1:length(sort_list))
  {
    do.call("<-",list(names(sort_list)[[i]], sort_list[[i]]));
  }
  #Calculating outcome variables observable parts
  y_treatment_observed=as.matrix(y_variables[[group_treatment]])%*%
    matrix(x0[coef_y[[outcome_treatment]]],ncol=1);
  y_control_observed=as.matrix(y_variables[[group_control]])%*%
    matrix(x0[coef_y[[outcome_control]]],ncol=1);
  #Individual ATET_observed
  ATET_observed=y_treatment_observed-y_control_observed;
  #Calculating lambdas
  lambda_treatment=lambdaCalculate(x0 = x0,sort_list = sort_list,
                                   new_rules = rule_treatment, 
                                   groups_indices_ascending = group_treatment)$lambda[[1]];
  lambda_control=lambdaCalculate(x0 = x0,sort_list = sort_list,
                                 new_rules = rule_control, 
                                 groups_indices_ascending = group_treatment)$lambda[[1]];
  #Calculating epsilons
  epsilon_treatment_for_treatment=x0[sigma_last_index-outcome_treatment+1]*
    matrix(x0[rho_y_indices[[outcome_treatment]]],ncol=1)*
    rule_treatment;
  epsilon_control_for_treatment=x0[sigma_last_index-outcome_treatment+1]*
    matrix(x0[rho_y_indices[[outcome_treatment]]],ncol=1)*
    rule_treatment;
  epsilon_control_for_control=x0[sigma_last_index-outcome_control+1]*
    matrix(x0[rho_y_indices[[outcome_control]]],ncol=1)*
    rule_control;
  #Get converted selection equation index
  selection_equation_index_converted=rules_converter[[group[group_treatment_index]]][selection_equation_index];
  #Calculating outcome variable for treatment group
  y_treatment=as.matrix(y_variables[[group_treatment_index]])%*%
    matrix(x0[coef_y[[outcome_index]]],ncol=1)+
    x0[sigma_last_index-outcome_index+1]*
    as.matrix(lambda_no_ommit[[group_treatment_index]])%*%
    matrix(x0[rho_y_indices[[groups[group_treatment_index]]]],ncol=1);
  #Calculating outcome variable for control group
    #calculating alternative lambda
  rules_control=rules[group_treatment_index,];
  rules_control[selection_equation_index]=-rules_control[selection_equation_index];
  lambda_control=lambdaCalculate(x0=x0,sort_list=sort_list,new_rules = rules_control,
    groups_indices_ascending = group_treatment_index,is_2D = FALSE,
    selection_equation_remove_name=selection_equations_names[selection_equation_index],
    rules_control_selection_equation_index=rules_control[selection_equation_index])$lambda_no_ommit[[1]];
    #changing outcome variables that include selection
  index_control=grepl(selection_equations_names[selection_equation_index],
          colnames(y_variables[[group_treatment_index]]));
  y_variables_control=y_variables[[group_treatment_index]];
      #change 1 to 0 and vice versa for selection variable
  y_variables_control[,index_control]=y_variables_control[,index_control]*
    (rules_control[selection_equation_index]+1)/2;#from -1,1 to 0,1
  coef_y_control=x0[coef_y[[outcome_index]]];
    #calculating outcome variable itself
  x0[rho_y_indices[[outcome_index]]][selection_equation_index_converted]=
    -x0[rho_y_indices[[outcome_index]]][selection_equation_index_converted]
  print(x0[rho_y_indices[[outcome_index]]][selection_equation_index_converted])
  y_control=(as.matrix(y_variables_control)%*%
    matrix(coef_y_control,ncol=1)+
    x0[sigma_last_index-outcome_index+1]*
    as.matrix(lambda_control)%*%
    matrix(x0[rho_y_indices[[outcome_index]]],ncol=1));
  return(y_treatment-y_control)
}
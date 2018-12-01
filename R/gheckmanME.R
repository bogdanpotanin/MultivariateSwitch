gheckmanME<-function(sort_list, selection_equation_index, group_treatment_index, outcome_index=NULL)
{
  #Make sort_list elements variables
  for (i in 1:length(sort_list))
  {
    do.call("<-",list(names(sort_list)[[i]], sort_list[[i]]));
  }
  #Set group outcome index
  if(is.null(outcome_index))
  {
    outcome_index=groups[group_treatment_index];
  }
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
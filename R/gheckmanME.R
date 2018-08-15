gheckmanME<-function(sort_list, selection_equation_index, group_treatment_index, group_control_index, interactions_list=NULL)
{
  #Make sort_list elements variables
  for (i in 1:length(sort_list))
  {
    do.call("<-",list(names(sort_list)[[i]], sort_list[[i]]));
  }
}
 
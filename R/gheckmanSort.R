gheckmanSort<-function(y, y_variables, z, z_variables, groups=NULL, rules=NULL, remove_zero_columns=FALSE)
{
  #PHASE 1: Preparing data
  #Assigning main variables
  n_observations_total=NROW(y);#sample size
  n_selection_equations_max=NCOL(z);#total number of selection equations
  z[is.na(z)]=0;#ommited values conversion to 0 code
  y_and_z=data.frame(y,z);#all dependend variables data frame
    #decoding y into observed and not
  y_and_z_rules=y_and_z;
  y_and_z_rules$y[!is.na(y_and_z_rules$y)]=1;
  y_and_z_rules$y[is.na(y_and_z_rules$y)]=0;
  #If user has not predefined selection mechanism
  if (is.null(groups) || is.null(rules))
  {
    rules_with_y=as.matrix(unique(y_and_z_rules));#left just unique combinations
    rules_with_y[rowSums(abs(rules_with_y))==0,]=NA;#remove groupss without observations
    rules_with_y=na.omit(rules_with_y);
    groups=matrix(1,1,NROW(rules_with_y));#creating groups variable
    groups[rules_with_y[,1]==0]=0;#code 0 where is no observations for y
    rules=rules_with_y[,-1]#removing column corresponding to y
    #Here may be the problem when for the same rules y could be either observable or not
  }
  else
    {
      groups=matrix(as.vector(groups),1,NROW(rules));
      rules_with_y=as.matrix(cbind(t(groups),rules));
      rules=as.matrix(rules_with_y[,-1]);
    }
  rules_with_y[rules_with_y>1]=1
  n_groups=length(groups);#number of groups
  n_outcome=max(groups);#number of outcomes
  if (n_outcome==1) #if no switch regression
    {
      y_variables=matrix(list(y_variables))
    }
  else
  {
    y_variables_copy=y_variables;
    y_variables=matrix(list(),1,n_outcome);
    for (i in 1:n_outcome)
    {
      #each group contains all observations
      y_variables[[i]]=y_variables_copy;
    }
  }
  #Creating rules_no_ommit variable
  rules_no_ommit=matrix(list(), n_groups, 1);
  rules_converter=matrix(list(), n_groups, 1);#rules_converter[[i]] gives rules_no_ommit i-th element index in rules
  rules_index=1:n_selection_equations_max;#index columns of rules
  for (i in 1:n_groups)
  {
    rules_no_ommit[[i]]=matrix(rules[i,rules[i,]!=0],nrow=1);#removing 0 from rules_no_ommit
    rules_converter[[i]]=rules_index[which(rules[i,] != 0)];#determining converter mapping
  }
  n_selection_equations=as.matrix(rep(0,n_groups));#number of selection equations for groups i
  for (i in 1:n_groups)
  {
    n_selection_equations[i]=NCOL(rules_no_ommit[[i]]);#equals to the number of presisting equations
  }
  #PHASE 2: Indexing data
  #Defining sigma_last_index
  sigma_last_index<-sum(1:(n_selection_equations_max-1))*(n_selection_equations_max>1)+
    n_outcome*(n_selection_equations_max+1);#index of the last sigma in parameter vector
  #Indexing rho in parameter vector
  start_coef<-sum(1:(n_selection_equations_max-1))*(n_selection_equations_max>1)+1;
  rho_y_indices=matrix(list(), n_outcome, 1);
  for (i in 1:n_outcome)
  {
    end_coef=start_coef+n_selection_equations_max-1;
    rho_y_indices[[i]]=start_coef:end_coef;
    start_coef=end_coef+1;
  }
  #Names for different y_variables regressors
  y_variables_names=matrix(list(), n_outcome, 1);
  for (i in 1:n_outcome) #for each outcome type
  {
    y_variables_names[[i]]=paste('y_variables',toString(i));
    y_and_z$new=y_variables[[i]];
    colnames(y_and_z)[ncol(y_and_z)]=y_variables_names[[i]];
  }
  #Names for different z_variables regressors
  z_variables_names=matrix(list(), n_selection_equations_max, 1);
  for (i in 1:n_selection_equations_max) #for each selection eqation
  {
    z_variables_names[[i]]=paste('z_variables',toString(i));
    y_and_z$new=z_variables[[i]];
    colnames(y_and_z)[ncol(y_and_z)]=z_variables_names[[i]];
  }
  n_z_variables=as.matrix(rep(0,n_selection_equations_max));#number of variables in z_variables[[i]]
  for (i in (1:n_selection_equations_max))
  {
    n_z_variables[i]=NCOL(z_variables[[i]]);
  }
  #PHASE 3: grouping variables
  #Clearing in order to split among groups
  y_variables=matrix(list(), n_groups, 1)
  z_variables=matrix(list(), n_groups, n_selection_equations_max)
  z=matrix(list(), n_groups, 1)
  z_with_ommit=matrix(list(), n_groups, 1)
  y=matrix(list(), n_groups, 1)
  #Sorting according to rules
  y_and_z$sort_rank=rep(0,n_observations_total);#creating variables ranking
  for (i in 1:n_groups)
  {
    #Takes only each n_selection_equations_max element of ismember because of duplicates
    #y_and_z$sort_rank[apply(y_and_z_rules[,2:(n_selection_equations_max+1)],1,function(x) all(x==rules_with_y[i,]))]=i;
    y_and_z$sort_rank[apply(y_and_z_rules,1,function(x) all(x==rules_with_y[i,]))]=i;
  }
  y_and_z_rules<-NULL#free memmory
  #Sorting according to ranking variable
  y_and_z=y_and_z[with(y_and_z, order(sort_rank)),];
  #Division of z and y_variables according to groups
  n_y_variables=rep(0,n_outcome);#number of variables in y_variables[[i]]
    #division of y_variables
  for (i in 1:n_groups)
  {
    y[[i]]=as.matrix(y_and_z$y[y_and_z$sort_rank==i],ncol=1);
    colnames(y[[i]])=colnames(y_and_z)[1];
    if (groups[i]!=0)
    {
      y_variables[[i]]=y_and_z[y_variables_names[[groups[i]]]][y_and_z$sort_rank==i,];
      if(remove_zero_columns)
        {
          y_variables[[i]]=y_variables[[i]][,which(colSums(abs(y_variables[[i]]), na.rm = TRUE)>0)]
        }
      n_y_variables[groups[i]]=NCOL(y_variables[[i]])
    }
      #division of z_variables
    counter_z=0;
    if (n_selection_equations_max>1)
    {
      z_with_ommit[[i]]=y_and_z[,2:(n_selection_equations_max+1)][y_and_z$sort_rank==i,];
      z[[i]]=matrix(ncol=n_selection_equations[i],nrow=dim(z_with_ommit[[i]])[1]);#preinitialization
    }
    else
    {
      z_with_ommit[[i]]=matrix(y_and_z[,2:(n_selection_equations_max+1)][y_and_z$sort_rank==i],ncol=1);
      z[[i]]=matrix(ncol=n_selection_equations[i],nrow=length(z_with_ommit[[i]]));#preinitialization
    }
    for (j in 1:n_selection_equations_max)
    {
      if (rules[i,j]!=0)
      {
        counter_z=counter_z+1;
        z[[i]][,counter_z]<-z_with_ommit[[i]][,j];
      }
      z_variables[[i,j]]=as.matrix(y_and_z[z_variables_names[[j]]][y_and_z$sort_rank==i,]);#as matrix because else treated as list
    }
  }
  #Assigning parameter vector indexes
    #coefficients for y predictors depend on outcome
  coef_y=matrix(list(), n_outcome, 1);#set n_y_variables+1 if constant include automatically
  start_coef=sigma_last_index+1;
  for (i in (1:n_outcome))
  {
    end_coef=start_coef+n_y_variables[i]-1;
    coef_y[[i]]=as.matrix(start_coef:end_coef);
    start_coef=end_coef+1;
  }
    #coefficients for z predictors
  coef_z=matrix(list(), n_selection_equations_max, 1);#+1 greater then amount of parameters
  coef_z[[1]]=as.matrix((end_coef+1):(end_coef+n_z_variables[1]));
  if (n_selection_equations_max>1)
  {
  for (i in 2:n_selection_equations_max)
  {
    end_value=dim(coef_z[[i-1]])[1];
    coef_z[[i]]=as.matrix((coef_z[[i-1]][end_value]+1):(coef_z[[i-1]][end_value]+n_z_variables[i]));
  }
  }
  #Calculating groups sizes
  groups_observations=rep(0,n_groups);
  for (i in 1:n_groups)
  {
    groups_observations[i]=NROW(y[[i]]);
  }
  #Number of correlations between selection equations
  rho_z_n=sum(1:(n_selection_equations_max-1));
  #Combining y and y_variables by outcome
  y_outcome=matrix(list(), n_outcome, 1);
  y_variables_outcome=matrix(list(), n_outcome, 1);
  outcome_observations=matrix(0, n_outcome, 1);
  for (i in 1:n_groups)
  {
    if (groups[i]!=0)
    {
      y_outcome[[groups[i]]]=rbind(y_outcome[[groups[i]]],y[[i]]);
      y_variables_outcome[[groups[i]]]=rbind(y_variables_outcome[[groups[i]]],y_variables[[i]]);
      outcome_observations[groups[i]]=length(y_outcome[[groups[i]]]);
    }
  }
  #Combining z and z_variables by outcome
  z_outcome=matrix(list(), n_outcome, 1);
  z_variables_outcome=matrix(list(), n_outcome, n_selection_equations);
  for (i in 1:n_groups)
  {
    if (groups[i]!=0)
    {
      z_pseudo=matrix(0,groups_observations[i],n_selection_equations_max);
      z_pseudo[,rules[i,]!=0]=z[[i]];
      z_outcome[[groups[i]]]=rbind(z_outcome[[groups[i]]],z_pseudo);
      z_variables_outcome[[groups[i]]]=rbind(z_variables_outcome[[groups[i]]],z_variables[[i]]);
    }
  }
  #Rho z indices
  #Indexing rho in parameter vector
  return(list('y'=y, 'y_variables'=y_variables, 'z'=z, 'z_variables'=z_variables,
              'y_outcome'=y_outcome, 'y_variables_outcome'=y_variables_outcome,
              'z_outcome'=z_outcome, 'z_variables_outcome'=z_variables_outcome,
              'z_with_ommit'=z_with_ommit,
              'outcome_observations'=outcome_observations,
              'groups'=groups, 'rules'=rules, 'n_groups'=n_groups, 'n_outcome'=n_outcome,
              'rules_converter'=rules_converter, 'rules_no_ommit'=rules_no_ommit, 'n_selection_equations'=n_selection_equations, 'sigma_last_index'=sigma_last_index, 'n_selection_equations_max'=n_selection_equations_max, 'coef_y'=coef_y,
              'coef_z'=coef_z, 'groups_observations'=groups_observations, 'n_y_variables'=n_y_variables, 'n_z_variables'=n_z_variables, 'rho_y_indices'= rho_y_indices, 'rho_z_n'=rho_z_n))
}
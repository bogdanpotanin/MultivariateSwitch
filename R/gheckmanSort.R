gheckmanSort<-function(y, yh, z, zh, group=NULL, zo3=NULL, removeZeroColumns=FALSE)
{
  #PHASE 1: Preparing data
  #Assigning main variables
  n=NROW(y);#sample size
  nsMax=NCOL(z);#total number of selection equations
  z[is.na(z)]=0;#ommited values conversion to 0 code
  h=data.frame(y,z);#dependend variables data frame
  h1=h; h1$y[!is.na(h1$y)]=1; h1$y[is.na(h1$y)]=0;#decoding y into observed and not
  #If user has not predefined selection mechanism
  if (is.null(group) || is.null(zo3))
  {
    zo3Full=as.matrix(unique(h1));#left just unique combinations
    zo3Full[rowSums(abs(zo3Full))==0,]=NA;#remove groups without observations
    zo3Full=na.omit(zo3Full);
    group=matrix(1,1,NROW(zo3Full));#creating group variable
    group[zo3Full[,1]==0]=0;#code 0 where is no observations for y
    zo3=zo3Full[,-1]#removing column corresponding to y
    #Here may be the problem when for the same zo3 y could be either observable or not
  }
  else
    {
      group=matrix(as.vector(group),1,NROW(zo3));
      zo3Full=as.matrix(cbind(t(group),zo3));
      zo3=as.matrix(zo3Full[,-1]);
    }
  print(zo3)
  zo3Full[zo3Full>1]=1
  print(zo3Full)
  print(zo3)
  print(group)
  ngroup=length(group);#number of groups
  noutcome=max(group);#number of outcomes
  if (noutcome==1) {yh=matrix(list(yh))}#if no switch regression
  else
  {
    yh1=yh; yh=matrix(list(),1,noutcome)
    for (i in 1:noutcome)
    {
      yh[[i]]=yh1;
    }
  }
  #Creating zo variable
  zo=matrix(list(), ngroup, 1);
  zo3Converter=matrix(list(), ngroup, 1);#Zo3Converter[[i]] gives zo i-th element index in zo3
  zo3index=1:nsMax;#index columns of zo3
  for (i in 1:ngroup)
  {
    zo[[i]]=matrix(zo3[i,zo3[i,]!=0],nrow=1);#removing 0 from zo
    zo3Converter[[i]]=zo3index[which(zo3[i,] != 0)];
  }
  ns=as.matrix(rep(0,ngroup));#number of selection equations for group i
  for (i in 1:ngroup)
  {
    ns[i]=NCOL(zo[[i]]);#equals to the number of presisting equations
  }
  #PHASE 2: Indexing data
  #Defining nSigma
  nSigma<-sum(1:(nsMax-1))*(nsMax>1)+noutcome*(nsMax+1);#index of the last sigma in parameter vector
  #indexing rho in parameter vector
  startCoef<-sum(1:(nsMax-1))*(nsMax>1)+1;
  nrhoY=matrix(list(), noutcome, 1);
  for (i in 1:noutcome)
  {
    endCoef=startCoef+nsMax-1;
    nrhoY[[i]]=startCoef:endCoef;
    startCoef=endCoef+1;
  }
  #Names for different yh predictors
  yhstr=matrix(list(), noutcome, 1);
  for (i in 1:noutcome)#for each outcome type
  {
    yhstr[[i]]=paste('yh',toString(i));
    #MAYBE SHORTEN SOMEHOW
    h$new=yh[[i]];
    colnames(h)[ncol(h)]=yhstr[[i]];
  }
  #names for different zh predictors
  zhstr=matrix(list(), nsMax, 1);
  for (i in 1:nsMax)#for each selection eqation
  {
    zhstr[[i]]=paste('zh',toString(i));
    h$new=zh[[i]];
    colnames(h)[ncol(h)]=zhstr[[i]];
  }
  nzh=as.matrix(rep(0,nsMax));#number of variables in zh[[i]]
  for (i in (1:nsMax))
  {
    nzh[i]=NCOL(zh[[i]]);
  }
  #PHASE 3: Grouping variables
  #Clearing in order to split among groups
  yh=matrix(list(), ngroup, 1)
  zh=matrix(list(), ngroup, nsMax)
  z=matrix(list(), ngroup, 1)
  y=matrix(list(), ngroup, 1)
  #Sorting according to zo3
  h$sortRank=rep(0,n);#Creating variable ranking
  for (i in 1:ngroup)
  {
    #Takes only each nsMax element of ismember because of duplicates
    #h$sortRank[apply(h1[,2:(nsMax+1)],1,function(x) all(x==zo3Full[i,]))]=i;
    h$sortRank[apply(h1,1,function(x) all(x==zo3Full[i,]))]=i;
  }
  h1<-NULL#free memmory
  #Sorting according to ranking variable
  h=h[with(h, order(sortRank)),];
  #Division of z and yh according to groups
  nyh=rep(0,noutcome);#number of variables in yh[[i]]
  for (i in 1:ngroup)
  {
    y[[i]]=as.matrix(h$y[h$sortRank==i]);
    if (group[i]!=0)
    {
      yh[[i]]=h[yhstr[[group[i]]]][h$sortRank==i,];
      if(removeZeroColumns) {yh[[i]]=yh[[i]][,which(colSums(abs(yh[[i]]), na.rm = TRUE)>0)]}
      nyh[group[i]]=NCOL(yh[[i]])
    }
    counterz=0;
    if (nsMax>1)
    {
      zPseudo=h[,2:(nsMax+1)][h$sortRank==i,];
      z[[i]]=matrix(ncol=ns[i],nrow=dim(zPseudo)[1]);#Preinitialization
    }

    else
    {
      zPseudo=matrix(h[,2:(nsMax+1)][h$sortRank==i],ncol=1);
      z[[i]]=matrix(ncol=ns[i],nrow=length(zPseudo));#Preinitialization
    }
    for (j in 1:nsMax)
    {
      if (zo3[i,j]!=0)
      {
        counterz=counterz+1;
        z[[i]][,counterz]<-zPseudo[,j];
        zh[[i,j]]=as.matrix(h[zhstr[[j]]][h$sortRank==i,]);#as matrix because else treated as list
      }
    }
  }
  #Assigning parameter vector indexes
  #coefficients for y predictors depend on outcome
  coef=matrix(list(), noutcome, 1);#set nyh+1 if constant include automatically
  startCoef=nSigma+1;
  for (i in (1:noutcome))
  {
    endCoef=startCoef+nyh[i]-1;
    coef[[i]]=as.matrix(startCoef:endCoef);
    startCoef=endCoef+1;
  }
  #coefficients for z predictors
  ndz=matrix(list(), nsMax, 1);#+1 greater then amount of parameters
  ndz[[1]]=as.matrix((endCoef+1):(endCoef+nzh[1]));
  if (nsMax>1)
  {
  for (i in 2:nsMax)
  {
    endValue=dim(ndz[[i-1]])[1];
    print(nzh[i])
    print(nzh)
    ndz[[i]]=as.matrix((ndz[[i-1]][endValue]+1):(ndz[[i-1]][endValue]+nzh[i]));
  }
  }
  #Calculating groups sizes
  groupsize=rep(0,ngroup);
  for (i in 1:ngroup)
  {
    groupsize[i]=NROW(y[[i]]);
  }
  return(list('y'=y, 'yh'=yh, 'z'=z, 'zh'=zh, 'group'=group, 'zo3'=zo3, 'ngroup'=ngroup, 'noutcome'=noutcome,
              'zo3Converter'=zo3Converter, 'zo'=zo, 'ns'=ns, 'nSigma'=nSigma, 'nsMax'=nsMax, 'coef'=coef,
              'ndz'=ndz, 'groupsize'=groupsize, 'nyh'=nyh, 'nzh'=nzh, 'nrhoY'=nrhoY))
}

gheckmanLikelihood1=memoise(gheckmanLikelihood)
gheckmanF=function(x, y, zh, yh, zo, ns, ndz, nSigma, coef, group, ngroup, nsMax, zo3Converter, noutcome, zo3, groupsize, nrhoY, ShowInfo=TRUE, maximization=TRUE)
{
  return(gheckmanLikelihood1(x, y, zh, yh, zo, ns, ndz, nSigma, coef, group, ngroup, nsMax, zo3Converter, noutcome, zo3, groupsize, nrhoY, ShowInfo, maximization)[[1]])
}
gheckmanG=function(x, y, zh, yh, zo, ns, ndz, nSigma, coef, group, ngroup, nsMax, zo3Converter, noutcome, zo3, groupsize, nrhoY, ShowInfo=TRUE, maximization=TRUE)
{
  return(gheckmanLikelihood1(x, y, zh, yh, zo, ns, ndz, nSigma, coef, group, ngroup, nsMax, zo3Converter, noutcome, zo3, groupsize, nrhoY, ShowInfo, maximization)[[2]])
}
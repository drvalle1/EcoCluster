#' Samples z
#' 
#' This function samples the cluster assignment for each location (z)
#' 
#' @param dat matrix with L rows (e.g., locations) and S columns (e.g., species),
#'            containing the presence-absence data 
#' @param one.minus.dat matrix with L rows (e.g., locations) and S columns (e.g., species),
#'                      calculated as 1-dat           
#' @param phi K x S matrix with the probability of observing each species in each group
#' @param theta vector of length K with the probability of each location group
#' @param ngroup  maximum number of location groups (K)
#' @param nloc number of locations (L)
#' @param nspp number of species (S)
#' @param z vector of size L containing the current cluster assignment for each location
#' @param return this function returns a vector of size L with the cluster assignment of each location
#' @export

update.z=function(dat,one.minus.dat,phi,theta,ngroup,nloc,nspp,z){
  #pre-calculate some useful quantities
  log.theta=log(theta)
  log.phi=log(phi)
  log.one.minus.phi=log(1-phi)
  
  #determine the number of locations in each group
  tab=rep(0,ngroup)
  tmp=table(z)
  tab[as.numeric(names(tmp))]=tmp
  
  #calculate the log probability for each group that already exists
  tmp=matrix(NA,nloc,ngroup)
  for (i in 1:ngroup){
    rasc=dat*matrix(log.phi[i,],nloc,nspp,byrow=T)+
         one.minus.dat*matrix(log.one.minus.phi[i,],nloc,nspp,byrow=T)
    tmp[,i]=rowSums(rasc)+log.theta[i] #sum log of prior probability
  }
  
  #sample z
  for (i in 1:nloc){
    tab[z[i]]=tab[z[i]]-1
    prob=rep(NA,ngroup)
    cond=tab==0
    prob[ cond]=-nspp*log(2)+log.theta[cond] #log probability for a new group
    prob[!cond]=tmp[i,!cond]

    #get normalized probs
    tmp1=prob-max(prob) #for numerical stability
    tmp2=exp(tmp1) #exponentiate log probability
    prob=tmp2/sum(tmp2) #normalize to sum to 1

    #draw from multinomial distrib
    ind=rmultinom(1,size=1,prob=prob)
    ind1=which(ind==1)
    z[i]=ind1
    tab[ind1]=tab[ind1]+1
  }
  z  
}
#--------------------------------------------

#' Samples theta and v parameters
#' 
#' This function samples the v parameters, which are then used to calculate the theta parameters
#' 
#' @param z vector of size L with cluster assignment of each location 
#' @param ngroup maximum number of location groups (K)
#' @param gamma1 this is the truncated stick-breaking prior parameter for the 
#'                number of location groups. This value should be between 0 and 1, and
#'                small values enforce more parsimonius results (i.e., fewer groups)
#' @param burnin number of iterations to drop as part of burn-in phase
#' @param gibbs.step current iteration of the gibbs sampler
#' @param phi K x S matrix with the probability of observing each species in each group
#' @param theta vector of length K with the probability of each location group
#' @param return this function returns a list of 4 items (theta, z, v, and phi)
#' @export
#' 
update.theta=function(z,ngroup,gamma1,burnin,gibbs.step,theta,phi){
  #re-order thetas in decreasing order if in burn-in phase. 
  #Based on this re-ordering, re-order z's and phi's
  if(gibbs.step<burnin & gibbs.step%%50==0){
    ind=order(theta,decreasing=T)
    theta=theta[ind]
    phi=phi[ind,]
    
    #get z.new
    z.new=z; z.new[]=NA
    for (i in 1:ngroup){
      cond=z==ind[i]
      z.new[cond]=i
    }
    z=z.new
  }

  #calculate the number of locations assigned to each group
  tmp=table(z)
  nk=rep(0,ngroup)
  nk[as.numeric(names(tmp))]=tmp

  #sample v from a beta distribution
  n.greater.k=cumsum(nk[ngroup:1])[ngroup:1]
  v=rbeta(ngroup-1,nk[-ngroup]+1,n.greater.k[-1]+gamma1)
  v1=c(v,1)
  
  #get theta from v1 using the stick-breaking equation 
  theta=rep(NA,ngroup)
  tmp=1
  for (i in 1:ngroup){
    theta[i]=v1[i]*tmp
    tmp=tmp*(1-v1[i])
  }
  
  #to avoid numerical issues
  cond=v>0.99999
  v[cond]=0.99999
  
  #output results
  list(theta=theta,z=z,v=v,phi=phi)
}

#----------------------------
#' Sample the TSB prior parameter
#' 
#' This function samples the truncated stick breaking (TSB) prior parameter gamma
#' 
#' @param v vector of length L with probabilities 
#' @param ngroup maximum number of location groups (K)
#' @param gamma.possib vector of possible gamma parameter values
#' @param return this function returns a real number corresponding to gamma
#' @export
#' 
#' 
sample.gamma=function(v,ngroup,gamma.possib){
  #calculate the log probability associated with each possible value of gamma
  ngamma=length(gamma.possib)
  soma=sum(log(1-v[-ngroup]))
  k=(ngroup-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma
  #check this code: sum(dbeta(v[-ngroup],1,gamma.possib[5],log=T))
  
  #exponentiate and normalize probabilities
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  
  #sample from a categorical distribution
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}
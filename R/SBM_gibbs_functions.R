#' Sample psi parameters
#' 
#' Sample the presence probability for each location and species group (psi) from 
#' it's full conditional distribution
#' 
#' @param dat L x S matrix containing the presence-absence data 
#' @param ngroup.loc maximum number of location groups (KL)
#' @param ngroup.spp maximum number of species groups (KS)
#' @param z vector of length L, storing the current membership of each location
#' @param w vector of length S, storing the current membership of each species
#' @param return a KL x KS matrix of psi parameters
#' @export

sample.psi=function(z,w,dat,ngroup.loc,ngroup.spp){
  #summarize the data by calculating how many observations were assigned to each group
  #for which dat[i,j]=1 and for which dat[i,j]=0
  tmp=getql(z=z-1,w=w-1,dat=dat,ngrloc=ngroup.loc,ngrspp = ngroup.spp)
  
  #generate psi from beta distribution
  tmp1=rbeta(ngroup.loc*ngroup.spp,tmp$nql1+1,tmp$nql0+1)
  matrix(tmp1,ngroup.loc,ngroup.spp)
}

#' Sample theta and vk parameters
#' 
#' Sample vk parameters from their full conditional distributions and convert 
#' these vk parameters into theta parameters
#'
#' @param ngroup.loc maximum number of groups for locations (KL)
#' @param gamma.v truncated stick-breaking prior parameter for the 
#'                location groups. This value should be between 0 and 1, and
#'                small values enforce more parsimonius results (i.e., fewer groups)
#' @param burnin number of MCMC samples that are going to be thrown out as 
#'               part of the burn-in phrase
#' @param gibbs.step current iteration of the gibbs sampler
#' @param theta vector of size KL containing the current estimate of theta (i.e., probability of each location group)                  
#' @param psi KL x KS matrix  containing the current estimate of psi
#' @param z this is a vector of length L, storing the current membership of each location
#' @param return this function returns a list containing 4 items (theta, vk, psi, and z)
#' @export

sample.theta=function(ngroup.loc,gamma.v,burnin,gibbs.step,theta,psi,z){
  #re-order thetas in decreasing order if we are still in the burn-in phrase. 
  #based on this re-ordering, re-order z's and psi's
  if(gibbs.step<burnin & gibbs.step%%50==0){
    ind=order(theta,decreasing=T)
    theta=theta[ind]
    psi=psi[ind,]
    
    #get z.new
    z.new=z; z.new[]=NA
    for (i in 1:ngroup.loc){
      cond=z==ind[i]
      z.new[cond]=i
    }
    z=z.new
  }
  
  #get the number of locations in each group
  nk=rep(0,ngroup.loc)
  tmp=table(z)
  nk[as.numeric(names(tmp))]=tmp
  
  #sample vk from a beta distribution
  ind=ngroup.loc:1
  invcumsum=cumsum(nk[ind])[ind]
  vk=rbeta(ngroup.loc,nk+1,invcumsum-nk+gamma.v)
  vk[vk>0.99999]=0.99999 #to avoid numerical issues
  vk[ngroup.loc]=1
  
  #convert from vk to theta
  theta=convertSBtoNormal(vk)
  
  #output vk, theta, z, and psi
  list(vk=vk,theta=theta,z=z,psi=psi)
}

#' Sample phi and uk parameters
#' 
#' Sample uk parameters from their full conditional distributions and 
#' then calculate the implied phi parameters
#'
#' @param ngroup.spp maximum number of groups for species (KS)
#' @param gamma.u the truncated stick-breaking prior parameter for 
#'                species groups. This value should be between 0 and 1, and
#'                small values enforce more parsimonius results (i.e., fewer groups)
#' @param burnin number of MCMC samples that are going to be thrown out as 
#'               part of the burn-in phrase
#' @param gibbs.step current iteration of the gibbs sampler
#' @param phi vector of size KS containing the current estimate of phi (i.e., probability of each species group)                  
#' @param psi matrix of size KL x KS containing the current estimate  of psi
#' @param w vector of length S, storing the current membership of each species
#' @param return this function returns a list with 4 items (phi, uk, w, psi)
#' @export
#' 

sample.phi=function(ngroup.spp,gamma.u,burnin,gibbs.step,phi,psi,w){
  #re-order phi in decreasing order if we are still in burn-in phase 
  #Based on this re-ordering, re-order w's and psi's
  if(gibbs.step<burnin & gibbs.step%%50==0){
    ind=order(phi,decreasing=T)
    phi=phi[ind]
    psi=psi[,ind]
    
    #get z.new
    w.new=w; w.new[]=NA
    for (i in 1:ngroup.spp){
      cond=w==ind[i]
      w.new[cond]=i
    }
    w=w.new
  }
  
  #calculate the number of species in each group
  mk=rep(0,ngroup.spp)
  tmp=table(w)
  mk[as.numeric(names(tmp))]=tmp
  
  #sample uk from a beta distribution
  ind=ngroup.spp:1
  invcumsum=cumsum(mk[ind])[ind]
  uk=rbeta(ngroup.spp,mk+1,invcumsum-mk+gamma.u)
  uk[uk>0.9999999]=0.9999999 #for numerical issues
  uk[ngroup.spp]=1

  #convert uk to phi
  phi=convertSBtoNormal(uk)
  
  #retunr uk, phi, w, and psi
  list(uk=uk,phi=phi,w=w,psi=psi)
}
#--------------------------
#' Sample gamma.u
#' 
#' Sample the TSB prior parameter for species groups (gamma.u) from its full conditional distribution
#'
#' @param uk vector of size KS with probabilities
#' @param ngroup.spp maximum number of species groups (KS)
#' @param gamma.possib this is a vector containing the possible values that gamma.u can take
#' @param return this function returns a real number (gamma.u) 
#' @export
#' 
sample.gamma.u=function(uk,gamma.possib,ngroup.spp){
  #calculate the stick-breaking probabilities for different values of gamma.u
  ngamma=length(gamma.possib)
  soma=sum(log(1-uk[-ngroup.spp]))
  k=(ngroup.spp-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma
  #to check this code: sum(dbeta(v[-ngroup.spp],1,gamma.possib[5],log=T))
  
  #exponentiate and normalize these probabilities to draw from categorical distribution
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}
#--------------------------
#' Sample gamma.v
#' 
#' Sample the TSB prior parameter for location groups (gamma.v) from its full conditional distribution
#'
#' @param vk vector of size KL with probabilities
#' @param ngroup.loc maximum number of location groups (KL)
#' @param gamma.possib vector containing the possible values that gamma.u can take
#' @param return this function returns a real number (gamma.v) 
#' @export
#' 

sample.gamma.v=function(vk,gamma.possib,ngroup.loc){
  #calculate the stick-breaking probabilities for different values of gamma.v
  ngamma=length(gamma.possib)
  soma=sum(log(1-vk[-ngroup.loc]))
  k=(ngroup.loc-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma
  #to check this code: sum(dbeta(v[-ngroups],1,gamma.possib[5],log=T))
  
  #exponentiate and normalize these probabilities to draw from categorical distribution
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}
#--------------------------
#' Sample z
#' 
#' Sample the vector z containing the membership of each location
#'
#' @param ltheta equal to log(theta)
#' @param dat L x S matrix, containing the presence-absence data
#' @param dat1m L x S matrix, calculated as 1-dat
#' @param lpsi equal to log(psi)
#' @param l1mpsi equal to log(1-psi)
#' @param ngroup.loc maximum number of location groups (KL)
#' @param ngroup.spp maximum number of species groups (KS)
#' @param nloc total number of locations (L)
#' @param nspp total number of species (S)
#' @param w vector of length S, storing the current membership of each species
#' @param z vector of length L, storing the current membership of each locations
#' @param return this function returns a vector of length L containing z 
#' @export
#' 

sample.z=function(ltheta,dat,dat1m,lpsi,l1mpsi,ngroup.loc,ngroup.spp,nloc,nspp,w,z){
  #calculation of log probability for groups that already exist
  lprob.exist=matrix(NA,nloc,ngroup.loc)
  for (i in 1:ngroup.loc){
    lpsi1=matrix(lpsi[i,w],nloc,nspp,byrow=T)
    l1mpsi1=matrix(l1mpsi[i,w],nloc,nspp,byrow=T)
    lprob.exist[,i]=ltheta[i]+rowSums(dat*lpsi1+dat1m*l1mpsi1)
  }
  
  #calculation of log probability for LOCATION groups that do not exist yet
  tmp=rep(0,nloc)
  
  for (i in 1:ngroup.spp){
    cond=w==i #species assigned to group i
    soma=sum(cond)
    if (soma==0) n1ik=rep(0,nloc)
    if (soma==1) n1ik=dat[,cond]
    if (soma> 1) n1ik=rowSums(dat[,cond])
    tmp=tmp+lgamma(n1ik+1)+lgamma(soma-n1ik+1)-lgamma(2+soma)
  }
  lprob.nexist=matrix(ltheta,nloc,ngroup.loc,byrow=T)+matrix(tmp,nloc,ngroup.loc)
  
  #calculate the number of locations in each group
  tab=rep(0,ngroup.loc)
  tmp=table(z)
  tab[as.numeric(names(tmp))]=tmp
  
  #sample z
  for (i in 1:nloc){
    tab[z[i]]=tab[z[i]]-1
    lprob=rep(NA,ngroup.loc)
    cond=tab==0
    lprob[ cond]=lprob.nexist[i,cond]
    lprob[!cond]=lprob.exist[i,!cond]
    
    #get normalized probs
    tmp1=lprob-max(lprob) #for numerical stability
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
#--------------------------
#' Sample w
#' 
#' Sample the vector w containing the membership of each species
#'
#' @param lphi equal to log(phi)
#' @param dat L x S matrix, containing the presence-absence data
#' @param dat1m L x S matrix, calculated as 1-dat
#' @param lpsi equal to log(psi)
#' @param l1mpsi equal to log(1-psi)
#' @param ngroup.loc maximum number of location groups (KL)
#' @param ngroup.spp thismaximum number of species groups (KS)
#' @param nloc total number of locations (L)
#' @param nspp total number of species (S)
#' @param w vector of length S, storing the current membership of each species
#' @param z vector of length L, storing the current membership of each locations
#' @param return this function returns a vector of length S containing w 
#' @export
#' 
sample.w=function(lphi,dat,dat1m,lpsi,l1mpsi,ngroup.spp,ngroup.loc,nloc,nspp,w,z){
  #calculate log probability for groups that already exist
  lprob.exist=matrix(NA,nspp,ngroup.spp)
  for (i in 1:ngroup.spp){
    lpsi1=matrix(lpsi[z,i],nloc,nspp)
    l1mpsi1=matrix(l1mpsi[z,i],nloc,nspp)
    lprob.exist[,i]=lphi[i]+colSums(dat*lpsi1+dat1m*l1mpsi1)
  }
  
  #calculate log probability for SPECIES groups that do not exist yet
  tmp=rep(0,nspp)
  for (i in 1:ngroup.loc){
    cond=z==i
    soma=sum(cond)
    if (soma==0) n1jk=rep(0,nspp)
    if (soma==1) n1jk=dat[cond,]
    if (soma> 1) n1jk=colSums(dat[cond,])
    tmp=tmp+lgamma(n1jk+1)+lgamma(soma-n1jk+1)-lgamma(2+soma)
  }
  lprob.nexist=matrix(lphi,nspp,ngroup.spp,byrow=T)+matrix(tmp,nspp,ngroup.spp)
  
  #get the number of species assigned to each group
  tab=rep(0,ngroup.spp)
  tmp=table(w)
  tab[as.numeric(names(tmp))]=tmp
  
  #sample w
  for (i in 1:nspp){
    tab[w[i]]=tab[w[i]]-1
    lprob=rep(NA,ngroup.spp)
    cond=tab==0
    lprob[ cond]=lprob.nexist[i,cond]
    lprob[!cond]=lprob.exist[i,!cond]
    
    #get normalized probs
    tmp1=lprob-max(lprob) #for numerical stability
    tmp2=exp(tmp1) #exponentiate log probability
    prob=tmp2/sum(tmp2) #normalize to sum to 1
    
    #draw from multinomial distrib
    ind=rmultinom(1,size=1,prob=prob)
    ind1=which(ind==1)
    w[i]=ind1
    tab[ind1]=tab[ind1]+1
  }
  w
}
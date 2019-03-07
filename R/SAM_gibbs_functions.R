#' Sample from a truncated normal distribution
#' 
#' This function generates n samples from a truncated normal distribution using the CDF of the normal distribution 
#' 
#' @param n number of samples to generate 
#' @param lo lower truncation point. If there is no lower truncation, lo should be set to -Inf 
#' @param hi upper truncation point. If there is no upper truncation, hi should be set to Inf 
#' @param mu mode of the truncated normal distribution
#' @param sig dispersion parameter for truncated normal
#' @return this function returns a vector of size n
#' @export

tnorm <- function(n,lo,hi,mu,sig){   
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  
  #because qnorm can sometimes give imprecise results due to numerical issues, we need to fix this here
  cond=z<lo;    z[cond] = lo[cond]
  cond=z==-Inf; z[cond] = lo[cond]
  cond=z>hi;    z[cond] = hi[cond]
  cond=z==Inf;  z[cond] = hi[cond]
  z
}
#----------------------------------------------------------------------------------------------
#' Sample omega
#' 
#' This function generates samples of omega (the underlying latent continuous variable)
#' 
#' @param n number of samples to generate 
#' @param y L x S matrix, containing the presence-absence data 
#' @param xmat L x P matrix containing the covariates (columns) for each location (rows). 
#'             Notice that this matrix does not contain a column of 1's for the intercept            
#' @param nspp overall number of species S
#' @param nloc overall number of locations L
#' @param betas P x KS matrix containing slope parameters for each covariate (row) for each group (column)
#' @param cs group assignment for each species
#' @param alpha intercept for each species
#' @return this function returns a L x S matrix containing the underlying latent variables omega
#' @export

sample.omega=function(y,nspp,nloc,xmat,betas,cs,alpha){
  #set up empty object to hold omega
  omega=matrix(NA,nloc,nspp)
  
  #calculate overall mean
  media=xmat%*%betas
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  media1=alpha1+media[,cs]
  
  #generate omega from truncated normal distributions
  cond=y==1
  omega[ cond]=tnorm(sum(cond) ,lo=0,hi=Inf ,mu=media1[cond] ,sig=1)
  omega[!cond]=tnorm(sum(!cond),lo=-Inf,hi=0,mu=media1[!cond],sig=1)
  omega
}
#-----------------------------------------------
#' Sample alpha
#' 
#' This function generates samples alpha, the intercept for each species
#' 
#' @param nloc overall number of locations (L)
#' @param nspp overall number of species (S)
#' @param xmat this is the L x P design matrix containing the covariates (columns) for each
#'             location (rows). Notice that this matrix does not contain a column of 1's for the intercept            
#' @param betas P x KS matrix containing slope parameters for each covariate (row) and each group (column)
#' @param omega a matrix L x S containing the underlying latent variables omega
#' @param cs vector of size S containing the group assignment for each species
#' @return this function returns a vector of size S with the intercept of each species
#' @export
#' 

sample.alpha=function(nloc,xmat,betas,omega,cs,nspp){
  #calculate variance
  prec=nloc+(1/10)
  var1=1/prec
  
  #calculate part of the mean pmedia
  media=xmat%*%betas
  media1=media[,cs]
  pmedia=colSums(omega-media1)
  
  #draw alpha
  rnorm(nspp,mean=var1*pmedia,sd=sqrt(var1))
}
#-----------------------------------------------
#' Sample betas
#' 
#' This function generates samples of betas, the regression slope parameters
#' 
#' @param ngroups maximum number of species groups (KS)
#' @param cs vector of size S with group assignment for each species
#' @param nparam overall number of covariates (P)
#' @param xtx equal to t(X)*X
#' @param t.xmat equal to tranpose of the design matrix t(X)
#' @param alpha vector of size S with intercept for each species
#' @param nloc overall number of locations (L)
#' @param nspp overall number of species (S)
#' @param omega a matrix L x S containing the underlying latent variables omega
#' @return this function returns a P x KS matrix with regression slope parameters betas
#' @export
#' 

sample.betas=function(ngroups,cs,nparam,xtx,t.xmat,alpha,nloc,nspp,omega){
  #pre-calculate some stuff we will need
  i1=diag(1,nparam)
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  err=omega-alpha1
  
  betas=matrix(NA,nparam,ngroups)
  for (i in 1:ngroups){
    #calculate number of species in group i
    cond=cs==i
    nc=sum(cond)
    
    #if no species in group i, just sample from prior
    if (nc==0) betas[,i]=rmvnorm(1,mean=rep(0,nparam),sigma=i1) 
    #if there are species in group i, then sample from multivariate normal
    if (nc> 0) {
      prec=nc*xtx+i1
      var1=solve(prec)
      if (nc==1) soma.err=err[,cond]
      if (nc!=1) soma.err=rowSums(err[,cond])
      pmedia=t.xmat%*%soma.err
      betas[,i]=rmvnorm(1,mean=var1%*%pmedia,sigma=var1)
    }      
  }
  betas
}
#-----------------------------------------------
#' Sample cs
#' 
#' This function generates samples of the group assignment variable cs
#' 
#' @param ngroups maximum number of species groups (KS)
#' @param omega a matrix L x S containing the underlying latent variables omega
#' @param xmat the design matrix of dimension L x P 
#' @param alpha vector of size S with intercept for each species
#' @param betas P x KS matrix containing slope parameters for each covariate (row) and each group (column)
#' @param theta vector of size KS with the probability of each species group
#' @param nloc overall number of locations (L)
#' @param nspp overall number of species (S)
#' @param invSigma.precalc precision matrix for new group
#' @param Sigma.precalc covariance matrix for new group
#' @param lds log of determinant of covariance matrix for new group
#' @param cs vector of size S with group assignment for each species
#' @param nparam overall number of covariates (P)
#' @return a vector of size S with the group assignment variable cs
#' @export

sample.cs=function(ngroups,omega,xmat,alpha,betas,theta,nspp,nloc,
                   InvSigma.precalc,Sigma.precalc,lds,cs,nparam){
  #pre-calculate stuff for new group when beta is integrated out
  p1=nloc*log(2*pi)
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  err=omega-alpha1
  p2s=colSums(err^2)
  q1s=t(xmat)%*%err  
  mu=Sigma.precalc%*%q1s
  p3s=-t(mu)%*%InvSigma.precalc%*%mu; 
  p3s=as.numeric(diag(p3s))
  p4=-lds
  ltheta=log(theta)
  p5k=-2*ltheta
  
  #determine the number of species assigned to each group
  tab=rep(0,ngroups)
  tmp=table(cs)
  tab[as.numeric(names(tmp))]=tmp

  lprob.exist=matrix(NA,nspp,ngroups)
  lprob.not.exist=matrix(NA,nspp,ngroups)
  k=log(2*pi)
  resto.media=xmat%*%betas #mean for each location based on betas for each group 
  for (i in 1:ngroups){
    #calculate log probability if group already exists
    resto.media1=matrix(resto.media[,i],nloc,nspp)
    err2=(omega-(alpha1+resto.media1))^2
    dnorm1=-0.5*(k+err2)
    lprob.exist[,i]=colSums(dnorm1)+ltheta[i]
    #check results: z=colSums(dnorm(omega,mean=alpha1+resto.media1,sd=1,log=T))
    
    #calculate log probability for a group that does not yet exist
    lprob.not.exist[,i]=-0.5*(p1+p2s+p3s+p4+p5k[i])
  }

  #for each species, sample the cs's
  for (i in 1:nspp){
    tab[cs[i]]=tab[cs[i]]-1 #all but this species
    
    #get probabilities for each group
    lprob=rep(NA,ngroups)
    cond=tab==0
    lprob[ cond]=lprob.not.exist[i,cond] #groups that don't exist yet (i.e., empty groups)
    lprob[!cond]=lprob.exist[i,!cond] #groups that already exist 

    #normalize probabilities
    max1=max(lprob)
    tmp=lprob-max1
    tmp1=exp(tmp)
    prob=tmp1/sum(tmp1)
    
    #sample cs
    ind=rmultinom(1,size=1,prob)
    ind1=which(ind==1)
    cs[i]=ind1
    tab[ind1]=tab[ind1]+1
  }
  cs
}
#-----------------------------------------------
#' Sample theta
#' 
#' This function generates a sample of the vector of probabilities associated with each species group
#' 
#' @param cs vector of size S with group assignment for each species
#' @param ngroups maximum number of species groups (KS)
#' @param gamma the parameter for the stick-breaking prior
#' @param burnin number of iterations to discard as part of burn-in phase
#' @param gibbs.step current iteration of the gibbs sampler
#' @param betas P x KS matrix containing slope parameters for each covariate (row) and each group (column)
#' @param theta vector of size KS with the probability of each species group
#' @return a list of 4 items (theta,betas,cs,v)
#' @export

sample.theta.sam=function(cs,ngroups,gamma,burnin,gibbs.step,betas,theta){
  #re-order thetas in decreasing order. 
  #Based on this re-ordering, re-order betas and cs
  if(gibbs.step<burnin & gibbs.step%%50==0){
    ind=order(theta,decreasing=T)
    theta=theta[ind]
    betas=betas[,ind]
    
    #get cs.new
    cs.new=cs; cs.new[]=NA
    for (i in 1:ngroups){
      cond=cs==ind[i]
      cs.new[cond]=i
    }
    cs=cs.new
  }
  
  #count the number of species assigned to each group
  n=rep(0,ngroups)
  tmp=table(cs)
  n[as.numeric(names(tmp))]=tmp
  
  #sample v and calculate the corresponding theta
  v=theta=rep(NA,ngroups)
  prod=1
  for (i in 1:(ngroups-1)){
    n.greater.k=n[-(1:i)]
    v[i]=rbeta(1,n[i]+1,sum(n.greater.k)+gamma)
    theta[i]=v[i]*prod
    prod=prod*(1-v[i])
  }
  theta[ngroups]=prod

  #to avoid numerical issues
  cond=v>0.99999
  v[cond]=0.99999
  
  #output results
  list(theta=theta,betas=betas,cs=cs,v=v)
}
#----------------------------
#' Calculate the log-likelihood
#' 
#' This function calculates the log-likelihood based on the current parameter estimates
#' 
#' @param y L x S matrix, with the presence-absence data 
#' @param omega L x S matrix with the underlying latent variables omega
#' @param nloc overall number of locations (L)
#' @param nspp overall number of species (S)
#' @param xmat L x P design matrix containing the covariates (columns) for each location (rows). 
#'             Notice that this matrix does not contain a column of 1's for the intercept            
#' @param betas P x KS matrix containing slope parameters for each covariate (row) and each group (column)
#' @param cs vector of length S wth group assignment for each species
#' @param alpha vector of length S with intercept for each species
#' @return log-likelihood of probit regression models
#' @export

get.logl=function(y,omega,nspp,nloc,xmat,betas,cs,alpha){
  #get overall mean
  media=xmat%*%betas
  alpha1=matrix(alpha,nloc,nspp,byrow=T)
  media1=alpha1+media[,cs]
  
  #calculate log-likelihood
  p=y*pnorm(media1,log=T)+(1-y)*pnorm(media1,lower.tail=F,log.p=T) 
  #numerically, log(1-pnorm(media1)) is worse than pnorm(media1,lower.tail=F,log.p=T) 
  sum(p)
}
#----------------------------
#' Sample gamma
#' 
#' This function generates a sample of the truncated stick breaking prior parameter gamma
#' 
#' @param v vector of size KS with the stick-breaking parameters 
#' @param ngroups maximum number of species groups (KS)
#' @param gamma.possib set of possible values for gamma
#' @return a sample of gamma
#' @export

sample.gamma.sam=function(v,ngroups,gamma.possib){
  #calculate the log probability for each possible gamma value
  ngamma=length(gamma.possib)
  soma=sum(log(1-v[-ngroups]))
  k=(ngroups-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma
  #to check code: sum(dbeta(v[-ngroups],1,gamma.possib[5],log=T))
  
  #exponentiate and normalize probabilities
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)
  
  #sample from categorical distribution
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}
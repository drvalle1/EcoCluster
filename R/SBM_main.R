#' @useDynLib tsbp
#' @importFrom Rcpp sourceCpp
NULL

#' Main function of the Stochastic Block Model (SBM)
#' 
#' Runs the Gibbs sampler and returns samples from the posterior distribution
#' 
#' @param dat matrix has L rows (e.g., locations) and S columns (e.g., species),
#'            containing the presence-absence data 
#' @param ngroup.loc maximum number of location groups (KL)
#' @param ngroup.spp maximum number of species groups (KS)
#' @param ngibbs  number of Gibbs sampler iterations   
#' @param burnin  number of iterations to discard as part of burn-in phase
#' @param return this function returns samples from the posterior distribution in the form of a 
#'               list containing several matrices. These matrices have ngibbs-burnin rows and 
#'               contain the posterior distribution for:
#'               - theta: probability of each location group
#'               - phi:   probability of each species group
#'               - llk:   log-likelihood for each iteration
#'               - psi:   presence probability for each location group and species group                                
#'               - z:     cluster assignment of each location
#'               - w:     cluster assignment of each species
#'               - gamma: TSB prior parameters (one for location groups and the other for species groups)
#' @export

SBM=function(dat,ngroup.loc,ngroup.spp,ngibbs,burnin){
  #basic settings
  nloc=nrow(dat)
  nspp=ncol(dat)
  
  #get initial values for parameters
  theta=rep(1/ngroup.loc,ngroup.loc)
  phi=rep(1/ngroup.spp,ngroup.spp)
  z=sample(1:ngroup.loc,size=nloc,replace=T)
  w=sample(1:ngroup.spp,size=nspp,replace=T)
  tmp=runif(ngroup.loc*ngroup.spp)
  psi=matrix(tmp,ngroup.loc,ngroup.spp)
  gamma.v=0.1
  gamma.u=0.1
  gamma.possib=seq(from=0.1,to=1,by=0.05) #discretization of possible gamma values
  
  #pre-calculate useful quantities
  dat1m=1-dat
  
  #vectors to store results from iterations of gibbs sampler
  vec.llk=rep(NA,ngibbs)
  vec.theta=matrix(NA,ngibbs,ngroup.loc)
  vec.phi=matrix(NA,ngibbs,ngroup.spp)
  vec.psi=matrix(NA,ngibbs,ngroup.spp*ngroup.loc)
  vec.z=matrix(NA,ngibbs,nloc)
  vec.w=matrix(NA,ngibbs,nspp)
  vec.gamma=matrix(NA,ngibbs,2); colnames(vec.gamma)=c('gamma.u','gamma.v')
  
  #start gibbs sampler
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)
    
    #sample psi
    psi=sample.psi(z=z,w=w,dat=dat,ngroup.loc=ngroup.loc,ngroup.spp=ngroup.spp)
    
    #sample theta, vk (in some cases, this will lead to changes in z and psi as well)
    tmp=sample.theta(ngroup.loc=ngroup.loc,gamma.v=gamma.v,burnin=burnin,
                     gibbs.step=i,theta=theta,psi=psi,z=z)
    theta=tmp$theta
    vk=tmp$vk
    z=tmp$z
    psi=tmp$psi
    
    #sample phi, uk (in some cases, this will lead to changes in w and psi as well)
    tmp=sample.phi(ngroup.spp=ngroup.spp,gamma.u=gamma.u,burnin=burnin,
                   gibbs.step=i,phi=phi,psi=psi,w=w)
    phi=tmp$phi
    uk=tmp$uk
    w=tmp$w
    psi=tmp$psi
    
    #pre-calculate useful quantities
    lpsi=log(psi)
    l1mpsi=log(1-psi)
    
    #sample z
    z=sample.z(ltheta=log(theta),dat=dat,dat1m=dat1m,lpsi=lpsi,l1mpsi=l1mpsi,
               ngroup.loc=ngroup.loc,ngroup.spp=ngroup.spp,nloc=nloc,nspp=nspp,
               w=w,z=z)

    #sample w
    w=sample.w(lphi=log(phi),dat=dat,dat1m=dat1m,lpsi=lpsi,l1mpsi=l1mpsi,
               ngroup.spp=ngroup.spp,ngroup.loc=ngroup.loc,nloc=nloc,nspp=nspp,
               w=w,z=z)
    
    #sample gammas
    gamma.u=sample.gamma.u(uk=uk,gamma.possib=gamma.possib,ngroup.spp=ngroup.spp)
    gamma.v=sample.gamma.v(vk=vk,gamma.possib=gamma.possib,ngroup.loc=ngroup.loc)
    
    #get log-likelihood
    psi1=psi[z,]
    psi2=psi1[,w]
    loglikel=dat*log(psi2)+dat1m*log(1-psi2)
    
    #store Gibbs sampler results
    vec.llk[i]=sum(loglikel)    
    vec.theta[i,]=theta
    vec.phi[i,]=phi
    vec.psi[i,]=psi
    vec.z[i,]=z
    vec.w[i,]=w
    vec.gamma[i,]=c(gamma.u,gamma.v)
  }
  
  #output a list containing several matrices after removing burn-in
  seq1=burnin:ngibbs
  list(theta=vec.theta[seq1,],
       phi=vec.phi[seq1,],
       llk=vec.llk[seq1],
       psi=vec.psi[seq1,],
       z=vec.z[seq1,],
       w=vec.w[seq1,],
       gamma=vec.gamma[seq1,])
}

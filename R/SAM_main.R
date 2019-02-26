#' @useDynLib tsbp
#' @importFrom Rcpp sourceCpp
#' @importFrom mvtnorm rmvnorm
NULL

#' Main function of the Species Archetype Model (SAM)
#' 
#' Runs the Gibbs sampler and returns samples from the posterior distribution
#' 
#' @param y L X S matrix, containing the presence-absence data 
#' @param xmat L X P design matrix containing the covariates (columns) for each
#'             location (rows). Notice that this matrix does not contain a column of 1's for the intercept            
#' @param ngroups maximum number of species groups (KS)
#' @param ngibbs  number of Gibbs sampler iterations          
#' @param burnin  number of iterations to discard as part of burn-in phase
#' @param return this function returns a list containing several matrices, all of which have 
#'               ngibbs-burnin rows, containing samples from the posterior distribution of:
#'               - theta: probability associated with each species group
#'               - logl:  log-likelihood for each iteration
#'               - betas: slope parameters for each group
#'               - cs:    cluster assignment of each species
#'               - alpha: intercept of each species
#'               - gamma: TSB prior parameter
#' @export

gibbs.SAM=function(y,xmat,ngroups,ngibbs,burnin){
#pre-calculate useful quantities
t.xmat=t(xmat)
xtx=t.xmat%*%xmat
nloc=nrow(y) #number of locations
nspp=ncol(y) #number of species
nparam=ncol(xmat) #number of slope parameters
gamma.possib=seq(from=0.1,to=1,by=0.05) #possible gamma values
InvSigma.precalc=xtx+diag(1,nparam) #precision matrix when beta is integrated out for new group
Sigma.precalc=solve(InvSigma.precalc) #covariance matrix when beta is integrated out for new group
lds=determinant(Sigma.precalc)$modulus[1] #log of determinant of covariance matrix when beta is integrated out for new group

#initial parameter values
betas=matrix(0,nparam,ngroups)
alpha=rep(0,nspp)
omega=matrix(-1,nloc,nspp)
omega[y==1]=1
cs=sample(1:ngroups,size=nspp,replace=T)
theta=rep(1/ngroups,ngroups)
gamma=0.1

#matrices to store MCMC results
vec.betas=matrix(NA,ngibbs,nparam*ngroups)
vec.alpha=matrix(NA,ngibbs,nspp)
vec.theta=matrix(NA,ngibbs,ngroups)
vec.logl=matrix(NA,ngibbs,1)
vec.gamma=rep(NA,ngibbs,1)
vec.cs=matrix(NA,ngibbs,nspp)

#main Gibbs sampler loop
for (i in 1:ngibbs){
  print(i)
  
  #sample species assignment variable cs
  cs=sample.cs(ngroups=ngroups,omega=omega,xmat=xmat,alpha=alpha,betas=betas,theta=theta,
               nspp=nspp,nloc=nloc,InvSigma.precalc=InvSigma.precalc,
               Sigma.precalc=Sigma.precalc,lds=lds,cs=cs,nparam=nparam)

  #sample regression slope parameters betas
  betas=sample.betas(ngroups=ngroups,cs=cs,nparam=nparam,xtx=xtx,t.xmat=t.xmat,alpha=alpha,
                     nloc=nloc,nspp=nspp,omega=omega)

  #sample omega (underlying latent normally distributed variables)
  omega=sample.omega(y=y,nspp=nspp,nloc=nloc,xmat=xmat,betas=betas,cs=cs,alpha=alpha)

  #sample species specific intercepts alpha
  alpha=sample.alpha(nloc=nloc,xmat=xmat,betas=betas,omega=omega,cs=cs,nspp=nspp)
  
  #sample v and the corresponding probability of each species group theta
  #betas and cs might change order due to this function
  tmp=sample.theta.sam(cs=cs,ngroups=ngroups,gamma=gamma,burnin=burnin,gibbs.step=i,betas=betas,theta=theta)
  betas=tmp$betas
  cs=tmp$cs
  theta=tmp$theta
  v=tmp$v[-ngroups]
  
  #sample the TSB prior parameter gamma
  gamma=sample.gamma.sam(v=v,ngroups=ngroups,gamma.possib=gamma.possib)
  
  #get loglikelihood
  logl=get.logl(y=y,omega=omega,nspp=nspp,nloc=nloc,xmat=xmat,betas=betas,cs=cs,alpha=alpha)
  
  #store results
  vec.betas[i,]=betas
  vec.alpha[i,]=alpha
  vec.theta[i,]=theta
  vec.logl[i]=logl
  vec.cs[i,]=cs
  vec.gamma[i]=gamma
}
  #output results after eliminating iterations from burn-in period
  seq1=burnin:ngibbs
  list(theta=vec.theta[seq1,],
       logl=vec.logl[seq1],
       betas=vec.betas[seq1,],
       alpha=vec.alpha[seq1,],
       cs=vec.cs[seq1,],
       gamma=vec.gamma[seq1])
}
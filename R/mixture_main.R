#' @useDynLib EcoCluster
#' @importFrom Rcpp sourceCpp
NULL

#' Runs the Gibbs sampler and returns samples from the posterior distribution
#' 
#' @param dat this matrix has L rows (locations) and S columns (species)
#'            and contains the presence-absence data (i.e., number of times 
#'            a given species was observed at a given location)
#' @param nl  this vector has L elements (locations) and contains the number of observation 
#'            opportunities at each location
#' @param ngroup maximum number of location groups (K)
#' @param ngibbs number of Gibbs sampler iterations 
#' @param burnin number of iterations to discard as burn-in  
#' @param a.prior 'a' parameter for prior beta distribution
#' @param b.prior 'b' parameter for prior beta distribution                      
#' @return this function returns a list containing several matrices.
#'         These matrices have ngibbs-burnin rows and contain samples from the posterior distribution for:
#'         \itemize{
#'            \item phi:   probability of observing each species in each group
#'            \item theta: probability of each location group
#'            \item logl:  log-likelihood
#'            \item z:     cluster assignment of each location
#'            \item gamma: TSB prior parameter
#'         }
#' @export

mixture.gibbs=function(dat,ngroup,nl,ngibbs,burnin,a.prior,b.prior){
    
    #useful pre-calculated quantities 
    nloc=nrow(dat)
    nspp=ncol(dat)
    nl.mat=matrix(nl,nloc,nspp)
    n.minus.y=nl.mat-dat
    constant=nspp*(lgamma(a.prior+b.prior)-lgamma(a.prior)-lgamma(b.prior))
    
    #initial parameter values
    z=sample(1:ngroup,size=nloc,replace=T)
    theta=rep(1/ngroup,ngroup)
    tmp=runif(ngroup*nspp)
    phi=matrix(tmp,ngroup,nspp)
    gamma1=0.1
    gamma.possib=seq(from=0.1,to=1,by=0.05) #possible values for gamma
    
    #to store results from gibbs sampler
    store.phi=matrix(NA,ngibbs,nspp*ngroup)
    store.theta=matrix(NA,ngibbs,ngroup)
    store.z=matrix(NA,ngibbs,nloc)
    store.gamma=matrix(NA,ngibbs,1)
    store.logl=rep(NA,ngibbs)
    
    #run gibbs sampler
    norder=50
    for (i in 1:ngibbs){
        print(i)
        
        #sample group allocation vector z
        z=update.z(dat=dat,nl=nl,n.minus.y=n.minus.y,phi=phi,theta=theta,
                   ngroup=ngroup,nloc=nloc,nspp=nspp,z=z,
                   a.prior=a.prior,b.prior=b.prior,constant=constant)
        
        #re-order groups if necessary
        if (i%%norder==0 & i<burnin){
            order1=order(theta,decreasing=T)
            theta=theta[order1]
            phi=phi[order1,]
            znew=rep(NA,nloc)
            for (j in 1:ngroup){
                cond=z==order1[j]
                znew[cond]=j
            }
            z=znew
        }
        
        #summarize the data
        #determine the number of times a particular species was observed in each location group. This is ncs1
        #determine the number of times a particular species was not observed in each location group. This is ncs0
        tmp=ncs(dat=dat,nminusy=n.minus.y,z=z-1,nspp=nspp,nloc=nloc,ngroup=ngroup)
        
        #sample the phi matrix containing the spp composition characterization of each group of locations
        phi=matrix(rbeta(ngroup*nspp,tmp$ncs1+a.prior,tmp$ncs0+b.prior),ngroup,nspp)
        
        #sample v and theta (this might also result in changes for z and phi)
        tmp=update.theta(z=z,ngroup=ngroup,gamma1=gamma1,burnin=burnin,gibbs.step=i,theta=theta,phi=phi)
        theta=tmp$theta
        z=tmp$z
        v=tmp$v
        phi=tmp$phi
        
        #sample the TSB parameter gamma
        gamma1=sample.gamma(v=v,ngroup=ngroup,gamma.possib=gamma.possib)
        
        #get loglikelihood
        prob=phi[z,]
        logl=dbinom(dat,size=nl.mat,prob=prob,log=T)
        
        #store results
        store.phi[i,]=phi
        store.theta[i,]=theta
        store.logl[i]=sum(logl)
        store.z[i,]=z
        store.gamma[i]=gamma1
    }
    
    #output MCMC results after discarding the burn-in phase
    seq1=burnin:ngibbs
    list(phi=store.phi[seq1,],theta=store.theta[seq1,],logl=store.logl[seq1],
         z=store.z[seq1,],gamma=store.gamma[seq1])
}
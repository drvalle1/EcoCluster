## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----eval=TRUE-----------------------------------------------------------
rm(list=ls(all=TRUE))
set.seed(1)

#general settings
nloc=300 #number of locations
nspp=50  #number of species
ngroup=5 #number of groups

#set parameters
z=sample(1:ngroup,size=nloc,replace=T)         #cluster assignment of each location
phi=matrix(rbeta(ngroup*nspp,1,1),ngroup,nspp) #species composition of each cluster  

#generate data
y=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  phi1=phi[z[i],]
  y[i,]=rbinom(nspp,size=1,prob=phi1)
}
colnames(y)=paste0('spp',1:nspp)
rownames(y)=paste0('loc',1:nloc)

head(y[,1:13])

## ----eval=TRUE,results='hide'--------------------------------------------
library('EcoCluster')
ngibbs=1000
res=mixture.gibbs.main.func(dat=y,ngroup=50,ngibbs=ngibbs,burnin=ngibbs/2)

## ----eval=TRUE-----------------------------------------------------------
str(res)

## ---- eval=TRUE, fig.align='center'--------------------------------------
plot(res$logl,type='l',main='convergence assessment',xlab='Gibbs iterations',ylab='log-likelihood')

## ---- eval=TRUE, fig.align='center'--------------------------------------
theta=colMeans(res$theta)
plot(theta[1:20],type='h',xlab='clusters',main='number of groups',ylab='theta',lwd=2)
abline(v=5,col='red',lty=3,lwd=0.1)

## ----eval=TRUE-----------------------------------------------------------
rm(list=ls(all=TRUE))
set.seed(1)

#general settings
nloc=300      #number of locations
nspp=50       #number of species
ngroup.loc=5  #number of location groups
ngroup.spp=3  #number of species groups

#set parameters
theta=rep(1/ngroup.loc,ngroup.loc)  #probability of each location group
phi=rep(1/ngroup.spp,ngroup.spp)    #probability of each species group

#probabilities associated with each location and species group
psi=matrix(c(0.05,0.5,0.95,
             0.5,0.05,0.95,
             0.05,0.95,0.5,
             0.5,0.95,0.05,
             0.1,0.5,0.05),ngroup.loc,ngroup.spp,byrow=T)

#get location group and species group assignments
z=sample(1:ngroup.loc,size=nloc,replace=T); 
w=sample(1:ngroup.spp,size=nspp,replace=T); 

#generate data
y=matrix(NA,nloc,nspp)
for (i in 1:nloc){
  for (j in 1:nspp){
    y[i,j]=rbinom(1,size=1,prob=psi[z[i],w[j]])    
  }
}
colnames(y)=paste0('spp',1:nspp)
rownames(y)=paste0('loc',1:nloc)
head(y[,1:13])

## ----eval=TRUE,results='hide'--------------------------------------------
library('EcoCluster')
ngibbs=1000
res=SBM(dat=y,ngroup.loc=50,ngroup.spp=50,ngibbs=ngibbs,burnin=ngibbs/2)

## ----eval=TRUE-----------------------------------------------------------
str(res)

## ----eval=TRUE, fig.align='center'---------------------------------------
plot(res$llk,type='l',main='convergence assessment',xlab='Gibbs iterations',ylab='log-likelihood')

## ----eval=TRUE, fig.align='center'---------------------------------------
theta=colMeans(res$theta)
plot(theta[1:20],type='h',xlab='clusters',main='number of location groups',ylab='theta',lwd=2)
abline(v=5,col='red',lty=3,lwd=0.1)

phi=colMeans(res$phi)
plot(phi[1:20],type='h',xlab='clusters',main='number of species groups',ylab='phi',lwd=2)
abline(v=3,col='red',lty=3,lwd=0.1)

## ----eval=TRUE-----------------------------------------------------------
rm(list=ls(all=TRUE))
library('mvtnorm')
set.seed(1)

#general settings
nloc=300   #number of locations
nspp=50    #number of species
nparam=6   #number of covariates
ngroup1=5  #number of species groups

#create design matrix with covariates
xmat=matrix(rnorm(nparam*nloc),nloc,nparam)

#set parameters
alpha=rnorm(nspp,mean=0,sd=0.4) #intercept of each species
tmp=rnorm(nparam*ngroup1)
betas=matrix(tmp,nparam,ngroup1) #slope parameters for each group
cs=sample(1:ngroup1,size=nspp,replace=T) #cluster assignment for each species

#generate data assuming a probit formulation    
omega=matrix(NA,nloc,nspp)
for (i in 1:nspp){
  media=alpha[i]+xmat%*%betas[,cs[i]]
  omega[,i]=rnorm(nloc,mean=media,sd=1)
}
y=omega.true=omega
y[omega>0]=1
y[omega<0]=0

colnames(y)=paste0('spp',1:nspp)
rownames(y)=paste0('loc',1:nloc)
head(y[,1:13])

## ---- eval=TRUE, results='hide'------------------------------------------
library('EcoCluster')
ngibbs=1000

res=gibbs.SAM(y=y,xmat=xmat,ngroups=50,ngibbs=ngibbs,burnin=ngibbs/2)

## ----eval=TRUE-----------------------------------------------------------
str(res)

## ---- eval=TRUE, fig.align='center'--------------------------------------
plot(res$logl,type='l',main='convergence assessment',xlab='Gibbs iterations',ylab='log-likelihood')

## ---- eval=TRUE, fig.align='center'--------------------------------------
theta=colMeans(res$theta)
plot(theta[1:20],type='h',xlab='clusters',main='number of species groups',ylab='theta',lwd=2)
abline(v=5,col='red',lty=3,lwd=0.1)


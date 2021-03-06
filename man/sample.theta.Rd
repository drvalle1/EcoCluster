% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SBM_gibbs_functions.R
\name{sample.theta}
\alias{sample.theta}
\title{Sample theta and vk parameters}
\usage{
sample.theta(ngroup.loc, gamma.v, burnin, gibbs.step, theta, psi, z)
}
\arguments{
\item{ngroup.loc}{maximum number of groups for locations (KL)}

\item{gamma.v}{truncated stick-breaking prior parameter for the
location groups. This value should be between 0 and 1, and
small values enforce more parsimonius results (i.e., fewer groups)}

\item{burnin}{number of MCMC samples that are going to be thrown out as
part of the burn-in phrase}

\item{gibbs.step}{current iteration of the gibbs sampler}

\item{theta}{vector of size KL containing the current estimate of theta (i.e., probability of each location group)}

\item{psi}{KL x KS matrix  containing the current estimate of psi}

\item{z}{this is a vector of length L, storing the current membership of each location}
}
\value{
this function returns a list containing 4 items (theta, vk, psi, and z)
}
\description{
Sample vk parameters from their full conditional distributions and convert
these vk parameters into theta parameters
}

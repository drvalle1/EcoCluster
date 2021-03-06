% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAM_gibbs_functions.R
\name{sample.alpha}
\alias{sample.alpha}
\title{Sample alpha}
\usage{
sample.alpha(nloc, xmat, betas, omega, cs, nspp)
}
\arguments{
\item{nloc}{overall number of locations (L)}

\item{xmat}{this is the L x P design matrix containing the covariates (columns) for each
location (rows). Notice that this matrix does not contain a column of 1's for the intercept}

\item{betas}{P x KS matrix containing slope parameters for each covariate (row) and each group (column)}

\item{omega}{a matrix L x S containing the underlying latent variables omega}

\item{cs}{vector of size S containing the group assignment for each species}

\item{nspp}{overall number of species (S)}
}
\value{
this function returns a vector of size S with the intercept of each species
}
\description{
This function generates samples alpha, the intercept for each species
}

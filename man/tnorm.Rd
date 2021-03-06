% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAM_gibbs_functions.R
\name{tnorm}
\alias{tnorm}
\title{Sample from a truncated normal distribution}
\usage{
tnorm(n, lo, hi, mu, sig)
}
\arguments{
\item{n}{number of samples to generate}

\item{lo}{lower truncation point. If there is no lower truncation, lo should be set to -Inf}

\item{hi}{upper truncation point. If there is no upper truncation, hi should be set to Inf}

\item{mu}{mode of the truncated normal distribution}

\item{sig}{dispersion parameter for truncated normal}
}
\value{
this function returns a vector of size n
}
\description{
This function generates n samples from a truncated normal distribution using the CDF of the normal distribution
}

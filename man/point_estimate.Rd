% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract.R
\name{point_estimate}
\alias{point_estimate}
\title{Extract point estimates of parameters from a fit object}
\usage{
point_estimate(fit, pars = c("pi", "theta", "z"), ...)
}
\arguments{
\item{fit}{A rater fit object}

\item{pars}{A character vector of parameter names to return. By default
\code{c("pi", "theta", "z")}.}

\item{...}{Extra arguments}
}
\value{
A named list of the parameter estimates.
}
\description{
Extract point estimates of parameters from a fit object
}
\details{
If the passed fit object was fit using MCMC then the posterior
means are returned. If it was fit through optimisation the maximum a
priori (MAP) estimates are returned. The z parameter returned is the
value of class probabilities which is largest. To return the full
posterior distributions of the latent class use \code{class_probabilities()}.
}
\examples{

\donttest{
# A model fit using MCMC.
mcmc_fit <- rater(anesthesia, "dawid_skene")

# This will return the posterior mean (except for z)
post_mean_estimate <- point_estimate(mcmc_fit)

# A model fit using optimisation.
optim_fit <- rater(anesthesia, dawid_skene(), method = "optim")

# This will output MAP estimates of the parameters.
map_estimate <- point_estimate(optim_fit)

}

}
\seealso{
\code{class_probabilities()}
}

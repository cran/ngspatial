\name{summary.autologistic}
\alias{summary.autologistic}
\title{Print a summary of a centered autologistic model fit.}
\usage{
  \method{summary}{autologistic} (object, alpha = 0.05,
    digits = 4, ...)
}
\arguments{
  \item{object}{an object of class \code{autologistic},
  typically the result of a call to
  \code{\link{autologistic}}.}

  \item{alpha}{the significance level for the quantile/HPD
  intervals. The default is 0.05.}

  \item{digits}{the number of significant digits to
  display. The default is 4.}

  \item{\dots}{additional arguments.}
}
\description{
  Print a summary of a centered autologistic model fit.
}
\details{
  This function displays (1) the call to
  \code{\link{autologistic}}, (2) a table of estimates, and
  (3) the size of the bootstrap/posterior samples. Each row
  of the table of estimates shows the estimated regression
  coefficient, the \eqn{(\alpha/2)100\%}{(\alpha/2)100\%}
  and \eqn{(1-\alpha/2)100\%}{(1-\alpha/2)100\%} bootstrap
  quantiles for the MPLE coefficient or
  \eqn{(1-\alpha)100\%}{(1-\alpha)100\%} HPD interval for
  the Bayesian coefficient, and the Monte Carlo standard
  error.
}
\seealso{
  \code{\link{autologistic}}
}


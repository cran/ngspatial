\name{summary.sparse.sglmm}
\alias{summary.sparse.sglmm}
\title{Print a summary of a sparse SGLMM fit.}
\usage{
  \method{summary}{sparse.sglmm} (object, alpha = 0.05,
    digits = 4, ...)
}
\arguments{
  \item{object}{an object of class \code{sparse.sglmm},
  typically the result of a call to
  \code{\link{sparse.sglmm}}.}

  \item{alpha}{the significance level used to compute the
  HPD intervals. The default is 0.05.}

  \item{digits}{the number of significant digits to
  display. The default is 4.}

  \item{\dots}{additional arguments.}
}
\description{
  Print a summary of a sparse SGLMM fit.
}
\details{
  This function displays (1) the call to
  \code{\link{sparse.sglmm}}, (2) a table of estimates, (3)
  the DIC value for the fit, and (4) the size of the
  posterior sample. Each row of the table of estimates
  shows the estimated regression coefficient, the
  \eqn{(1-\alpha)100\%}{(1-\alpha)100\%} HPD interval for
  the coefficient, and the Monte Carlo standard error.
}
\seealso{
  \code{\link{sparse.sglmm}}
}


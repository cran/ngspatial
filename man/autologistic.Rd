\name{autologistic}
\alias{autologistic}
\title{Fit a centered autologistic model using maximum pseudolikelihood estimation or MCMC for Bayesian inference.}
\usage{
  autologistic(formula, data, A, method = c("pl", "bayes"),
    optit = 1000, model = TRUE, x = FALSE, y = FALSE,
    type = c("SOCK", "PVM", "MPI", "NWS"), bootit = 1000,
    parallel = TRUE, nodes, trainit = 1e+05, tol = 0.01,
    minit = 10000, maxit = 1e+06, sigma = 1e+06,
    eta.max = 2)
}
\arguments{
  \item{formula}{an object of class \code{\link{formula}}:
  a symbolic description of the model to be fitted.}

  \item{data}{an optional data frame, list, or environment
  (or object coercible by \code{\link{as.data.frame}} to a
  data frame) containing the variables in the model. If not
  found in \code{data}, the variables are taken from
  \code{environment(formula)}, typically the environment
  from which \code{autologistic} is called.}

  \item{A}{the adjacency matrix for the underlying graph,
  which is assumed to be undirected and free of loops and
  parallel edges.}

  \item{method}{the method to use for inference.
  \dQuote{\code{pl}} (the default) enables maximum
  pseudolikelihood estimation, and \dQuote{\code{bayes}}
  enables Bayesian inference.}

  \item{optit}{the maximum number of iterations to be used
  by \code{\link{optim}} in obtaining the MPLE estimate of
  \eqn{\theta}. Defaults to 1,000.}

  \item{parallel}{(PL) a boolean variable indicating
  whether to use parallel bootstrapping, which requires the
  \pkg{snow} package. Defaults to \code{TRUE}, in which
  case the number of nodes must be supplied.}

  \item{nodes}{(PL) the number of nodes to use for parallel
  bootstrapping.}

  \item{type}{(PL) the type of cluster to use for parallel
  bootstrapping. The available types are
  \dQuote{\code{SOCK}}, \dQuote{\code{PVM}},
  \dQuote{\code{MPI}}, and \dQuote{\code{NWS}}. The default
  type is \dQuote{\code{SOCK}}.}

  \item{bootit}{(PL) the desired size of the bootstrap
  sample. Defaults to 1,000.}

  \item{trainit}{(Bayes) the number of iterations to use
  for estimating the posterior covariance matrix. Defaults
  to 100,000.}

  \item{tol}{(Bayes) a tolerance. If all Monte Carlo
  standard errors are smaller than \code{tol}, no more
  samples are drawn from the posterior. Defaults to 0.01.}

  \item{minit}{(Bayes) the minimum sample size. This should
  be large enough to permit accurate estimation of Monte
  Carlo standard errors. Defaults to 10,000.}

  \item{maxit}{(Bayes) the maximum sample size. Sampling
  from the posterior terminates when all Monte Carlo
  standard errors are smaller than \code{tol} or when
  \code{maxit} samples have been drawn, whichever happens
  first. Defaults to 1,000,000.}

  \item{sigma}{(Bayes) a scalar or a \eqn{(p-1)}-vector
  providing the variance(s) of the spherical normal prior
  for \eqn{\beta}. Defaults to 1,000,000.}

  \item{eta.max}{(Bayes) the upper limit for \eqn{\eta}.
  Defaults to 2. The lower limit is 0.}

  \item{model}{a logical value indicating whether the model
  frame should be included as a component of the returned
  value.}

  \item{x}{a logical value indicating whether the model
  matrix used in the fitting process should be returned as
  a component of the returned value.}

  \item{y}{a logical value indicating whether the response
  vector used in the fitting process should be returned as
  a component of the returned value.}
}
\value{
  \code{autologistic} returns an object of class
  \dQuote{\code{autologistic}}, which is a list containing
  the following components.  \item{coefficients}{the point
  estimate of \eqn{\theta}.} \item{fitted.values}{the
  fitted mean values, obtained by transforming the linear
  predictors by the inverse of the link function.}
  \item{linear.predictors}{the linear fit on link scale.}
  \item{residuals}{the response residuals.} \item{iter}{the
  size of the bootstrap/posterior sample.} \item{sample}{an
  \code{iter} by \eqn{p} matrix containing the
  bootstrap/posterior samples.} \item{mcse}{a
  \eqn{p}-vector of Monte Carlo standard errors.}
  \item{V}{(Bayes) the estimated posterior covariance
  matrix from the training run.} \item{accept}{(Bayes) the
  acceptance rate for the MCMC sampler.} \item{y}{if
  requested (the default), the \code{y} vector used.}
  \item{X}{if requested, the model matrix.} \item{model}{if
  requested (the default), the model frame.}
  \item{call}{the matched call.} \item{formula}{the formula
  supplied.} \item{method}{the method used for inference.}
  \item{convergence}{an integer code. The code has value 0
  if \code{\link{optim}} succeeded in optimizing the
  pseudolikelihood. Possible error codes are 1 and 10. The
  former indicates that the iteration limit was reached
  before optimization completed. The latter indicates that
  the Nelder-Mead simplex degenerated.} \item{message}{a
  character string to go along with \code{convergence}
  equal to 1 or 10.} \item{terms}{the \code{\link{terms}}
  object used.} \item{data}{the \code{data} argument.}
  \item{xlevels}{(where relevant) a record of the levels of
  the factors used in fitting.}
}
\description{
  Fit a centered autologistic model using maximum
  pseudolikelihood estimation or MCMC for Bayesian
  inference.
}
\details{
  This function fits the centered autologistic model of
  Caragea and Kaiser (2009) using maximum pseudolikelihood
  estimation or MCMC for Bayesian inference.  The joint
  distribution for the centered autologistic model is
  \deqn{\pi(Z\mid\theta)=c(\theta)^{-1}\exp\left(Z^\prime
  X\beta - \eta Z^\prime A\mu + \frac{\eta}{2}Z^\prime
  AZ\right),}{\pi(Z | \theta)=c(\theta)^{-1} exp(Z'X\beta -
  \eta Z'A\mu + 0.5 \eta Z'AZ),} where \eqn{\theta =
  (\beta^\prime, \eta)^\prime}{\theta = (\beta', \eta)'} is
  the parameter vector, \eqn{c(\theta)} is an intractable
  normalizing function, \eqn{Z} is the response vector,
  \eqn{X} is the design matrix, \eqn{\beta} is a
  \eqn{(p-1)}-vector of regression coefficients, \eqn{A} is
  the adjacency matrix for the underlying graph, \eqn{\mu}
  is the vector of independence expectations, and
  \eqn{\eta} is the spatial dependence parameter.  \cr \cr
  Maximum pseudolikelihood estimation sidesteps the
  intractability of \eqn{c(\theta)} by maximizing the
  product of the conditional likelihoods.  Confidence
  intervals are then obtained using a parametric bootstrap.
  The bootstrap datasets are generated by perfect sampling
  (\code{\link{rautologistic}}).  The bootstrap samples can
  be generated in parallel using the \pkg{snow} package.
  \cr \cr Bayesian inference is obtained using the
  auxiliary variable algorithm of Moller et al. (2006).
  The auxiliary variables are generated by perfect
  sampling.  \cr \cr The prior distributions are (1)
  zero-mean normal with independent coordinates for
  \eqn{\beta}, and (2) uniform for \eqn{\eta}.  The
  variance(s) for the normal prior can be supplied by the
  user. The default is a common variance of 1,000,000. The
  uniform prior has support [0, 2] by default, but the
  right endpoint can be supplied (as \code{eta.max}) by the
  user.  \cr \cr The posterior covariance matrix of
  \eqn{\theta} is estimated using samples obtained during a
  training run. The default number of iterations for the
  training run is 100,000, but this can be controlled by
  the user (via argument \code{trainit}). The estimated
  covariance matrix is then used as the proposal variance
  for a Metropolis-Hastings random walk. The proposal
  distribution is normal. The posterior samples obtained
  during the second run are used for inference. The length
  of the run can be controlled by the user via arguments
  \code{minit}, \code{maxit}, and \code{tol}. The first
  determines the minimum number of iterations. If
  \code{minit} has been reached, the sampler will terminate
  when \code{maxit} is reached or all Monte Carlo standard
  errors are smaller than \code{tol}, whichever happens
  first.
}
\references{
  Caragea, P. and Kaiser, M. (2009) Autologistic models
  with interpretable parameters. \emph{Journal of
  Agricultural, Biological, and Environmental Statistics},
  \bold{14}(3), 281--300.

  Hughes, J., Haran, M. and Caragea, P. C. (2011)
  Autologistic models for binary data on a lattice.
  \emph{Environmetrics}, \bold{22}(7), 857--871.

  Moller, J., Pettitt, A., Berthelsen, K., and Reeves, R.
  (2006) An efficient Markov chain Monte Carlo method for
  distributions with intractable normalising constants.
  \emph{Biometrika}, \bold{93}(2), 451--458.
}
\seealso{
  \code{\link{rautologistic}},
  \code{\link{residuals.autologistic}},
  \code{\link{summary.autologistic}},
  \code{\link{vcov.autologistic}}
}


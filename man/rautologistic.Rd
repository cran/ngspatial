\name{rautologistic}
\alias{rautologistic}
\title{Return a perfect sample from a centered autologistic model.}
\usage{
  rautologistic(X, A, theta)
}
\arguments{
  \item{X}{the design matrix.}

  \item{A}{the adjacency matrix for the underlying graph,
  which is assumed to be undirected and free of loops and
  parallel edges.}

  \item{theta}{the vector of parameter values: \eqn{\theta
  = (\beta^\prime, \eta)^\prime}{\theta = (\beta',
  \eta)'}.}
}
\value{
  A vector that is distributed exactly according to the
  centered autologistic model with the given design matrix
  and parameter values.
}
\description{
  Return a perfect sample from a centered autologistic
  model.
}
\details{
  This function implements a perfect sampler for the
  centered autologistic model. The sampler employs coupling
  from the past.
}
\references{
  Moller, J. (1999) Perfect simulation of conditionally
  specified models. \emph{Journal of the Royal Statistical
  Society, Series B, Methodological}, \bold{61}, 251--264.

  Propp, J. G. and Wilson, D. B. (1996) Exact sampling with
  coupled Markov chains and applications to statistical
  mechanics. \emph{Random Structures and Algorithms},
  \bold{9}(1-2), 223--252.
}


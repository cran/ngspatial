\name{rautologistic}
\alias{rautologistic}
\title{Return a perfect sample from a centered autologistic model.}
\usage{
  rautologistic(X, A, theta)
}
\arguments{
  \item{X}{the design matrix.}

  \item{A}{the adjacency matrix for the underlying graph.}

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
\examples{
\dontrun{

# Use the 20 x 20 square lattice as the underlying graph.

n = 20
A = adjacency.matrix(n)

# Assign coordinates to each vertex such that the coordinates are restricted to the unit square
# centered at the origin.

x = rep(0:(n - 1) / (n - 1), times = n) - 0.5
y = rep(0:(n - 1) / (n - 1), each = n) - 0.5
X = cbind(x, y)                                 # Use the vertex locations as spatial covariates.
beta = c(2, 2)                                  # These are the regression coefficients.
mu = exp(X \%*\% beta)
mu = mu / (1 + mu)                              # Compute the independence expectations.
col1 = "white"
col2 = "black"
colfunc = colorRampPalette(c(col1, col2))

# Now produce a level plot of the independence expectations. This shows that the large-scale
# structure corresponds to a probability gradient that increases as one moves from southwest
# to northeast.

dev.new()
levelplot(mu ~ x * y, aspect = "iso", col.regions = colfunc(n^2))

# Simulate a dataset with the above mentioned regression component and eta equal to 0.6. This
# value of eta corresponds to dependence that is moderate in strength.

theta = c(beta, eta = 0.6)
Z = rautologistic(X, A, theta)

# Create a level plot of the simulated data.

dev.new()
levelplot(Z ~ x * y, aspect = "iso", col.regions = c("white", "black"), colorkey = FALSE)
}
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



#' Return a perfect sample from a centered autologistic model.
#'
#' @details This function implements a perfect sampler for the centered autologistic model. The sampler employs coupling from the past. 
#' @param X the design matrix.
#' @param A the adjacency matrix for the underlying graph, which is assumed to be undirected and free of loops and parallel edges.
#' @param theta the vector of parameter values: \eqn{\theta = (\beta^\prime, \eta)^\prime}{\theta = (\beta', \eta)'}.
#' @return A vector that is distributed exactly according to the centered autologistic model with the given design matrix and parameter values.
#' @references
#' Moller, J. (1999) Perfect simulation of conditionally specified models. \emph{Journal of the Royal Statistical Society, Series B, Methodological}, \bold{61}, 251--264.
#' @references
#' Propp, J. G. and Wilson, D. B. (1996) Exact sampling with coupled Markov chains and applications to statistical mechanics. \emph{Random Structures and Algorithms}, \bold{9}(1-2), 223--252.
#' @export

rautologistic = function(X, A, theta)
{
	perfsampler$rautologistic(X, A, theta)
}

autologistic.bmse = function(mat)
{
	if (sum(is.na(mat)) > 0)
	{
		warning("The sample contains NAs.")
		bmse = NA
	}
	else
	{
		mat = as.matrix(mat)
		p = ncol(mat)
		bmse = c()
		for (i in 1:p)
		{
			temp = perfsampler$bmse(mat[, i])
			if (temp == -1) 
				temp = NA
			bmse = c(bmse, temp)
		}
	}
	bmse
}

autologistic.quantile = function(p, Xbeta_i, etasumN_i)
{
    q_i = 1 / (1 + exp(Xbeta_i + etasumN_i))
    if (p > q_i)
        return(1)
    0   
}

autologistic.objective = function(params, X, A, Z)
{
    p = length(params)
    -autologistic.logPL(X, A, Z, c(params[-p], exp(params[p])))
}

autologistic.logPL = function(X, A, Z, theta)
{
    p = length(theta)
    beta = theta[-p]
    eta = theta[p]
    Xbeta = X %*% beta
    mu = exp(Xbeta) / (1 + exp(Xbeta))
    logPL = Xbeta + eta * A %*% (Z - mu)
    logPL = t(Z) %*% logPL - sum(log(1 + exp(logPL)))
    logPL
}

autologistic.fit = function(X, A, Z, optit = 1000)
{
	start = glm(Z ~ X - 1, family = binomial)$coef
    opt = try(optim(c(start, 0), autologistic.objective, X = X, A = A, Z = Z, control = list(maxit = optit)),
              silent = TRUE)
    if (class(opt) == "try-error")
    {
    	coefficients = NULL
    	fitted.values = NULL
    	linear.predictors = NULL
    	residuals = NULL
    	convergence = NULL
    	message = opt[1]
    }
    else
    {
        convergence = opt$convergence
		if (convergence == 1)
			message = "optim iteration limit 'optit' was reached."
		else if (convergence == 10)
			message = "The Nelder-Mead simplex degenerated."
		else
			message = NULL
		p = ncol(X) + 1
        coefficients = opt$par
        coefficients[p] = exp(coefficients[p])
        names(coefficients)[p] = "eta"
        Xbeta = X %*% coefficients[-p]
        mu = exp(Xbeta)
        mu = mu / (1 + mu)
        autocovariate = A %*% (Z - mu)
        linear.predictors = Xbeta + coefficients[p] * autocovariate
        fitted.values = exp(linear.predictors)
        fitted.values = fitted.values / (1 + fitted.values)
        residuals = Z - fitted.values
    }
    object = list(coefficients = coefficients, fitted.values = fitted.values,
                  linear.predictors = linear.predictors, residuals = residuals, 
				  convergence = convergence, message = message, optit = optit)
    class(object) = "autologistic"
    object
}

autologistic.boothelper = function(dummy, X, A, theta, optit)
{
    Z = rautologistic(X, A, theta)
    fit = autologistic.fit(X, A, Z, optit)
    fit
}

autologistic.bootstrap = function(X, A, theta, type, bootit, parallel, nodes, optit)
{
	boot.sample = data.frame(matrix(, bootit, length(theta)))
	if (! parallel)
	{
		for (j in 1:bootit)
		{
			fit = autologistic.boothelper(NULL, X, A, theta, optit)
			temp = rep(NA, length(theta))
			if (is.null(fit$convergence) || fit$convergence != 0)
			    warning(fit$message)
			else
			    temp = fit$coef
			boot.sample[j, ] = temp
	    }
	}
	else
	{
		cl = makeCluster(nodes, type)
		clusterSetupRNG(cl, seed = 1000000 * runif(6))
		clusterEvalQ(cl, library(ngspatial))
		gathered = clusterApplyLB(cl, 1:bootit, autologistic.boothelper, X, A, theta, optit)
		stopCluster(cl)
		for (j in 1:bootit)
		{
			fit = gathered[[j]]
			temp = rep(NA, length(theta))
			if (is.null(fit$convergence) || fit$convergence != 0)
			    warning(fit$message)
			else
			    temp = fit$coef
			boot.sample[j, ] = temp			
		}
   	}
	boot.sample
}

#' Fit a centered autologistic model using maximum pseudolikelihood estimation or MCMC for Bayesian inference.
#'
#' @details This function fits the centered autologistic model of Caragea and Kaiser (2009) using maximum pseudolikelihood estimation or MCMC for Bayesian inference. 
#'			The joint distribution for the centered autologistic model is 
#'			\deqn{\pi(Z\mid\theta)=c(\theta)^{-1}\exp\left(Z^\prime X\beta - \eta Z^\prime A\mu + \frac{\eta}{2}Z^\prime AZ\right),}{\pi(Z | \theta)=c(\theta)^{-1} exp(Z'X\beta - \eta Z'A\mu + 0.5 \eta Z'AZ),}
#'			where \eqn{\theta = (\beta^\prime, \eta)^\prime}{\theta = (\beta', \eta)'} is the parameter vector, \eqn{c(\theta)} is an intractable normalizing function, \eqn{Z} is the response vector, \eqn{X} is the design matrix, 
#'			\eqn{\beta} is a \eqn{(p-1)}-vector of regression coefficients, \eqn{A} is the adjacency matrix for the underlying graph, \eqn{\mu} is the vector of independence expectations, 
#'			and \eqn{\eta} is the spatial dependence parameter. 
#' 			\cr
#'			\cr
#'			Maximum pseudolikelihood estimation sidesteps the intractability of \eqn{c(\theta)} by maximizing the product of the conditional likelihoods.
#'			Confidence intervals are then obtained using a parametric bootstrap. The bootstrap datasets are generated by perfect sampling (\code{\link{rautologistic}}).
#'          The bootstrap samples can be generated in parallel using the \pkg{snow} package.
#' 			\cr
#'			\cr
#' 			Bayesian inference is obtained using the auxiliary variable algorithm of Moller et al. (2006).
#'			The auxiliary variables are generated by perfect sampling.
#'          \cr
#'          \cr
#'          The prior distributions are (1) zero-mean normal with independent coordinates for \eqn{\beta}, and (2) uniform for \eqn{\eta}.
#'          The variance(s) for the normal prior can be supplied by the user. The default is a common variance of 1,000,000. The uniform prior has support [0, 2] by default, but the right endpoint can be supplied (as \code{eta.max}) by the user. 
#'			\cr
#'			\cr
#'			The posterior covariance matrix of \eqn{\theta} is estimated using samples obtained during a training run. The default number of iterations for the training run is 100,000, but this can be controlled by the user (via argument \code{trainit}). The estimated covariance matrix is then used as the proposal variance for a Metropolis-Hastings random walk. The proposal distribution is normal. The posterior samples obtained during the second run are used for inference. The length of the run can be controlled by the user via arguments \code{minit}, \code{maxit}, and \code{tol}. The first determines the minimum number of iterations. If \code{minit} has been reached, the sampler will terminate when \code{maxit} is reached or all Monte Carlo standard errors are smaller than \code{tol}, whichever happens first.
#'
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model to be fitted.
#' @param data an optional data frame, list, or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from \code{environment(formula)}, typically the environment from which \code{autologistic} is called.
#' @param A the adjacency matrix for the underlying graph, which is assumed to be undirected and free of loops and parallel edges.
#' @param method the method to use for inference. \dQuote{\code{pl}} (the default) enables maximum pseudolikelihood estimation, and \dQuote{\code{bayes}} enables Bayesian inference.
#' @param optit the maximum number of iterations to be used by \code{\link{optim}} in obtaining the MPLE estimate of \eqn{\theta}. Defaults to 1,000.
#' @param parallel (PL) a boolean variable indicating whether to use parallel bootstrapping, which requires the \pkg{snow} package. Defaults to \code{TRUE}, in which case the number of nodes must be supplied.
#' @param nodes (PL) the number of nodes to use for parallel bootstrapping.
#' @param type (PL) the type of cluster to use for parallel bootstrapping. The available types are \dQuote{\code{SOCK}}, \dQuote{\code{PVM}}, \dQuote{\code{MPI}}, and \dQuote{\code{NWS}}. The default type is \dQuote{\code{SOCK}}.
#' @param bootit (PL) the desired size of the bootstrap sample. Defaults to 1,000.
#' @param trainit (Bayes) the number of iterations to use for estimating the posterior covariance matrix. Defaults to 100,000. 
#' @param tol (Bayes) a tolerance. If all Monte Carlo standard errors are smaller than \code{tol}, no more samples are drawn from the posterior. Defaults to 0.01.
#' @param minit (Bayes) the minimum sample size. This should be large enough to permit accurate estimation of Monte Carlo standard errors. Defaults to 10,000.
#' @param maxit (Bayes) the maximum sample size. Sampling from the posterior terminates when all Monte Carlo standard errors are smaller than \code{tol} or when \code{maxit} samples have been drawn, whichever happens first. Defaults to 1,000,000.
#' @param sigma (Bayes) a scalar or a \eqn{(p-1)}-vector providing the variance(s) of the spherical normal prior for \eqn{\beta}. Defaults to 1,000,000.
#' @param eta.max (Bayes) the upper limit for \eqn{\eta}. Defaults to 2. The lower limit is 0.
#' @param model a logical value indicating whether the model frame should be included as a component of the returned value.
#' @param x a logical value indicating whether the model matrix used in the fitting process should be returned as a component of the returned value.
#' @param y a logical value indicating whether the response vector used in the fitting process should be returned as a component of the returned value.
#' @return \code{autologistic} returns an object of class \dQuote{\code{autologistic}}, which is a list containing the following components.
#'         \item{coefficients}{the point estimate of \eqn{\theta}.}
#'         \item{fitted.values}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
#'         \item{linear.predictors}{the linear fit on link scale.}
#'         \item{residuals}{the response residuals.}
#'         \item{iter}{the size of the bootstrap/posterior sample.}
#'         \item{sample}{an \code{iter} by \eqn{p} matrix containing the bootstrap/posterior samples.}
#'         \item{mcse}{a \eqn{p}-vector of Monte Carlo standard errors.}
#'         \item{V}{(Bayes) the estimated posterior covariance matrix from the training run.}
#'         \item{accept}{(Bayes) the acceptance rate for the MCMC sampler.}
#'         \item{y}{if requested (the default), the \code{y} vector used.}
#'         \item{X}{if requested, the model matrix.}
#'         \item{model}{if requested (the default), the model frame.}
#'         \item{call}{the matched call.}
#'         \item{formula}{the formula supplied.}
#'         \item{method}{the method used for inference.}
#'         \item{convergence}{an integer code. The code has value 0 if \code{\link{optim}} succeeded in optimizing the pseudolikelihood. Possible error codes are 1 and 10. The former indicates that the iteration limit was reached before optimization completed. The latter indicates that the Nelder-Mead simplex degenerated.}
#'         \item{message}{a character string to go along with \code{convergence} equal to 1 or 10.}
#'         \item{terms}{the \code{\link{terms}} object used.}
#'         \item{data}{the \code{data} argument.}
#'         \item{xlevels}{(where relevant) a record of the levels of the factors used in fitting.}
#' @references
#' Caragea, P. and Kaiser, M. (2009) Autologistic models with interpretable parameters. \emph{Journal of Agricultural, Biological, and Environmental Statistics}, \bold{14}(3), 281--300.
#' @references
#' Hughes, J., Haran, M. and Caragea, P. C. (2011) Autologistic models for binary data on a lattice. \emph{Environmetrics}, \bold{22}(7), 857--871. 
#' @references
#' Moller, J., Pettitt, A., Berthelsen, K., and Reeves, R. (2006) An efficient Markov chain Monte Carlo method for distributions with intractable normalising constants. \emph{Biometrika}, \bold{93}(2), 451--458.
#' @seealso \code{\link{rautologistic}}, \code{\link{residuals.autologistic}}, \code{\link{summary.autologistic}}, \code{\link{vcov.autologistic}}
#' @export

autologistic = function(formula, data, A, method = c("pl", "bayes"), optit = 1000, model = TRUE, x = FALSE, y = FALSE,
						type = c("SOCK", "PVM", "MPI", "NWS"), bootit = 1000, parallel = TRUE, nodes,
						trainit = 100000, tol = 0.01, minit = 10000, maxit = 1000000, sigma = 1000000, eta.max = 2)
{
	cl = match.call()	
	if (missing(formula))
		stop("You must supply a formula.")
	if (missing(A) || ! isSymmetric(A) || ! (A == 0 || A == 1) || ! is.matrix(A))
		stop("You must supply a symmetric binary adjacency matrix.")
    diag(A) = rep(0, nrow(A))
	if (! is.numeric(optit) || ! is.wholenumber(optit) || optit < 0)
		stop("'optit' must be a positive whole number.")
	default.args.pl = list(bootit = bootit, type = match.arg(type), parallel = parallel, nodes = NULL)
	default.args.bayes = list(trainit = trainit, tol = tol, minit = minit, maxit = maxit, sigma = sigma, eta.max = eta.max)
	mf = match.call(expand.dots = FALSE)
	m = match(c("formula", "data"), names(mf), 0)
	mf = mf[c(1, m)]
	mf[[1]] = as.name("model.frame")
	mf = eval(mf, parent.frame())
	mt = attr(mf, "terms")
	Z = model.response(mf, "numeric")
	X = model.matrix(mt, mf)
	if (! is.vector(Z) || ! (Z == 0 || Z == 1))
		stop("The response should be a binary vector.")
	if (sum(c(nrow(X), nrow(A),ncol(A)) != length(Z)) > 0)
		stop("The supplied response vector/design matrix/adjacency matrix are not conformable.")
	cl$method = match.arg(method)
	cl$optit = optit
	fit = autologistic.fit(X, A, Z, optit)
	if (fit$convergence != 0)
		stop(fit$message, " Consider increasing the value of 'optit'.")
	if (cl$method == "pl")
	{
		cl = ngspatial.silent(cl, names(default.args.bayes))
		if (! is.numeric(bootit) || ! is.wholenumber(bootit) || bootit < 1)
			stop("'bootit' must be a positive whole number.")
		if (is.null(parallel))
			cat("Parallel bootstrapping is enabled by default.\n")
		if (parallel)
		{
			if (missing(nodes))
				stop("You must supply the number of nodes to use for the parallel bootstrap.")
			if (is.null(nodes) || ! is.numeric(nodes) || ! is.wholenumber(nodes) || nodes < 1)
				stop("'nodes' must be a positive whole number.")
			if (! require(snow))
				stop("Parallel bootstrapping requires the snow package.")
		}
		else
		{
			if (! is.null(cl$nodes) && cl$nodes != 1)
				cat("Parallel bootstrapping is disabled. Argument 'nodes' will be ignored.\n")
			nodes = 1
		}
		default.args.pl$nodes = nodes
		cl = ngspatial.existOrFill(cl, names(default.args.pl), default.args.pl)
		if (cl$nodes == 1)
			cl$parallel = FALSE
		fit$sample = autologistic.bootstrap(X, A, fit$coef, cl$type, cl$bootit, cl$parallel, cl$nodes, cl$optit)
		fit$mcse = autologistic.bmse(fit$sample)
		fit$iter = cl$bootit
		if (! cl$parallel)
			cl$nodes = NULL
	}
	else # cl$method == "bayes"
	{
		cl = ngspatial.silent(cl, names(default.args.pl))
		if (! is.numeric(trainit) || ! is.wholenumber(trainit) || trainit < 0)
			stop("'trainit' must be a positive whole number.")
		if (trainit < 10000)
			 warning("Consider increasing the value of 'trainit'.")
		if (! is.numeric(tol) || tol <= 0 || tol >= 1)
			stop("'tol' must be between 0 and 1.")
		if (! is.numeric(minit) || ! is.wholenumber(minit) || minit < 1)
			stop("'minit' must be a positive whole number.")
		if (minit < 10000)
		    warning("Consider increasing the value of 'minit'.")
		if (! is.numeric(maxit) || ! is.wholenumber(maxit) || maxit < minit)
			stop("'maxit' must be a positive whole number greater than 'minit'.")
		if (! is.numeric(sigma) || sum(sigma < 0) > 0)
			stop("'sigma' must be a positive scalar or a vector with only positive elements.")
		if (! is.numeric(eta.max) || eta.max <= 0)
			stop("'eta.max' must be a positive number.")
	    if (eta.max < 1)
	        warning("Consider increasing the value of 'eta.max'.")
		cl = ngspatial.existOrFill(cl, names(default.args.bayes), default.args.bayes)
		p = length(fit$coef) - 1
		if (length(cl$sigma) == 1)
			cl$sigma = rep(cl$sigma, p)
		else if (length(cl$sigma) != p)
			stop("'sigma' must be a scalar or a vector of length ", p)
		cat("\nWarning: MCMC may be time consuming.\n\n")
		flush.console()
		Sys.sleep(0.5)
		temp = Moller.run(X, A, Z, fit$coef, cl$trainit, cl$tol, cl$minit, cl$maxit, cl$sigma, cl$eta.max)
		fit$coefficients = NULL
		fit = c(fit, temp)
		fit$accept = sum(diff(fit$sample[, 1]) != 0) / fit$iter
		class(fit) = c("autologistic")
	}
	fit$xlevels = .getXlevels(mt, mf)
	fit$call = cl
	fit$terms = mt
	fit$method = cl$method
	if (model)
		fit$model = mf
	if (x)
		fit$x = X
	if (y)
		fit$y = Z
	fit
}

#' Extract model residuals.
#'
#' @param object an object of class \code{autologistic}, typically the result of a call to \code{\link{autologistic}}.
#' @param type the type of residuals that should be returned. The alternatives are \dQuote{\code{deviance}} (default), \dQuote{\code{pearson}}, and \dQuote{\code{response}}.
#' @param \dots additional arguments.
#' @return A vector of residuals.
#' @seealso \code{\link{autologistic}}, \code{\link{residuals.glm}}
#' @method residuals autologistic
#' @export

residuals.autologistic = function(object, type = c("deviance", "pearson","response"), ...)
{
	type = match.arg(type)
    if (type == "response")
        return(object$residuals)
    else if (type == "deviance")
    {
    	if (is.null(object$y))
    	    y = object$residuals + object$fitted.values
    	phat = object$fitted.values
    	d = numeric(length(y))
    	zero = which(y == 0)
    	d[zero] = -2 * log(1 - phat[zero])
    	one = which(y == 1)
    	d[one] = -2 * log(phat[one])
        return(sqrt(d) * sign(object$residuals))
    }
    else # type == "pearson"
    {
    	phat = object$fitted.values
    	se = sqrt(phat * (1 - phat))
    	return(object$residuals / se)
    }
}

#' Return the estimated covariance matrix for an \code{autologistic} model object.
#'
#' @param object a fitted \code{autologistic} model object.
#' @param \dots additional arguments.
#' @return An estimate of the covariance matrix of the parameters (in a Bayesian setting), or an estimate of the covariance matrix of the maximum pseudolikelihood estimator of the parameters. The latter requires a bootstrap sample. If there is no such sample, \code{NULL} is returned.
#' @method vcov autologistic
#' @export

vcov.autologistic = function(object, ...)
{
	if (is.null(object$sample))
	    V = NULL
	else if (sum(is.na(object$sample)) > 0)
		stop("The sample contains NAs.")
	else
	{
		V = cov(object$sample)
		rownames(V) = colnames(V) = names(object$coef)
	}
	V
}

#' Print a summary of a centered autologistic model fit.
#'
#' @details This function displays (1) the call to \code{\link{autologistic}}, (2) a table of estimates, and (3) the size of the bootstrap/posterior samples. Each row of the table of estimates shows the estimated regression coefficient, the \eqn{(\alpha/2)100\%}{(\alpha/2)100\%} and \eqn{(1-\alpha/2)100\%}{(1-\alpha/2)100\%} bootstrap quantiles for the MPLE coefficient or \eqn{(1-\alpha)100\%}{(1-\alpha)100\%} HPD interval for the Bayesian coefficient, and the Monte Carlo standard error.
#' @param object an object of class \code{autologistic}, typically the result of a call to \code{\link{autologistic}}.
#' @param alpha the significance level for the quantile/HPD intervals. The default is 0.05.
#' @param digits the number of significant digits to display. The default is 4.
#' @param \dots additional arguments.
#' @seealso \code{\link{autologistic}}
#' @method summary autologistic
#' @export

summary.autologistic = function(object, alpha = 0.05, digits = 4, ...)
{
	cat("\nCall:\n")
	print(object$call)
	p = length(object$coef)
    if (is.null(object$sample) || (sum(is.na(object$sample)) > 0))
	{
		warning("The sample is NULL or contains NAs.")
    	coef.table = cbind(object$coef, NA, NA, NA)
	}
	else
	{
		ci = matrix(, p, 2)
		if (object$method == "pl")
		{
			for (j in 1:p)
				ci[j, ] = quantile(object$sample[, j], c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
		}
		else if(object$method == "bayes")
    	{
			for (j in 1:p)
				ci[j, ] = hpd(object$sample[, j], alpha)
    	}
		coef.table = cbind(object$coef, ci, object$mcse)
	}
    colnames(coef.table) = c("Estimate", "Lower", "Upper", "MCSE")
	rownames(coef.table) = names(object$coef)
	cat("\nCoefficients:\n")
	print(signif(coef.table, digits))
	cat("\nNumber of iterations:", object$iter, "\n\n")
}





rand.walk = function(X, A, Z, theta, tol, minit, maxit, sigma, eta.max, V)
{	
	p = length(theta) - 1
	if (is.null(sigma))
		sigma = rep(1000000, p)
	else if (length(sigma) != p)
		stop("'sigma' must be a vector of length ", p)
	eta.range = c(0, eta.max)
	moller$randWalk(X, A, Z, theta, tol, minit, maxit, sigma, eta.range, V)
}

rand.walk.train = function(X, A, Z, theta, trainit, sigma, eta.max, V)
{	
	p = length(theta) - 1
	if (is.null(sigma))
		sigma = rep(1000000, p)
	else if (length(sigma) != p)
		stop("'sigma' must be a vector of length ", p)
	eta.range = c(0, eta.max)
	moller$randWalkTrain(X, A, Z, theta, trainit, sigma, eta.range, V)
}

Moller.train = function(X, A, Z, theta, trainit, sigma, eta.max)
{
	p = length(theta)
	V.init = diag(0.01, p)
	V = matrix(0, p, p)
    if (trainit >= 1)
	{
		estimates = rand.walk.train(X, A, Z, theta, trainit, sigma, eta.max, V.init)
		V = cov(estimates)
    }
	if (is.zero(V) || is.na(V))
		V = V.init
	V
}

Moller.run = function(X, A, Z, theta, trainit, tol, minit, maxit, sigma, eta.max)
{
    p = length(theta)
	mcse = numeric(p)
	coefficients = numeric(p)
	result = list()
	V = Moller.train(X, A, Z, theta, trainit, sigma, eta.max)
	estimates = data.frame(rand.walk(X, A, Z, theta, tol, minit, maxit, sigma, eta.max, V))
    for (i in 1:p)
    {
    	temp = bm(estimates[, i])
    	coefficients[i] = temp$est
    	mcse[i] = temp$se
    }
	iter = nrow(estimates)
	list(sample = estimates, coefficients = coefficients, mcse = mcse, iter = iter, V = V)
}


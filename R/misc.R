
#' Return an adjacency matrix for a square lattice.
#'
#' @details This function builds the adjacency matrix for the \code{m} by \code{n} square lattice.
#' @param m the number of rows in the lattice.
#' @param n the number of columns in the lattice. Defaults to \code{NULL}. If missing, the lattice is assumed to be \code{m} by \code{m}. 
#' @return A matrix \eqn{A} of 0s and 1s, where \eqn{A_{ij}} is equal to 1 iff vertices \eqn{i} and \eqn{j} are adjacent.
#' @export

adjacency.matrix = function(m, n = NULL)
{
    if (missing(n))
    {
        A = matrix(0, m^2, m^2)
        for (i in 1:m^2)
        {
            up = i - m
            down = i + m
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= m^2)
                A[i, down] = 1
            if (left %% m != 0)
                A[i, left] = 1
            if (i %% m != 0)
                A[i, right] = 1
        }
    }
    else
    {
        A = matrix(0, m * n, m * n)
        for (i in 1:(m * n))
        {
            up = i - n
            down = i + n
            left = i - 1
            right = i + 1
            if (up > 0)
                A[i, up] = 1
            if (down <= (m * n))
                A[i, down] = 1
            if (left %% n != 0)
                A[i, left] = 1
            if (i %% n != 0)
                A[i, right] = 1
        }
    }
    A
}

hpd = function(x, alpha = 0.05)
{
    n = length(x)
    m = round(n * alpha)
    x = sort(x)
    y = x[(n - m + 1):n] - x[1:m]
    z = min(y)
    k = which(y == z)[1]
    c(x[k], x[n - m + k])
}

is.wholenumber = function(x, tol = .Machine$double.eps^0.5)
{
    abs(x - round(x)) < tol
}

is.zero = function(x, tol = .Machine$double.eps^0.5)
{
    abs(x) < tol
}

matmult.par = function(cl, A, B)
{
    if (ncol(A) != nrow(B))
        stop("The supplied matrices are not conformable.")
    ind = parallel::splitIndices(nrow(A), length(cl))
	Alist = lapply(ind, function(ii) A[ii, , drop = FALSE])
	ans = parallel::clusterApply(cl, Alist, get("%*%"), B)
	do.call(rbind, ans)
}

#' Family function for negative binomial GLMs.
#'
#' @description Provides the information required to apply a sparse SGLMM to conditionally negative binomial outcomes.
#'
#' @usage negbinomial(theta = stop("'theta' must be specified."), link = "log")
#'
#' @param theta the dispersion parameter (must be positive).
#' @param link the link function, as a character string, name, or one-element character vector, specifying one of \code{log}, \code{sqrt}, or \code{identity}, or an object of class \dQuote{\code{\link[=family]{link-glm}}}
#'
#' @return An object of class \dQuote{\code{family}}, a list of functions and expressions needed to fit a negative binomial GLM.
#'
#' @export

negbinomial = function(theta = stop("'theta' must be specified."), link = "log") 
{
    if (! is.vector(theta, mode = "numeric") || length(theta) > 1)
        stop("You must supply a scalar value for 'theta'.")
    if (theta <= 0)
        stop("'theta' must be positive.")
    linktemp = substitute(link)
    if (! is.character(linktemp)) 
        linktemp = deparse(linktemp)
    if (linktemp %in% c("log", "identity", "sqrt")) 
        stats = make.link(linktemp)
    else if (is.character(link))
    {
        stats = make.link(link)
        linktemp = link
    }
    else
    {
        if (inherits(link, "link-glm"))
        {
            stats = link
            if (! is.null(stats$name)) 
                linktemp = stats$name
        }
        else stop(gettextf("\"%s\" link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"", linktemp))
    }
    variance = function(mu) mu + mu^2 / theta
    validmu = function(mu) all(mu > 0)
    dev.resids = function(y, mu, wt) 2 * wt * (y * log(pmax(1, y) / mu) - (y + theta) * log((y + theta) / (mu + theta)))
    aic = function(y, n, mu, wt, dev)
    {
        term = (y + theta) * log(mu + theta) - y * log(mu) + lgamma(y + 1) - theta * log(theta) + lgamma(theta) - lgamma(theta + y)
        2 * sum(term * wt)
    }
    initialize = expression({
        if (any(y < 0))
            stop("negative values not allowed for the negative binomial family")
        n = rep(1, nobs)
        mustart = y + (y == 0) / 6
    })
    famname = "negbinomial"
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun, 
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize, theta = theta,
        validmu = validmu, valideta = stats$valideta), 
        class = "family")
}


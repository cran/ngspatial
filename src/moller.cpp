#include "moller.h"
#include "perfsampler.h"
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace std;
using namespace arma;
 
mat randWalk(const mat& X, const mat& A, const colvec& Z, 
             const colvec& theta, double tol, int minit, int maxit,
             const colvec& sigma, const colvec& etaRange,
             const mat& V, bool verbose)
{
    double logAccept;
    int iterations = 0;
    int p = theta.size() - 1;
    int n = X.n_rows;
    colvec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, V);
    mat R = eigvec * diagmat(sqrt(eigval));
    colvec Y = rautologistic(X, A, theta);
    colvec Ynew(n);
    colvec thetaold = theta;
    colvec thetanew = theta;
    colvec normvec(p + 1);
    mat estimates(maxit + 1, p + 1);
    estimates.row(0) = theta.t();
    colvec beta = theta.subvec(0, p - 1);
    double eta = theta[p];
    colvec Xbeta = X * beta;
    colvec mu = exp(Xbeta);
    mu = mu / (1 + mu);
    double Zsum = as_scalar(Z.t() * A * Z);
    int fiveper = 0.05 * maxit;
	RNGScope scope;
	Function rnorm("rnorm"),
	         runif("runif");
    do
    {
        do
        {
            normvec = as<colvec>(rnorm(p + 1));
            thetanew = thetaold + R * normvec;
        }
        while (thetanew[p] < etaRange[0] || thetanew[p] > etaRange[1]);
        Ynew = rautologistic(X, A, thetanew);
        colvec betaold = thetaold.subvec(0, p - 1);
        colvec betanew = thetanew.subvec(0, p - 1);
        double etaold = thetaold[p];
        double etanew = thetanew[p];
        colvec Xbetaold = X * betaold;
        colvec Xbetanew = X * betanew;
        colvec muold = exp(Xbetaold);
        muold = muold / (1 + muold);
        colvec munew = exp(Xbetanew );
        munew = munew / (1 + munew);
        double Ynewsum = as_scalar(Ynew.t() * A * Ynew);
        double Ysum = as_scalar(Y.t() * A * Y);
        logAccept = as_scalar(0.5 * (eta * (Ynewsum - Ysum) + etanew * (Zsum - Ynewsum) + etaold * (Ysum - Zsum))
                    + trans(Ynew - Y) * Xbeta + trans(Z - Ynew) * Xbetanew + trans(Y - Z) * Xbetaold
                    + etaold * trans(Z - Y) * A * muold + etanew * trans(Ynew - Z) * A * munew + eta * trans(Y - Ynew) * A * mu)
                    + accu((square(betaold) - square(betanew)) / (2 * sigma));
        if (as_scalar(log(as<double>(runif(1)))) < logAccept)
        {
            estimates.row(iterations + 1) = thetanew.t();
            thetaold = thetanew;
            Y = Ynew;
        }
        else
            estimates.row(iterations + 1) = thetaold.t();
        iterations++;
        if (verbose && iterations % fiveper == 0)
        {
            Rcout << "Progress => " << double(iterations) / double(maxit) * 100 << "%" << std::endl;
            R_FlushConsole();
        }
    }
    while ((! isLessTol(estimates.rows(1, iterations), tol) || iterations < minit) && iterations < maxit); 
    return estimates.rows(1, iterations);
}

mat randWalkTrain(const mat& X, const mat& A, const colvec& Z,  
                  const colvec& theta, int trainit, 
                  const colvec& sigma, const colvec& etaRange,
                  const mat& V)
{
    double logAccept;
    int p = theta.size() - 1;
    int n = X.n_rows;
    colvec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, V);
    mat R = eigvec * diagmat(sqrt(eigval));
    colvec Y = rautologistic(X, A, theta);
    colvec Ynew(n);
    colvec thetaold = theta;
    colvec thetanew = theta;
    colvec normvec(p + 1);
    mat estimates(trainit + 1, p + 1);
    estimates.row(0) = theta.t();
    colvec beta = theta.subvec(0, p - 1);
    double eta = theta[p];
    colvec Xbeta = X * beta;
    colvec mu = exp(Xbeta);
    mu = mu / (1 + mu);
    double Zsum = as_scalar(Z.t() * A * Z);
	RNGScope scope;
	Function rnorm("rnorm"),
	         runif("runif");
    for (int i = 1; i < trainit + 1; i++)
    {
        do
        {
            normvec = as<colvec>(rnorm(p + 1));
            thetanew = thetaold + R * normvec;
        }
        while(thetanew[p] < etaRange[0] || thetanew[p] > etaRange[1]);
        Ynew = rautologistic(X, A, thetanew);
        colvec betaold = thetaold.subvec(0, p - 1);
        colvec betanew = thetanew.subvec(0, p - 1);
        double etaold = thetaold[p];
        double etanew = thetanew[p];
        colvec Xbetaold = X * betaold;
        colvec Xbetanew = X * betanew;
        colvec muold = exp(Xbetaold);
        muold = muold / (1 + muold);
        colvec munew = exp(Xbetanew );
        munew = munew / (1 + munew);
        double Ynewsum = as_scalar(Ynew.t() * A * Ynew);
        double Ysum = as_scalar(Y.t() * A * Y);
        logAccept = as_scalar(0.5 * (eta * (Ynewsum - Ysum) + etanew * (Zsum - Ynewsum) + etaold * (Ysum - Zsum))
                            + trans(Ynew - Y) * Xbeta + trans(Z - Ynew) * Xbetanew + trans(Y - Z) * Xbetaold
                            + etaold * trans(Z - Y) * A * muold + etanew * trans(Ynew - Z) * A * munew + eta * trans(Y - Ynew) * A * mu)
                            + accu((square(betaold) - square(betanew)) / (2 * sigma));
        if (as_scalar(log(as<double>(runif(1)))) < logAccept)
        {
            estimates.row(i) = thetanew.t();
            thetaold = thetanew;
            Y = Ynew;
        }
        else
            estimates.row(i) = thetaold.t();
    }
    return estimates.rows(1, trainit);
} 

RCPP_MODULE(moller)
{
    Rcpp::function("randWalk", &randWalk);
    Rcpp::function("randWalkTrain", &randWalkTrain);
}

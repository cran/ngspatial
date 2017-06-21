#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
using namespace std;
using namespace arma;

// [[Rcpp::export]]
double bmse(const vec& vals)
{   
    int N = vals.size();
    double result;
    if (N < 10)
        result = -1;
    else
    {
        int b = floor(sqrt(N));
        int a = floor(N / b);
        colvec Ys(a);
        colvec temp(b);
        for (int i = 0; i < a; i++)
        {
            temp = vals.subvec(i * b, (i + 1) * b - 1);
            Ys[i] = mean(temp);
        }
        double muhat = mean(Ys);
        double sigmahatsq = b * as_scalar(sum(pow((Ys - muhat), 2))) / (a - 1);
        result = sqrt(sigmahatsq / N);
    }
    return result;
}

bool isLessTol(const mat& vals, double tol)
{
    int p = vals.n_cols;
    bool lessTol = true;
    for (int i = 0; i < p; i++)
    {   
        if (bmse(vals.col(i)) > tol || bmse(vals.col(i)) < 0)
        {
            lessTol = false;
            break;
        }
    }
    return lessTol;
}

// [[Rcpp::export]]
vec rautologistic_(const mat& X, const mat& A, const vec& theta)
{
    colvec beta = theta.subvec(0, theta.size() - 2);
    double eta = theta[theta.size() - 1];
    int n = X.n_rows;
    colvec Xbeta = X * beta;
    colvec mu = exp(Xbeta);
    mu = mu / (1 + mu);
    colvec summu = A * mu;
    colvec L = zeros<colvec>(n);
    colvec U = ones<colvec>(n);
    int T = 2;
    mat R = zeros<mat>(T, n);
    int t = 0;
    bool restart = false;
    double sumL_i,
           sumU_i,
           q;
    RNGScope scope;		   
	Function runif("runif");
    while (true)
    {
        if (t == T && as_scalar(sum(U - L)) == 0)
            return L;
        t++;
        if (t > T)
        {
            L.zeros();
            U.ones();
            R.insert_rows(0, T);
            T *= 2;
            t = 1;
            restart = true;
        }
        if (! restart || (restart && t <= T / 2))
            R.row(t - 1) = as<rowvec>(runif(n)); //randu<rowvec>(n);
        for (int i = 0; i < n; i++)
        {
            sumL_i = as_scalar(A.row(i) * L);
            q = 1 / (1 + exp(Xbeta[i] + eta * (sumL_i - summu[i])));
            if (R(t - 1, i) > q)
                L[i] = 1;
            else
                L[i] = 0;
            sumU_i = as_scalar(A.row(i) * U);
            q = 1 / (1 + exp(Xbeta[i] + eta * (sumU_i - summu[i])));
            if (R(t - 1, i) > q)
                U[i] = 1;
            else
                U[i] = 0;
        }
    }
    return L;
}

// [[Rcpp::export]]
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
    colvec Y = rautologistic_(X, A, theta);
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
    int onepct = 0.01 * maxit;
	int pct = 1;
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
        Ynew = rautologistic_(X, A, thetanew);
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
        if (verbose && iterations % onepct == 0)
        {
            Rcout << "\r|";
			for (int i = 1; i <= pct; i++)
				Rcout << "+";
			for (int i = 1; i < (100 - pct); i++)
				Rcout << " ";
			Rcout << "| " << pct << "%";
			pct++;
            R_FlushConsole();
        }
    }
    while ((! isLessTol(estimates.rows(1, iterations), tol) || iterations < minit) && iterations < maxit);
	if (verbose)
		Rcout << std::endl;
    return estimates.rows(1, iterations);
}

// [[Rcpp::export]]
mat randWalkTrain(const mat& X, const mat& A, const colvec& Z,  
                  const colvec& theta, int trainit, 
                  const colvec& sigma, const colvec& etaRange,
                  const mat& V, bool verbose)
{
    double logAccept;
    int p = theta.size() - 1;
    int n = X.n_rows;
    colvec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, V);
    mat R = eigvec * diagmat(sqrt(eigval));
    colvec Y = rautologistic_(X, A, theta);
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
    int onepct = 0.01 * trainit;
	int pct = 1;
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
        Ynew = rautologistic_(X, A, thetanew);
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
        if (verbose && i % onepct == 0)
        {
            Rcout << "\rTraining progress: |";
			for (int j = 1; j <= pct; j++)
				Rcout << "+";
			for (int j = 1; j < (100 - pct); j++)
				Rcout << " ";
			Rcout << "| " << pct << "%";
			pct++;
            R_FlushConsole();
        }
    }
	if (verbose)
		Rcout << std::endl << std::endl;
    return estimates.rows(1, trainit);
} 

RCPP_MODULE(ngspatialmod)
{
    Rcpp::function("rautologistic_", &rautologistic_);
    Rcpp::function("bmse", &bmse);
    Rcpp::function("randWalk", &randWalk);
    Rcpp::function("randWalkTrain", &randWalkTrain);
}


#include "perfsampler.h"
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;
using namespace RcppArmadillo;
using namespace std;
using namespace arma;
    
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

vec rautologistic(const mat& X, const mat& A, const vec& theta)
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

RCPP_MODULE(perfsampler)
{
    Rcpp::function("rautologistic", &rautologistic);
    Rcpp::function("bmse", &bmse);
}


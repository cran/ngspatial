#ifndef _ngspatial_MOLLER_H
#define _ngspatial_MOLLER_H

#include <RcppArmadillo.h>

arma::mat randWalk(const arma::mat& X, const arma::mat& A, const arma::colvec& Z,
                   const arma::colvec& theta, double tol, int minit, int maxit, 
                   const arma::colvec& sigma, const arma::colvec& etaRange, const arma::mat& V,
                   bool verbose);
arma::mat randWalkTrain(const arma::mat& X, const arma::mat& A, const arma::colvec& Z,
                        const arma::colvec& theta, int trainit, const arma::colvec& sigma,
                        const arma::colvec& etaRange, const arma::mat& V);
                    
#endif

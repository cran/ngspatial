#ifndef _ngspatial_PERFSAMPLER_H
#define _ngspatial_PERFSAMPLER_H

#include <RcppArmadillo.h>

arma::vec rautologistic(const arma::mat& X, const arma::mat& A, const arma::vec& theta);
double bmse(const arma::vec& vals);
bool isLessTol(const arma::mat& vals, double tol);

#endif

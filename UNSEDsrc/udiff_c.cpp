#include <iostream>
#include <stdio.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


arma::mat fast_dummy(const arma::vec& x) {
  arma::vec vec_u = arma::unique(x);
  vec_u = vec_u(arma::span(1,vec_u.n_elem-1)); // Remove the baseline
  arma::mat mats = arma::zeros(x.n_elem, (vec_u.n_elem));
  for (arma::uword i=0; i<vec_u.n_elem; ++i) {
    mats.col(i) = x;
    double idx = vec_u.at(i);
    mats.col(i).for_each([idx](arma::vec::elem_type& val){
      (val == idx) ? (val = 1) : (val = 0);
    });
  }
  return mats;
}


arma::mat udiff_c(const arma::vec& theta, const arma::mat& Y,
                 const arma::mat& X, const bool& lg = false)
  // The group variable is the last variable of X.
{
  if (Y.n_cols != 1 ) {
    throw std::invalid_argument("Y should only have 1 column!");
  }
  if (X.n_cols != 2 ) {
    throw std::invalid_argument("X should have 2 columns!");
  }
  arma::vec uY = arma::unique(Y);
  arma::mat X_X = fast_dummy(X.col(0));
  arma::mat X_Z = fast_dummy(X.col(1));
  arma::mat Y_Y = fast_dummy(Y.col(0));
  arma::mat W = arma::ones(X_Z.n_rows, (X_Z.n_cols+1));
  W(arma::span::all,arma::span(1,W.n_cols-1)) = X_Z;
  arma::uword colX = X_X.n_cols;
  arma::uword colZ = X_Z.n_cols;
  arma::uword colY = Y_Y.n_cols;
  if (theta.n_elem != (colY*(colZ+1)+colY*colX+colZ)) {
    throw std::invalid_argument("Wrong size of theta!");
  }
  arma::uword end = colY*(colZ+1)-1;
  arma::uword end2 = colY*colX;
  arma::uword end3 = colZ;
  arma::mat theta_y = arma::mat(theta(arma::span(0, end)));
  theta_y.reshape(colZ+1, colY);
  arma::mat psi_y = arma::mat(theta(arma::span(end+1, end+end2)));
  psi_y.reshape(colX, colY);
  arma::mat phi = arma::mat(theta(arma::span(end+end2+1,end+end2+end3)));
  phi.reshape(colZ,1);
  arma::mat expZ = arma::exp(X_Z * phi);
  arma::mat phiX = X_X  * psi_y;
  arma::mat part2 = phiX.each_col() % expZ;
  arma::mat sum_part = W * theta_y + part2;
  arma::mat first = arma::sum(Y_Y % sum_part, 1);
  arma::mat second = log(1+arma::sum(arma::exp(sum_part), 1));
  arma::mat l = first - second;
  return l;
}

#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
//define the main function, MovementRcpp
arma::mat Fast_FOI_function(arma::uvec id,
                            const arma::mat& centroids,
                            const int cells,
                            const double F2_int,
                            const double F2_B,
                            const double F2i_int,
                            const double F2i_B,
                            const arma::mat& Imat,
                            const arma::mat& Cmat,
                            const double B1,
                            const double B2
) {


const int Ninf = id.n_elem;

//get col centroids
arma::mat X_centroids = centroids.col(0);
arma::mat Y_centroids = centroids.col(1);

//get submatrix with locs of infected individuals
arma::mat Xinf = X_centroids.rows(id);
arma::mat Yinf = Y_centroids.rows(id);

//format matrices for matrix multiplication
arma::mat X1 = repmat(Xinf, 1, cells);
arma::mat Y1 = repmat(Yinf, 1, cells);

//format matrices for matrix multiplication
arma::mat X2 = repmat(trans(X_centroids), Ninf, 1);
arma::mat Y2 = repmat(trans(Y_centroids), Ninf, 1);

//piecemeal pythagorean distance formula
//get diffs in x/y locations
arma::mat diffX = X2-X1;
arma::mat diffY = Y2-Y1;

//square diffs
arma::mat diffX_2 = diffX % diffX;
arma::mat diffY_2 = diffY % diffY;

//sqrt added squared diffs
arma::mat dist = sqrt(diffX_2+diffY_2);

//use distances in distanct/contact formula to get prob of contact
//infected live individuals
arma::mat plogit = F2_int + (F2_B*dist);
arma::mat prob = exp(plogit)/(1+exp(plogit));
//infected carcasses
arma::mat pilogit = F2i_int + (F2i_B*dist);
arma::mat probi = exp(pilogit)/(1+exp(pilogit));

//find where dist=0, set those indices in prob/probi to zero
//zero distance contact (same cell) is dealt with outside of Rcpp with simpler probability
arma::uvec dist_zero = find(dist == 0);
prob(dist_zero) = dist(dist_zero);
probi(dist_zero) = dist(dist_zero);

//Get scaled probabilities, combined live infected contact with carcass infected contact
//B1 and B2 are scaling parameters, determined via sensitivity analyses
arma::mat B = B1*(Imat%prob)+B2*(Cmat%probi);

return(B);

}
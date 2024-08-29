#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
//define the main function, MovementRcpp
arma::mat Fast_FOI_function(arma::uvec id,
                            const arma::mat& X_centroids,
                            const arma::mat& Y_centroids,
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

//convert I into transformed repmat version
//need I transformed to col, then repeated along rows by num infected
arma::mat I_cells = repmat(Imat, 1, cells);
arma::mat C_cells = repmat(Cmat, 1, cells);

//get submatrix with locs of infected individuals
arma::mat Xinf = X_centroids.rows(id);
arma::mat Yinf = Y_centroids.rows(id);

//format matrices for matrix multiplication
arma::mat X1 = repmat(Xinf, 1, cells);
arma::mat Y1 = repmat(Yinf, 1, cells);

//format matrices for matrix multiplication
arma::mat X2 = repmat(trans(X_centroids), Ninf, 1);
arma::mat Y2 = repmat(trans(Y_centroids), Ninf, 1);

//sqrt added squared diffs
arma::mat dist = sqrt(((X2-X1) % (X2-X1))+((Y2-Y1) % (Y2-Y1)));

//initiate empty matrix for dist
//arma::mat dist_empty zeros(dist.n_rows,1);

//use distances in distanct/contact formula to get prob of contact
//infected live individuals
//arma::mat plogit = F2_int + (F2_B*dist);

arma::mat prob = zeros(dist.n_rows,dist.n_cols);
arma::mat probi = zeros(dist.n_rows,dist.n_cols);
arma::mat B = zeros(dist.n_rows,dist.n_cols);

for(std::size_t c = 0; c < dist.n_cols; c++){
for(std::size_t r = 0; r < dist.n_rows; r++){

//vectorized versions:
//infected live individuals
//arma::mat prob = exp(F2_int + (F2_B*dist))/(1+exp(F2_int + (F2_B*dist)));
//infected carcasses
//arma::mat probi = exp(F2i_int + (F2i_B*dist))/(1+exp(F2i_int + (F2i_B*dist)));

//loop versions:
if(dist(r,c)<5 && dist(r,c)!=0){
prob(r,c) = exp(F2_int + (F2_B*dist(r,c)))/(1+exp(F2_int + (F2_B*dist(r,c))));
probi(r,c) = exp(F2_int + (F2i_B*dist(r,c)))/(1+exp(F2i_int + (F2i_B*dist(r,c))));
B(r,c) = B1*(I_cells(r,c)*prob(r,c))+B2*(I_cells(r,c)*probi(r,c));

}

}
}

//find where dist=0, set those indices in prob/probi to zero
//zero distance contact (same cell) is dealt with outside of Rcpp with simpler probability
//arma::uvec dist_zero = find(dist == 0);
//prob(dist_zero) = dist(dist_zero);
//probi(dist_zero) = dist(dist_zero);

//Get scaled probabilities, combined live infected contact with carcass infected contact
//B1 and B2 are scaling parameters, determined via sensitivity analyses

//arma::mat B = B1*(I_cells%prob)+B2*(I_cells%probi);

return(B);

}
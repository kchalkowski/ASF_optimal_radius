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

//initiate empty matrices for looping
//arma::mat prob = zeros(dist.n_rows,dist.n_cols);
//arma::mat probi = zeros(dist.n_rows,dist.n_cols);
arma::mat B = zeros(dist.n_rows,dist.n_cols);

//cols in the outer loop is slightly faster
//Note: below loop is structured this way to make use of loop unrolling optimization
for(std::size_t c = 0; c < dist.n_cols-2; c+=2){
for(std::size_t r = 0; r < dist.n_rows; r++){

//set limits on which probabilities need to be calculated
//far end of grid, end up with wildly low numbers-- below even machine epsilon precision possible
//also beyond likely pig movement distances
if(dist(r,c)<5 && dist(r,c)!=0){
double dval = dist(r,c);
double prob = exp(F2_int + (F2_B*dval));
double probi = exp(F2_int + (F2i_B*dval));

B(r,c) = B1*(I_cells(r,c)*prob)+B2*(C_cells(r,c)*probi);
}

if(dist(r,c+1)<5 && dist(r,c+1)!=0){
double dval1 = dist(r,c+1);
double prob1 = exp(F2_int + (F2_B*dval1));
double probi1 = exp(F2_int + (F2i_B*dval1));
B(r,c+1) = B1*(I_cells(r,c+1)*prob1)+B2*(C_cells(r,c+1)*probi1);
}

}
}

return(B);

}
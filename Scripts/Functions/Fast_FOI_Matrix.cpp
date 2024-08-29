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


///const arma::mat& Imat,
///                            const arma::mat& Cmat,
///                            const double B1,
///                            const double B2



//const arma::mat& M

//initiate output vector
//arma::mat dist(cells,1,fill::zeros);

const int Ninf = id.n_elem;

//arma::mat x_inf = centroids.submat(id,0);
arma::mat X_centroids = centroids.col(0);
arma::mat Y_centroids = centroids.col(1);

arma::mat Xinf = X_centroids.rows(id);
arma::mat Yinf = Y_centroids.rows(id);

//centroids.submat(id,0);
arma::mat X1 = repmat(Xinf, 1, cells);
arma::mat Y1 = repmat(Yinf, 1, cells);

arma::mat X2 = repmat(trans(X_centroids), Ninf, 1);
arma::mat Y2 = repmat(trans(Y_centroids), Ninf, 1);

arma::mat diffX = X2-X1;
arma::mat diffY = Y2-Y1;

arma::mat diffX_2 = diffX % diffX;
arma::mat diffY_2 = diffY % diffY;

arma::mat dist = sqrt(diffX_2+diffY_2);

//arma::mat dist = sqrt((X2-X1)^2+(Y2-Y1)^2);

arma::mat plogit = F2_int + (F2_B*dist);
arma::mat prob = exp(plogit)/(1+exp(plogit));

arma::mat pilogit = F2i_int + (F2i_B*dist);
arma::mat probi = exp(pilogit)/(1+exp(pilogit));

//find where dist=0, set those indices in prob/probi to zero
arma::uvec dist_zero = find(dist == 0);
//arma::ivec distmask(dist.n_rows);

prob(dist_zero) = dist(dist_zero);
probi(dist_zero) = dist(dist_zero);


//B=B1*t(repmat(I[id],cells,1))*prob+B2*t(repmat(C[id],cells,1))*probi
arma::mat B = B1*(Imat%prob)+B2*(Cmat%probi);

//return B;                        

//return Rcpp::List::create(Rcpp::Named("outer")=prob,
//                              Rcpp::Named("inner")=probi);

return(B);

}
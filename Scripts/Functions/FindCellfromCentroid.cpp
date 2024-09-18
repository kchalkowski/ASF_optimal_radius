
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::umat FindCellfromCentroid(const arma::mat& pts,
                            const arma::mat& centroids){

    arma::umat cellnums(pts.n_rows,1);

    for(arma::uword p = 0; p < pts.n_rows; p++) {
        double pt_y = pts(p,1);
        arma::uword y_first;
        arma::uword c_top = 0;
        double cent_y;

    while(c_top<centroids.n_rows && pt_y!=cent_y){
        cent_y = centroids(c_top,1);
        y_first=c_top;
        c_top=c_top+1;
    }

    double cent_y2;
    arma::uword y_last;
    arma::uword c_bottom = centroids.n_rows-1;

    while(c_bottom>=0 && pt_y!=cent_y2){
        cent_y2 = centroids(c_bottom,1);
        y_last=c_bottom;
        c_bottom=c_bottom-1;
    }

        for(arma::uword i = y_first; i < y_last+1; i++){
            if(pts(p,0)==centroids(i,0)){
            cellnums(p,0) = i;
            }
        }
    }


    return(cellnums);

}
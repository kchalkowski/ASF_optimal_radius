
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//Efficient algorithm for finding the cell number of the grid, given the x/y centroid coordinates
//useful for aggregation of disease prev info, etc.
//can be further optimized in future by integrating with RcppParallel.
//row of centroids is synonymous with cell number.

//[[Rcpp::export]]
arma::umat FindCellfromCentroid(const arma::mat& pts,
                            const arma::mat& centroids){

    //initialize output matrix
    arma::umat cellnums(pts.n_rows,1);

    //loop through x/y coordinate points
    for(arma::uword p = 0; p < pts.n_rows; p++) {

        //get py_y, value of y coordinate point p
        double pt_y = pts(p,1);

        //initialize y_first, first row index with y that matches point p
        arma::uword y_first;

        //initialize c_top, for looping through centroids and checking values
        arma::uword c_top = 0;

        //initialize cent_y, value of first matching y coordinate in centroids
        double cent_y;

    //loop through centroids top -> down until pt_y == cent_y
    while(c_top<centroids.n_rows && pt_y!=cent_y){
        cent_y = centroids(c_top,1);
        y_first=c_top;
        c_top=c_top+1;
    }

    //initialize cent_y2, value of last matching y_coordinate in centroids
    double cent_y2;

    //initialize y_last, row index of last matching y_coordinate in centroids
    arma::uword y_last;

    //initialize c_bottom, for looping through centroids and checking values
    arma::uword c_bottom = centroids.n_rows-1;

    //loop through centroids bottom -> up until pt_y == cent_y2
    while(c_bottom>=0 && pt_y!=cent_y2){
        cent_y2 = centroids(c_bottom,1);
        y_last=c_bottom;
        c_bottom=c_bottom-1;
    }

        //loop from rows y_first to y_last to find matching x coordinate
        for(arma::uword i = y_first; i < y_last+1; i++){

            //if pt p c coordinate matches centroids x, that is the row index/cellnumber
            if(pts(p,0)==centroids(i,0)){
            cellnums(p,0) = i; //fyi: need to add 1 if getting output in R to get R index, since cpp indices start from zero
            }
        }
    }

    //return vector of cell indices
    return(cellnums);

}
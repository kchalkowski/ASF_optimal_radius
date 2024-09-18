﻿#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//define func to determine if point is inside a polygon, given the polygon vertices
double calcSidecpp(const arma::vec& pt,
                         const arma::vec& v){

double s = sqrt(pow(pt(0)-v(0),2)+pow(pt(1)-v(1),2));
  return(s);

}

double calcAnglecpp(const arma::vec& pt,
                    const arma::vec& v1,
                    const arma::vec& v2){

    double s1=calcSidecpp(pt,v1);
    double s2=calcSidecpp(pt,v2);
    double s3=calcSidecpp(v1,v2);

    double theta = acos( (pow(s1,2)+pow(s2,2)-pow(s3,2)) / (2*s1*s2) );

    return(theta);
     
}

//[[Rcpp::export]]
//define func to determine if point is inside a polygon, given the polygon vertices
arma::mat determineInsidePoint(const arma::mat& verts,
                               const arma::mat& testpt){
    //verts is vertices of polygon
    //testpt is point testing to see if inside or outside polygon

  //#get all permutations of vertex indices




}


//define convex hull function here
//arma::mat fastConvexHull(const arma::mat& pts){

//}

//define the main function, findZone_fast
//arma::mat findZone_fast(const arma::mat& sl,
 //                       const arma::mat& grid,
 //                       const arma::mat& centroids
 //                       const arma::vec& weeks,
 //                       arma::uvec startloc) {

////Process
///Inputs:
    //sl- out.list$sounderlocs
    //grid- lookup table for x/y coords of grid
    //centroids- matrix of centroid x/y coords for whole grid
    //startloc- index of indicated starting location (for R input, subtract 1 from R index)
            //ie, R index of 20101 = Cpp index of 20100

//1-outside of this function, separate sounderlocs to get either I or C version
//2-loop through weeks (cpp col 2), using input weeks vector
//3-get submatrix of each week, of only positives
    //a-pull out points, sequentially growing vector with each loop
        //i-maybe could preset size-- should be total number of points, plus first case
        //ii-get convex hull of all points *
        //iii-add buffer to convex hull pts *
        //iv-get intersection of chull+buffer and all sounderlocs *
            //-label points in intersection as 'inZone'
        //get cell number for each x/y loc in zone

        
        

//return();

//}
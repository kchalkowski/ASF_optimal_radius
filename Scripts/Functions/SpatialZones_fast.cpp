﻿#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//define func to determine if point is inside a polygon, given the polygon vertices

//pt into calcSidecpp is matrix with rows for each point, two cols for x and y coords, respectively
//vec is matrix with rows for each vertex of a polygon, two cols for x and y coordinates, respectively
double calcSidecpp(const arma::vec& pt,
                         const arma::vec& v){

double s = sqrt(pow(pt(0)-v(0),2)+pow(pt(1)-v(1),2));
  return(s);

}

//pt into calcSidecpp is matrix with rows for each point, two cols for x and y coords, respectively
//v1 and v2 are also x,y points in same format as pt
double calcAnglecpp(const arma::vec& pt,
                    const arma::vec& v1,
                    const arma::vec& v2){

    double s1=calcSidecpp(pt,v1);
    double s2=calcSidecpp(pt,v2);
    double s3=calcSidecpp(v1,v2);

    double theta = acos( (pow(s1,2)+pow(s2,2)-pow(s3,2)) / (2*s1*s2) );

    return(theta);
     
}

//define func to determine if point is inside a polygon, given the polygon vertices
//verts is matrix of polygon vertices, two columns with x and y coords
//pt is x,y coordinates of query point
int determineInsidePoint(arma::mat& verts,
                         const arma::vec& pt){
}

//[[Rcpp::export]]
arma::imat determineInsidePoints(arma::mat& poly,
                          const arma::mat& pts){

    //initialize output matrix
    arma::imat inZone(pts.n_rows,1);

    //loop through all points, determine if located in center of polygon
    for(arma::uword p = 0; p < (pts.n_rows); p++) {


        //arma::vec testpt = pts.row(p);
        arma::vec testpt = arma::conv_to<arma::vec>::from(arma::rowvec(pts.row(p)));
        //arma::vec pt_intcheckv = arma::conv_to<arma::vec>::from(arma::rowvec(pt_intcheck));

        //int inZone_p = determineInsidePoint(poly,testpt);
        //Rcout << "poly" << poly << "\n";
        //Rcout << "testpt" << testpt << "\n";

        //int inZone_p = determineInsidePoint(poly,testpt);
        inZone.row(p) = determineInsidePoint(poly,testpt);

        

    }

    return(inZone);

}

arma::mat DetQuadInterior(const arma::mat& verts,
                          arma::mat pts){
//pts is a matrix where each point is a row. col1=x, col2=y, col3=quadrant assignment, col4=chull assignment


int orientation(arma::mat p,
                arma::mat i,
                arma::mat q){ //p, i, q
int is_counterclockwise;
double det = (q(0,0)-p(0,0))*(i(0,1)-p(0,1))-(q(0,1)-p(0,1))*(i(0,0)-p(0,0));
if(det<0){
is_counterclockwise = 1;
} else{
is_counterclockwise = 0;
}

return(is_counterclockwise);
}

arma::mat JarvisMarch(arma::mat pts){
    //Rcout << "entering jarvis" << "\n";
    uword l = 0;
    arma::uvec chull(pts.n_rows);
    arma::uword chullseq = 0;
    chull(chullseq) = l;
    
    arma::uword p = l;
    arma::uword q = p+1;
    //arma::uword q;
        //Rcout << "before while" << "\n";
        //Rcout << "q" << q << "\n";

    while(q!=l&p<(pts.n_rows-1)){
        q=p+1;
        //Rcout << "q" << q << "\n";
        //Rcout << "q" << q << "\n";
        
            //Rcout << "while top" << "\n";
        for(arma::uword i = 0; i < (pts.n_rows); i++){
        //Rcout << "i" << i << "\n";
                //Rcout << "inside for" << "\n";
                 //Rcout << "pts.row(p)" << pts.row(p) << "\n";
                 //Rcout << "pts.row(i)" << pts.row(i) << "\n";
                 //Rcout << "pts" << pts << "\n";
                 //Rcout << "q" << q << "\n";
                 //Rcout << "pts.row(q)" << pts.row(q) << "\n";

            int i_is_right_of_p = orientation(pts.row(p),pts.row(i),pts.row(q));
            //Rcout << "i_is_right_of_p" << i_is_right_of_p << "\n";
                //Rcout << "after orientation" << "\n";
            if(p!=i && i_is_right_of_p==1){
                q=i;
                //Rcout << "i_is_right_of_p" << i_is_right_of_p << "\n";
            }
      
        }

        chullseq = chullseq + 1;
        chull(chullseq) = q;
        p=q;
        //Rcout << "q" << q << "\n";
        //Rcout << "p" << p << "\n";

    }

    //Rcout << "after while loop" << "\n";


    arma::uvec chull_out=chull(span(0,chullseq));

    arma::mat pts_out = pts.rows(chull_out);

    return(pts_out);

}

//[[Rcpp::export]]
arma::mat FindChullcpp(arma::mat ptsmat){


//[[Rcpp::export]]
arma::mat BufferChullcpp(const arma::mat& pts,
                         arma::mat pt_chull,
                         const int buffer){
//input:
    //takes input from FindChullcpp output
    //pts matrix with nrow number of pts, ncol=4
        //col0-x, col1-y, col2-quadrant number, col3- chull designation
        //chull designation: 1- part of chull, 0-interior of chull

//output:matrix with same form as points, subset to chull points and scaled with input buffer.

//subset chull pts
    //arma::mat pt_chull = pts.rows(find(pts.col(3)==1));
    arma::vec pt_ctr(2);

//get centroid of all points in pts
    pt_ctr(0) = mean(pts.col(0));
    pt_ctr(1) = mean(pts.col(1));

//subtract centerpoint xy from each vertex xy
    pt_chull.col(0) = pt_chull.col(0)-pt_ctr(0);
    pt_chull.col(1) = pt_chull.col(1)-pt_ctr(1);

//for each centered vertex point, scale out
    //initialize alpha scaling parameter
    //formula created using rules for similar triangles
        double C = sqrt(pow(pt_ctr(0),2)+pow(pt_ctr(1),2));
        double alpha = (buffer+C)/C;
    //scale each vertex x and y by alpha, x=alpha*x
        pt_chull.col(0) = pt_chull.col(0)*alpha;
        pt_chull.col(1) = pt_chull.col(1)*alpha;

//shift back to center, re-add centerpoint xy to each vertex xy
    pt_chull.col(0) = pt_chull.col(0)+pt_ctr(0);
    pt_chull.col(1) = pt_chull.col(1)+pt_ctr(1);


return(pt_chull);

}
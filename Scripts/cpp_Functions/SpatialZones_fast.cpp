#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

////////////////////////////////////////
///// Purpose /////////////////////////
////////////////////////////////////////

//Series of functions to enable efficient spatial summarization of points
//Find convex hull from a given matrix of xy coordinates
//Add buffer to convex hull
//Currently used to summarize sounderlocs output
//May later be used in simulation process to more efficiently identify zone around infected individuals

////////////////////////////////////////////////
///// Helper functions /////////////////////////
///////////////////////////////////////////////

//Functions used inside cpp functions, not exported to R

/////////calcSidecpp:
//Calculates the length between two points using pythagorean theorum
//pt is a numeric vector with an x and y coordinate
//vec is a numeric vector with an x and y coordinate
//outputs double s, length between points pt and v
double calcSidecpp(const arma::vec& pt,
                         const arma::vec& v){
double s = sqrt(pow(pt(0)-v(0),2)+pow(pt(1)-v(1),2));
  return(s);
}

/////////calcAnglecpp
//Given three points of a triangle, calculates the angle at point pt using cosine law
//pt, v1 and v2 are all numeric vectors with x, y coordinates
//outputs double theta, angle at point pt.
double calcAnglecpp(const arma::vec& pt,
                    const arma::vec& v1,
                    const arma::vec& v2){

    double s1=calcSidecpp(pt,v1);
    double s2=calcSidecpp(pt,v2);
    double s3=calcSidecpp(v1,v2);

    double theta = acos( (pow(s1,2)+pow(s2,2)-pow(s3,2)) / (2*s1*s2) );

    return(theta);
     
}

/////////determineInsidePoint:
//Determines if point is inside a polygon, given a point and a matrix of polygon vertices. If sum of angles between point and each pair of vertices is not approx equal to 2*pi, the point is outside the set of vertices. 
//verts is numeric matrix of polygon vertices, two columns with x and y coords, ordered counter-clockwise from xmin vertex
//pt is numeric vector x,y coordinates of query point
int determineInsidePoint(arma::mat& verts,
                         const arma::vec& pt){  //initialize vertpairs unsigned integer matrix, for storing pairs of vertices in counter-clockwise order  arma::umat vertpairs(verts.n_rows,2);  //loop through vertices, get pairs  for(arma::uword i = 0; i < (verts.n_rows-1); i++) {    vertpairs.row(i)={i,i+1};    }  //get last pair (last vertex plus first)  vertpairs.row(verts.n_rows-1) = {(verts.n_rows-1),0};  //initialize angles matrix, store angles between point and each pair  arma::mat angles(verts.n_rows,1);   //loop through each pair   for(arma::uword j = 0; j < (verts.n_rows); j++) {     //get coordinates of first vertex in pair     arma::rowvec v1 = verts.row(vertpairs(j,0));     //get coordinates of second vertex in pair     arma::rowvec v2 = verts.row(vertpairs(j,1));     //convert format of vertices from rowvec to vector     arma::vec v1v = arma::conv_to<arma::vec>::from(arma::rowvec(v1));     arma::vec v2v = arma::conv_to<arma::vec>::from(arma::rowvec(v2));     //calculate angle at point pt     angles.row(j) = calcAnglecpp(pt,v1v,v2v);        }  //get sum of angles  double alltheta = accu(angles);  //if angle is near 2*pi (allow for 1e-6 wiggle room), point pt is inside polygon  //pt=1 is inside; pt=0 is not inside  if(alltheta<3.141593*2+0.000001&alltheta>3.141593*2-0.000001){  int out = 1;return(out);  } else {  int out = 0;return(out);}
}

////DetQuadInterior:
//Assigns each point to quadrants around a quadrilateral polygon. Allows for efficient determination of convex hull in FindChullcpp.
//Part of algorithm following method in preprint "CudaChain: Practical GPU-accelerated 2D Convex Hull Algorithm"; arxiv.org/pdf/1508.05488
//verts is numeric matrix with x y coordinates of 4 quadrilateral vertices
//pts is numeric matrix with 4 cols: col1=x, col2=y, col3=quadrant assignment, col4=chull assignment
//output is numeric matrix pts with third col populated with quadrant assignment.
arma::mat DetQuadInterior(const arma::mat& verts,
                          arma::mat pts){
//pts is a matrix where each point is a row. col1=x, col2=y, col3=quadrant assignment, col4=chull assignment
    arma::umat vertpairs(verts.n_rows,2);        for(arma::uword i = 0; i < (verts.n_rows-1); i++) {            vertpairs.row(i)={i,i+1};        }        vertpairs.row(verts.n_rows-1) = {(verts.n_rows-1),0};//for each line made by each combination of vertices, counter-clockwisearma::vec quads = linspace(1, verts.n_rows, verts.n_rows);   for(arma::uword l = 0; l < (verts.n_rows); l++) {        arma::rowvec v1 = verts.row(vertpairs(l,0));        arma::rowvec v2 = verts.row(vertpairs(l,1));        for(arma::uword p = 0; p < (pts.n_rows); p++) {            arma::rowvec pt = pts.row(p);            double det = (v2(0,0)-v1(0,0))*(pt(0,1)-v1(0,1))-(v2(0,1)-v1(0,1))*(pt(0,0)-v1(0,0));            if(det<0){                pts(p,2) = quads(l);            }        }  }return(pts);}

////orientation:
//given three points p, i and q, function identifies whether point i is counterclockwise to the line between points p and q
//uses determinant to do this-- if determinant is negative, point i is counterclockwise (0=on line, positive=clockwise).
//used in Jarvis March algorithm to determine convex hull.
//input p, i and q are all numeric matrices with one row and two cols made of x/y coordinates
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

////JarvisMarch:
//Algorithm for determining the convex hull from a set of points
//input is numeric matrix pts with two cols of x/y coordinates.
//output is numeric matrix, subset of pts x/y coordinates that make up convex hull.
arma::mat JarvisMarch(arma::mat pts){

    //initialize first index, l=0.
    //pts is already organized with first row of points corresponding to xmin,
    //so we know the first index is part of the convex hull.
    uword l = 0;

    //initialize unsigned integer vector chull-- will store indices of convex hull here
    //max length of chull is pts n_rows, will subset later when get actual number of chull pts
    arma::uvec chull(pts.n_rows);

    //chullseq iterator, start at 0
    arma::uword chullseq = 0;

    //first index of convex hull is 0, first point
    chull(chullseq) = l;

    //set point p to l, first point
    arma::uword p = l;

    //set point q to p+1, next point in sequence
    arma::uword q = p+1;

    //loop through points and see if there are any points counter-clockwise to p-q
    //loop through until q=l and before end of p=n_row-1
    //when reach that point, means we have gone around counter-clockwise through all points
    while(q!=l&p<(pts.n_rows-1)){
        //reset q to p+1, needs to happen at start of every loop
        q=p+1;

        //go through all other points to see if any are counter-clockwise to current p-q
        //if any are, set that point to q, and continue through points to see if any are counterclockwise of p-q again.
        for(arma::uword i = 0; i < (pts.n_rows); i++){
            int i_is_right_of_p = orientation(pts.row(p),pts.row(i),pts.row(q));
            //if point i is counter-clockwise of p-q, reset q to i
            if(p!=i && i_is_right_of_p==1){
                q=i;
            }
      
        }

        //after check all i for being right of p-q, reset chullseq iterator
        //will also be used later to get last index of chull for subsetting
        chullseq = chullseq + 1;

        //set next index of chull to q-- this is next index in chull
        chull(chullseq) = q;

        //reset p to q to try for next p-q
        p=q;


    }

    //subset chull, keep only indices found as part of convex hull
    arma::uvec chull_out=chull(span(0,chullseq));

    //use chull_out indices to get x y coordinates of convex hull points
    arma::mat pts_out = pts.rows(chull_out);

    return(pts_out);

}

/////////////////////////////////////////////////////
///// Exported R functions /////////////////////////
///////////////////////////////////////////////////

//These functions are exported for use inside R scripts

//determineInsidePoints
//same as determineInsidePoint, but allows for looping through series of points to identify any that are interior
//poly is numeric matrix with x y coordinates of polygon vertices, ordered counterclockwise from xmin vertex
//pts is numeric matrix with x y coordinates of points to loop through and check if any are interior of polygon
//output is integer matrix with 1s and 0s corresponding to indices of pts; 1=inside poly, 0=not inside poly
//[[Rcpp::export]]
arma::imat determineInsidePoints(arma::mat& poly,
                          const arma::mat& pts){

    //initialize output matrix
    arma::imat inZone(pts.n_rows,1);

    //loop through all points, determine if located in center of polygon
    for(arma::uword p = 0; p < (pts.n_rows); p++) {

        //convert pts p to numeric vector from row vector
        arma::vec testpt = arma::conv_to<arma::vec>::from(arma::rowvec(pts.row(p)));

        //use determineInsidePoint to get 1/0 indicating in/out of polygon
        inZone.row(p) = determineInsidePoint(poly,testpt);
    }

    return(inZone);

}


////FindChullcpp:
//Given set of xy coordinates, finds points forming the convex hull
//Uses algorithm following method in preprint "CudaChain: Practical GPU-accelerated 2D Convex Hull Algorithm"; arxiv.org/pdf/1508.05488
//ptsmat is numeric vector with 2 cols of x/y coordinates
//[[Rcpp::export]]
arma::mat FindChullcpp(arma::mat ptsmat){
    arma::mat pts(ptsmat.n_rows,4);    pts.col(0)=ptsmat.col(0); //add quadrant designation column    pts.col(1)=ptsmat.col(1); //add chull designation column    //initiate new matrix for ordered pts, will be for final output    arma::mat pts_out2(0,2);    arma::mat ptvt(0,pts.n_cols);  if(pts.n_rows>4){  //initialize min-max polygon matrix  arma::mat verts(4,2);  //subset x and y cols  arma::colvec vx = pts.col(0);  arma::colvec vy = pts.col(1);  //get minx/miny/maxx/maxy indices  arma::uword minx = vx.index_min();  arma::uword miny = vy.index_min();  arma::uword maxx = vx.index_max();  arma::uword maxy = vy.index_max();  //get max/min xy coords  verts.row(0) = pts.submat(minx,0,minx,1);  verts.row(1) = pts.submat(miny,0,miny,1);  verts.row(2) = pts.submat(maxx,0,maxx,1);  verts.row(3) = pts.submat(maxy,0,maxy,1);  //make vector of vertex indices  arma::uvec vertices = {minx, miny, maxx, maxy};  //remove vertices from points matrix, add back later into chull  pts.shed_rows(vertices);  //determine quadrants of each point  pts=DetQuadInterior(verts,pts);  //get vertpairs, ordering matrix  arma::umat vertpairs(verts.n_rows,2);        for(arma::uword i = 0; i < (verts.n_rows-1); i++) {            vertpairs.row(i)={i,i+1};        }        vertpairs.row(verts.n_rows-1) = {(verts.n_rows-1),0};  //order pts matrix-- need to go counter-clockwise around verts polygon to determine chull    arma::colvec quadnums = unique(pts.col(2));    //initiate new matrix for ordered pts    arma::mat opts(0,pts.n_cols);    for(arma::uword q = 0; q < (quadnums.size()); q++) {    arma::mat ptq = pts.rows(find(pts.col(2)==quadnums(q)));    //ifquadnums(q)==1, decreasing y THEN increasing x    if(quadnums(q)==1){    arma::uvec quad_order=sort_index(ptq.col(1),"descend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(0),"ascend");    ptq = ptq.rows(quad_order_2);    opts = join_cols(opts,ptq);    }    //ifquadnums(q)==2, increasing X THEN increasing y    if(quadnums(q)==2){    arma::uvec quad_order=sort_index(ptq.col(0),"ascend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(1),"ascend");    ptq = ptq.rows(quad_order_2);    opts = join_cols(opts,ptq);    }    //ifquadnums(q)==3, increasing y THEN decreasing X    if(quadnums(q)==3){    arma::uvec quad_order=sort_index(ptq.col(1),"ascend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(0),"descend");    ptq = ptq.rows(quad_order_2);    opts = join_cols(opts,ptq);    }    //ifquadnums(q)==4, decreasing X THEN decreasing y    if(quadnums(q)==4){    arma::uvec quad_order=sort_index(ptq.col(0),"descend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(1),"descend");    ptq = ptq.rows(quad_order_2);    opts = join_cols(opts,ptq);    }    if(quadnums(q)==0){    opts = join_cols(opts,ptq);    }    }    arma::vec quads = {1, 2, 3, 4};    for(arma::uword q = 0; q < 4; q++) {        //get vertices set of minmax polygon for quadrant q        arma::mat v1 = verts.row(vertpairs(q,0));        arma::mat v2 = verts.row(vertpairs(q,1));          //maximum nrow of qpoly is total number of points in quadrant        arma::uvec quadrant_pts_index = find(opts.col(2)==quads(q)&&opts.col(3)==0);        arma::mat qpoly((quadrant_pts_index.size()+2),2);        qpoly.row(0)=v1;        //set while iterator to zero        int w = 0;        while(any(opts.col(2)==quads(q)&&opts.col(3)==0)){        arma::uvec quadrant_pts = find(opts.col(2)==quads(q)&&opts.col(3)==0);        arma::uword chainlink_index = quadrant_pts(0);             opts(chainlink_index,3) = 1;        arma::mat pt_chull = opts.row(chainlink_index);                 qpoly.row(1+w)=pt_chull.cols(0,1);        qpoly.row(2+w)=v2;        arma::uword lenqpoly = 2+w;        arma::uvec quad_int_check = quadrant_pts(find(quadrant_pts!=chainlink_index));        for(arma::uword p = 0; p < quad_int_check.size(); p++) {            arma::rowvec pt_intcheck = opts.row(quad_int_check(p));    //convert pt_intcheck to vector for use in determineInsidePoint         arma::vec pt_intcheckv = arma::conv_to<arma::vec>::from(arma::rowvec(pt_intcheck));         arma::mat qpoly_subview = qpoly.rows(0,lenqpoly);          int inside = determineInsidePoint(qpoly_subview,pt_intcheckv.subvec(0,1));            if(inside==1){            //is inside            opts(quad_int_check(p),3) = 2;          }    } //for p closing bracket    //bump iterator    w=w+1;} //while closing bracket    } //for each quad closing loop    arma::mat verts_rb(verts.n_rows,4);    verts_rb.submat(0,0,verts.n_rows-1,1) = verts;    verts_rb.col(2)={1,2,3,4};    verts_rb.col(3)={1,1,1,1};    arma::vec quads2 = {1,2,3,4,0};    arma::mat apts=join_cols(opts,verts_rb);    for(arma::uword q = 0; q < (quads2.size()); q++) {    arma::mat ptq = apts.rows(find(apts.col(2)==quads2(q)));    //ifquadnums(q)==1, decreasing y THEN increasing x    if(quads2(q)==1){    arma::uvec quad_order=sort_index(ptq.col(1),"descend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(0),"ascend");    ptq = ptq.rows(quad_order_2);    ptvt = join_cols(ptvt,ptq);    }        //ifquadnums(q)==2, increasing X THEN increasing y    if(quads2(q)==2){    arma::uvec quad_order=sort_index(ptq.col(0),"ascend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(1),"ascend");    ptq = ptq.rows(quad_order_2);    ptvt = join_cols(ptvt,ptq);    }    //ifquadnums(q)==3, increasing y THEN decreasing X    if(quads2(q)==3){    arma::uvec quad_order=sort_index(ptq.col(1),"ascend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(0),"descend");    ptq = ptq.rows(quad_order_2);    ptvt = join_cols(ptvt,ptq);    }    //ifquadnums(q)==4, decreasing X THEN decreasing y    if(quads2(q)==4){    arma::uvec quad_order=sort_index(ptq.col(0),"descend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(1),"descend");    ptq = ptq.rows(quad_order_2);    ptvt = join_cols(ptvt,ptq);    }    if(quads2(q)==0){    ptvt = join_cols(ptvt,ptq);    }    }    arma::mat pts_out = ptvt.cols(0,1);    //1 value indicates unique, 0 is duplicates    arma::uvec uniques(pts_out.n_rows);    uniques(0) = 1;    arma::uword ucheck;    //check to remove any duplicated indices of pts_out    for(arma::uword u = 0; u < (pts_out.n_rows-1); u++) {            ucheck = u+1;            if(pts_out(u,0)==pts_out(ucheck,0) && pts_out(u,1)==pts_out(ucheck,1)){                uniques(ucheck) = 0;            } else{                uniques(ucheck) = 1;            }        //}    }    //probably can remove above unique check thing, not sure if needed    pts_out = pts_out.rows(find(uniques==1));        pts_out = JarvisMarch(pts_out);        pts_out2 = join_cols(pts_out2,pts_out);    } else{    pts_out2 = join_cols(pts_out2,pts.cols(0,1));        }    return(pts_out2);}

////BufferChullcpp:
//Creates buffer around series of points forming convex hull, outputs new vertices forming the buffered polygon
//centers the buffer around set of points
//pts is numeric matrix of x/y coordinates to center buffer around
//pt_chull is numeric matrix of points forming a convex hull around pts
//buffer is a double-- measure of the width of the buffer to draw around pt_chull polygon
//[[Rcpp::export]]
arma::mat BufferChullcpp(const arma::mat& pts,
                         arma::mat pt_chull,
                         const double buffer){
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

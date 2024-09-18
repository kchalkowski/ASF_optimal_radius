#include <RcppArmadillo.h>

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
                         const arma::vec& pt){//Rcout<<"testdip1"<< "\n";  arma::umat vertpairs(verts.n_rows,2);//Rcout<<"testdip1"<< "\n";//Rcout<<vertpairs<< "\n";//Rcout<<verts<< "\n";  for(arma::uword i = 0; i < (verts.n_rows-1); i++) {    vertpairs.row(i)={i,i+1};        }//Rcout<<"testdip1";   vertpairs.row(verts.n_rows-1) = {(verts.n_rows-1),0};//Rcout<<"testdip1";    arma::mat angles(verts.n_rows,1); //Rcout<<"testdip1";    for(arma::uword j = 0; j < (verts.n_rows); j++) {     arma::rowvec v1 = verts.row(vertpairs(j,0));     arma::rowvec v2 = verts.row(vertpairs(j,1));     arma::vec v1v = arma::conv_to<arma::vec>::from(arma::rowvec(v1));     arma::vec v2v = arma::conv_to<arma::vec>::from(arma::rowvec(v2));     angles.row(j) = calcAnglecpp(pt,v1v,v2v);        }//Rcout<<"testdip1";  //#get sum of angles  double alltheta = accu(angles);//Rcout<<"testdip1";    if(alltheta<3.141593*2+0.000001&alltheta>3.141593*2-0.000001){  int out = 1;return(out);  } else {  int out = 0;return(out);}
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
    arma::umat vertpairs(verts.n_rows,2);        for(arma::uword i = 0; i < (verts.n_rows-1); i++) {            vertpairs.row(i)={i,i+1};        }        vertpairs.row(verts.n_rows-1) = {(verts.n_rows-1),0};//for each line made by each combination of vertices, counter-clockwisearma::vec quads = linspace(1, verts.n_rows, verts.n_rows);   for(arma::uword l = 0; l < (verts.n_rows); l++) {        arma::rowvec v1 = verts.row(vertpairs(l,0));        arma::rowvec v2 = verts.row(vertpairs(l,1));        for(arma::uword p = 0; p < (pts.n_rows); p++) {            arma::rowvec pt = pts.row(p);            double det = (v2(0,0)-v1(0,0))*(pt(0,1)-v1(0,1))-(v2(0,1)-v1(0,1))*(pt(0,0)-v1(0,0));            if(det<0){                pts(p,2) = quads(l);            }        }  }return(pts);}

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
    arma::mat pts(ptsmat.n_rows,4);    pts.col(0)=ptsmat.col(0); //add quadrant designation column    pts.col(1)=ptsmat.col(1); //add chull designation column    //initiate new matrix for ordered pts, will be for final output    arma::mat pts_out2(0,2);    arma::mat ptvt(0,pts.n_cols);  if(pts.n_rows>4){  //initialize min-max polygon matrix  arma::mat verts(4,2);  //subset x and y cols  arma::colvec vx = pts.col(0);  arma::colvec vy = pts.col(1);  //get minx/miny/maxx/maxy indices  arma::uword minx = vx.index_min();  arma::uword miny = vy.index_min();  arma::uword maxx = vx.index_max();  arma::uword maxy = vy.index_max();  //Rcout << "minx" << minx << "\n";  //Rcout << "maxx" << maxx << "\n";  //get max/min xy coords  verts.row(0) = pts.submat(minx,0,minx,1);  verts.row(1) = pts.submat(miny,0,miny,1);  verts.row(2) = pts.submat(maxx,0,maxx,1);  verts.row(3) = pts.submat(maxy,0,maxy,1);  //Rcout << "verts" << verts << "\n";  //make vector of vertex indices  arma::uvec vertices = {minx, miny, maxx, maxy};  //Rcout << "vertices" << vertices << "\n"; //0 0 1 3  //remove vertices from points matrix, add back later into chull  pts.shed_rows(vertices);  //Rcout << "pts" << pts << "\n";  //determine quadrants of each point  pts=DetQuadInterior(verts,pts);  //Rcout << "pts" << pts << "\n";  //get vertpairs, ordering matrix  arma::umat vertpairs(verts.n_rows,2);        for(arma::uword i = 0; i < (verts.n_rows-1); i++) {            vertpairs.row(i)={i,i+1};        }        vertpairs.row(verts.n_rows-1) = {(verts.n_rows-1),0};  //order pts matrix-- need to go counter-clockwise around verts polygon to determine chull    arma::colvec quadnums = unique(pts.col(2));    //initiate new matrix for ordered pts    arma::mat opts(0,pts.n_cols);    //Rcout << "test1" << "\n";    //for(q in 1:length(quadnums){    for(arma::uword q = 0; q < (quadnums.size()); q++) {    arma::mat ptq = pts.rows(find(pts.col(2)==quadnums(q)));    //ifquadnums(q)==1, decreasing y THEN increasing x    if(quadnums(q)==1){    arma::uvec quad_order=sort_index(ptq.col(1),"descend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(0),"ascend");    ptq = ptq.rows(quad_order_2);    opts = join_cols(opts,ptq);    }    //ifquadnums(q)==2, increasing X THEN increasing y    if(quadnums(q)==2){    arma::uvec quad_order=sort_index(ptq.col(0),"ascend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(1),"ascend");    ptq = ptq.rows(quad_order_2);    opts = join_cols(opts,ptq);    }    //ifquadnums(q)==3, increasing y THEN decreasing X    if(quadnums(q)==3){    arma::uvec quad_order=sort_index(ptq.col(1),"ascend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(0),"descend");    ptq = ptq.rows(quad_order_2);    opts = join_cols(opts,ptq);    }    //ifquadnums(q)==4, decreasing X THEN decreasing y    if(quadnums(q)==4){    arma::uvec quad_order=sort_index(ptq.col(0),"descend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(1),"descend");    ptq = ptq.rows(quad_order_2);    opts = join_cols(opts,ptq);    }    if(quadnums(q)==0){    opts = join_cols(opts,ptq);    }    }    //Rcout << "test2" << "\n";    arma::vec quads = {1, 2, 3, 4};    for(arma::uword q = 0; q < 4; q++) {        //get vertices set of minmax polygon for quadrant q        arma::mat v1 = verts.row(vertpairs(q,0));        arma::mat v2 = verts.row(vertpairs(q,1));          //maximum nrow of qpoly is total number of points in quadrant        arma::uvec quadrant_pts_index = find(opts.col(2)==quads(q)&&opts.col(3)==0);        arma::mat qpoly((quadrant_pts_index.size()+2),2);        qpoly.row(0)=v1;        //set while iterator to zero        int w = 0;//Rcout << "test3" << "\n";        while(any(opts.col(2)==quads(q)&&opts.col(3)==0)){        arma::uvec quadrant_pts = find(opts.col(2)==quads(q)&&opts.col(3)==0);        arma::uword chainlink_index = quadrant_pts(0);             opts(chainlink_index,3) = 1;        arma::mat pt_chull = opts.row(chainlink_index);                 qpoly.row(1+w)=pt_chull.cols(0,1);        qpoly.row(2+w)=v2;        arma::uword lenqpoly = 2+w;        arma::uvec quad_int_check = quadrant_pts(find(quadrant_pts!=chainlink_index));        for(arma::uword p = 0; p < quad_int_check.size(); p++) {            arma::rowvec pt_intcheck = opts.row(quad_int_check(p));    //Rcout << "test4" << "\n";    //convert pt_intcheck to vector for use in determineInsidePoint         arma::vec pt_intcheckv = arma::conv_to<arma::vec>::from(arma::rowvec(pt_intcheck));         arma::mat qpoly_subview = qpoly.rows(0,lenqpoly);          int inside = determineInsidePoint(qpoly_subview,pt_intcheckv.subvec(0,1));            if(inside==1){            //is inside            opts(quad_int_check(p),3) = 2;          }    //Rcout << "test5" << "\n";} //for p closing bracket    //Rcout << "test6" << "\n";    //bump iterator    w=w+1;} //while closing bracket         //Rcout << "test7" << "\n";    } //for each quad closing loop    arma::mat verts_rb(verts.n_rows,4);    verts_rb.submat(0,0,verts.n_rows-1,1) = verts;    verts_rb.col(2)={1,2,3,4};    verts_rb.col(3)={1,1,1,1};    //Rcout << "test8" << "\n";    arma::vec quads2 = {1,2,3,4,0};    arma::mat apts=join_cols(opts,verts_rb);    //Rcout << "test9" << "\n";    for(arma::uword q = 0; q < (quads2.size()); q++) {    arma::mat ptq = apts.rows(find(apts.col(2)==quads2(q)));    //Rcout << "test10" << "\n";    //ifquadnums(q)==1, decreasing y THEN increasing x    if(quads2(q)==1){    arma::uvec quad_order=sort_index(ptq.col(1),"descend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(0),"ascend");    ptq = ptq.rows(quad_order_2);    ptvt = join_cols(ptvt,ptq);    }        //ifquadnums(q)==2, increasing X THEN increasing y    if(quads2(q)==2){    arma::uvec quad_order=sort_index(ptq.col(0),"ascend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(1),"ascend");    ptq = ptq.rows(quad_order_2);    ptvt = join_cols(ptvt,ptq);    }    //ifquadnums(q)==3, increasing y THEN decreasing X    if(quads2(q)==3){    arma::uvec quad_order=sort_index(ptq.col(1),"ascend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(0),"descend");    ptq = ptq.rows(quad_order_2);    ptvt = join_cols(ptvt,ptq);    }    //ifquadnums(q)==4, decreasing X THEN decreasing y    if(quads2(q)==4){    arma::uvec quad_order=sort_index(ptq.col(0),"descend");    ptq = ptq.rows(quad_order);    arma::uvec quad_order_2=stable_sort_index(ptq.col(1),"descend");    ptq = ptq.rows(quad_order_2);    ptvt = join_cols(ptvt,ptq);    }    if(quads2(q)==0){    ptvt = join_cols(ptvt,ptq);    }    }    //Rcout << "test11" << "\n";    arma::mat pts_out = ptvt.cols(0,1);    //Rcout << "ptvt" << ptvt << "\n";    //Rcout << "pts_out" << pts_out << "\n";    //1 value indicates unique, 0 is duplicates    arma::uvec uniques(pts_out.n_rows);    uniques(0) = 1;    arma::uword ucheck;    //check to remove any duplicated indices of pts_out    for(arma::uword u = 0; u < (pts_out.n_rows-1); u++) {        //for(arma::uword ucheck = u+1; ucheck < (pts_out.n_rows-1); ucheck++){            ucheck = u+1;            if(pts_out(u,0)==pts_out(ucheck,0) && pts_out(u,1)==pts_out(ucheck,1)){                uniques(ucheck) = 0;                //Rcout << "pts_out(u,0)" << pts_out(u,0) << "\n";                //Rcout << "pts_out(u,1)" << pts_out(u,1) << "\n";                //Rcout << "pts_out(ucheck,0)" << pts_out(ucheck,0) << "\n";                //Rcout << "pts_out(ucheck,1)" << pts_out(ucheck,1) << "\n";            } else{                uniques(ucheck) = 1;            }        //}    }    //Rcout << "uniques" << uniques << "\n";    //didn't work, probably remove above unique check thing    pts_out = pts_out.rows(find(uniques==1));        pts_out = JarvisMarch(pts_out);    //Rcout << "test13" << "\n";    //Rcout << "pts_out2" << pts_out2 << "\n";    pts_out2 = join_cols(pts_out2,pts_out);    //Rcout << "test14" << "\n";    } else{    //ptvt.submat(0,0,ptvt.n_rows,ptvt.n_cols) = pts.cols(0,1);    pts_out2 = join_cols(pts_out2,pts.cols(0,1));    //arma::vec chullsign;    //for(arma::uword row = 0; row < (ptvt.n_rows); row++){    //    ptvt(row,3) = 1;    //}        }    return(pts_out2);}

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

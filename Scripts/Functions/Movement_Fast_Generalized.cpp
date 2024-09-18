
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

/////////////////////////////////////////////////////
///////Set up wrapper to run with RcppParallel//////
///////////////////////////////////////////////////

//this wrapper makes the function available to RcppParallel
//the function that you want made parallel goes here in void operator
struct MoveLoop : public Worker {
   // input matrices
   const RMatrix<double> apop;
   const RMatrix<int> apopmat;
   const RMatrix<int> apoplocs;
   const RMatrix<double> acent;
   const int pref;

   // output matrix
   RMatrix<double> outpop;

   // initialize with input and output
   //kind of like calling the function with a header, but with slightly different format



MoveLoop(const NumericMatrix& apop,
           //const IntegerMatrix& apopabund,
           const IntegerMatrix& apopmat,
           const IntegerMatrix& apoplocs,
           const NumericMatrix& acent,
           const int pref,
           NumericMatrix outpop) 
           : apop(apop), apopmat(apopmat), apoplocs(apoplocs), acent(acent), pref(pref), outpop(outpop) {}

//Below conversion funcs are in place because we need them read in as NumericMatrix/IntegerMatrix format
//this is native format for Rcpp, and plays well with RcppParallel
//However, arma formats speeds up simulations considerably
//Reading them in as NumericMatrix/IntegerMatrix and converting them made it work
arma::mat convertpop()
 {
      RMatrix<double> tmp_mat = apop;
  const arma::mat apop2(tmp_mat.begin(), tmp_mat.nrow(), tmp_mat.ncol(), false);
  return apop2;
 }

arma::mat convertcent()
  {
    RMatrix<double> tmp_mat = acent;
    const arma::mat acent2(tmp_mat.begin(), tmp_mat.nrow(), tmp_mat.ncol(), false);
    return acent2;
  }

arma::imat convertapopmat()
  {
    RMatrix<int> tmp_mat = apopmat;
    const arma::imat apopmat2(tmp_mat.begin(), tmp_mat.nrow(), tmp_mat.ncol(), false);
    return apopmat2;
  }

arma::imat convertapoplocs()
  {
    RMatrix<int> tmp_mat = apoplocs;
    const arma::imat apoplocs2(tmp_mat.begin(), tmp_mat.nrow(), tmp_mat.ncol(), false);
    return apoplocs2;
  }

///////////////////////////////////////////////////////////////////////
/////// Parallelized loop through population matrix starts here //////
/////////////////////////////////////////////////////////////////////


void operator()(std::size_t begin, std::size_t end) {
    arma::mat apop3 = convertpop();
    arma::imat apopmat3 = convertapopmat();
    arma::imat apoplocs3 = convertapoplocs();
    arma::mat acent3 = convertcent();
    arma::mat diff(acent.nrow(),1);

    //loop through j rows of pop matrix
    for(std::size_t j = begin; j < end; j++) {

    //get pointers for centroids matrix
    double* cent_x = acent3.colptr(0);
    double* cent_y = acent3.colptr(1);

    //initialize integers for pig movement distance and abundance
    //slightly faster than grabbing the numbers each time
    double pop_j_3=apop3(j,3); //distancem pop[,4]
    int pop_j_0=apop3(j,0); //abundance, pop[,1]

    //loop through each element in centroids to get distance, then take the difference between that and assigned movement distance (pop_j_3)
    //if abundance (pop_j_0) and distance (pop_j_3) are greater than zero
    if(pop_j_3 > 0 & pop_j_0 > 0){
    //initialize distance matrix mask...
    arma::vec mask(diff.n_rows);

    double diffk_0; //initialize diffk_0 double

    //loop through each cell in centroids (acent3)
    for(std::size_t k = 0; k < acent.nrow(); k++) {
    //get distance between current sounder and each cell in centroids, using spatial distance function
    //get difference between that and assigned movement distance
    diffk_0=abs(sqrt(pow((cent_x[k]-apop3(j,4)),2)+pow((cent_y[k]-apop3(j,5)),2))-apop3(j,3));

    //assign difference between assigned and actual distance to diff matrix
    diff(k,0)=diffk_0;

    //find distances closest to assigned movement distance, set mask to isolate those
    if(diffk_0>=0 & diffk_0<=0.4){
    mask[k]=1;
    } else {
    mask[k]=0;    
    }

    } //going through centroids closing bracket

    //get the indices for set of selected cells with distance near the assigned movement distance (within 0.4)
    arma::uvec set = find(mask==1);

    //get size of set
    const int setsize = set.n_elem;

    //if some possible cells to move to... start next selection process
    if(setsize>0){

    //initialize truemin-- selected cellnumber to move to
    arma::uvec truemin;

    ///////////////////////////////////////
    /////// Distance-based movement //////
    /////////////////////////////////////

    //This is default-- if no other options selected (i.e., abundance or rsf), movement is determined by distance alone

    //if pref ==0 (distance-only), randomly sample a cell in set
    if(pref==0){

    truemin = Rcpp::RcppArmadillo::sample(set,1,false);
    }

    ///////////////////////////////////////////////////
    /////// Abundance-based movement preference //////
    /////////////////////////////////////////////////

    //Decide which cells in set to move to, based which cells in set have lowest abundance
    
    if(pref==1){
    //initialize vector for total abundance in each cell in set
    //set corresponds to indices of centroids/grids etc., i.e., cellnumber
    arma::imat abund = apopmat3.rows(set);

    //get minimum abundance value
    int minabundinset = abund.min();
    arma::ivec abundmask(abund.n_rows);

    //Return indices of set which are equal to the minimum abundance value
    for(std::size_t p = 0; p < abund.n_rows; ++p){
        if(abund(p)==minabundinset) abundmask[p]=1; //set poplocs mask to 1 if any loc matches cell in set
    else abundmask[p]=0;
    }

    //then, cellindarma shouldbe wherever mask==1 (ie, wherever the value was equal to the minimum)
    arma::uvec cellindarma = set.elem(find(abundmask==1));
    

    //if there are more than 1 cell with same minimum abundance, choose one at random
    if(cellindarma.size()>1){
    truemin = Rcpp::RcppArmadillo::sample(cellindarma,1,false);

    } else {
    truemin = cellindarma;
    } 

    } //if pref==1

    //////////////////////////////////
    /////// RSF-based movement //////
    ////////////////////////////////

    //base movement on RSF preference only

    if(pref==2){

    //get rsf vals for all cells and subset to cells selected in set
    //double* cent_rsf = acent3.colptr(2);
    //arma::colvec cent_rsf_cv = acent3.col(2);

    //convert to arma vec
    arma::vec cent_rsf = arma::conv_to<arma::vec>::from(arma::colvec(acent3.col(2)));

    //subset to get rsf_vals of cells in selected set
    arma::vec rsf_vals = cent_rsf(set);

    //use rsf_vals as probabilities for which cell to move to
    //truemin = Rcpp::RcppArmadillo::sample(cellindarma,1,false);
    truemin = Rcpp::RcppArmadillo::sample(set,1,0,rsf_vals);
    
    }

    //////////////////////////////////////////////////////////////////
    /////// Assign chosen location after movement pref options //////
    ////////////////////////////////////////////////////////////////

    //set location to selected cell in set with minimum abundance
    outpop(j,0)=truemin[0]+1; //+1 is to get appropriate index

    } else{ //else to 'if any cells in set'

    //there should always be a possible cell to move to (unless barriers introduced in model)
    //currently, if no cells in set, should generate error
    //this will output unrealistic location number that can be used in R script to generate error
    outpop(j,0)=acent.nrow()+1000;

    }


    } else{ //if movement distance is zero, or abundance is zero

    outpop(j,0)=apoplocs(j,0);


    } //else if movement distance/abundance is zero closing loop

    } //worker for loop

} //void operator closing loop

}; //worker closing loop


////////////////////////////////////////////
///////Call worker to run in parallel//////
//////////////////////////////////////////

//So, above, Movement_worker derives from RcppParallel::Worker
//this is required for function objects passed to parallelFor
//Now that the worker is described above, can call the Movement_worker worker we defined

//Here's a function that calls the SquareRoot worker defined above
// [[Rcpp::export]]
NumericMatrix parallelMovementRcpp_portion(const NumericMatrix& apop,
           const IntegerMatrix& apopmat,
           const IntegerMatrix& apoplocs,
           const NumericMatrix acent,
           const int pref){

  // allocate the output matrix
  //essentially just defining the size of the output
  NumericMatrix outpop(apop.nrow(), 1);
  
  // x is input matrix defined above
  //outpop is output matrix defined above
    MoveLoop moveloop(apop,apopmat,apoplocs,acent,pref,outpop);

  // call parallelFor to do the work
  // starting from 0 to length of apop, run the squareRoot function defined above in the worker
  //parallelFor(0, x.length(), squareRoot);
    parallelFor(0,apop.nrow(),moveloop);
  
  // return the output matrix
  return outpop;
//}
}

///////////////////////////
///////Run the thing//////
/////////////////////////

//[[Rcpp::export]]
//define the main function, MovementRcpp
NumericMatrix MovementRcppParallel(NumericMatrix& apop,
                       IntegerMatrix& apopmat,
                       IntegerMatrix& apoplocs,
					   NumericMatrix& acent,
                       const int pref) {

//define the main function, MovementRcpp

Rcpp::NumericMatrix popout = apop;

//set present locations to previous locations
popout(_,6)=popout(_,2);

//get new locations using parallel movement function
popout(_,2)=parallelMovementRcpp_portion(apop,apopmat,apoplocs,acent,pref);


return popout;

}
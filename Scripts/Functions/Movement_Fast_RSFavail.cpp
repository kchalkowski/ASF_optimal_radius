
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
   const RMatrix<int> RSF_mat0;
   const RMatrix<double> RSF_mat;
   const int pref;
   // output matrix
   RMatrix<double> outpop;

   // initialize with input and output
   //kind of like calling the function with a header, but with slightly different format



MoveLoop(const NumericMatrix& apop,
           const IntegerMatrix& apopmat,
           const IntegerMatrix& apoplocs,
           const NumericMatrix& acent,
           const IntegerMatrix& RSF_mat0,
           const NumericMatrix& RSF_mat,
           const int pref,
           NumericMatrix outpop) 
           : apop(apop), apopmat(apopmat), apoplocs(apoplocs), acent(acent), RSF_mat0(RSF_mat0), RSF_mat(RSF_mat), pref(pref), outpop(outpop) {}

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

arma::imat convertRSF_mat0()
  {
    RMatrix<int> tmp_mat = RSF_mat0;
    const arma::imat RSF_mat02(tmp_mat.begin(), tmp_mat.nrow(), tmp_mat.ncol(), false);
    return RSF_mat02;
  }

arma::mat convertRSF_mat()
  {
    RMatrix<double> tmp_mat = RSF_mat;
    const arma::mat RSF_mat2(tmp_mat.begin(), tmp_mat.nrow(), tmp_mat.ncol(), false);
    return RSF_mat2;
  }

///////////////////////////////////////////////////////////////////////
/////// Parallelized loop through population matrix starts here //////
/////////////////////////////////////////////////////////////////////


void operator()(std::size_t begin, std::size_t end) {
    arma::mat apop3 = convertpop();
    arma::imat apopmat3 = convertapopmat();
    arma::imat apoplocs3 = convertapoplocs();
    arma::mat acent3 = convertcent();
    arma::imat RSF_mat03 = convertRSF_mat0();
    arma::mat RSF_mat3 = convertRSF_mat();
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

    //if pref ==0 (distance-only), randomly sample a cell in set
    //if(pref==0){

    //truemin = Rcpp::RcppArmadillo::sample(set,1,false);
    //}

    ///////////////////////////////////////////////////
    /////// Abundance-based movement preference //////
    /////////////////////////////////////////////////

    //Decide which cells in set to move to, based which cells in set have lowest abundance
    
    //if(pref==1){
    //initialize vector for total abundance in each cell in set
    //set corresponds to indices of centroids/grids etc., i.e., cellnumber
    //arma::imat abund = apopmat3.rows(set);

    //get minimum abundance value
    //int minabundinset = abund.min();
    //arma::ivec abundmask(abund.n_rows);

    //Return indices of set which are equal to the minimum abundance value
    //for(std::size_t p = 0; p < abund.n_rows; ++p){
    //    if(abund(p)==minabundinset) abundmask[p]=1; //set poplocs mask to 1 if any loc matches cell in set
    //else abundmask[p]=0;
    //}

    //then, cellindarma shouldbe wherever mask==1 (ie, wherever the value was equal to the minimum)
    //arma::uvec cellindarma = set.elem(find(abundmask==1));
    

    //if there are more than 1 cell with same minimum abundance, choose one at random
    //if(cellindarma.size()>1){
    //truemin = Rcpp::RcppArmadillo::sample(cellindarma,1,false);

    //} else {
    //truemin = cellindarma;
    //} 

    //} //if pref==1

    ////////////////////////////////////////////////////
    /////// RSF-probabilities from grid movement //////
    //////////////////////////////////////////////////

    //base movement on RSF preference
    //this option utilizes a grid where RSF probabilities are coded into landscape
    //faster option if RSF preference is consistent and doesn't change much with availability

    //if(pref==2){

    //get rsf vals for all cells and subset to cells selected in set
    //double* cent_rsf = acent3.colptr(2);
    //arma::colvec cent_rsf_cv = acent3.col(2);

    //convert rsf val col from centroids to arma vec
    //arma::vec cent_rsf = arma::conv_to<arma::vec>::from(arma::colvec(acent3.col(2)));

    //subset to get rsf_vals of cells in selected set
    //arma::vec rsf_vals = cent_rsf(set);

    //use rsf_vals as probabilities for which cell to move to
    //truemin = Rcpp::RcppArmadillo::sample(cellindarma,1,false);
    //truemin = Rcpp::RcppArmadillo::sample(set,1,0,rsf_vals);
    
    //}

    ////////////////////////////////////////////////
    /////// RSF-availability matrix movement //////
    //////////////////////////////////////////////

    //base movement on RSF preference
    //this option utilizes a grid where RSF probabilities are coded into landscape
    //faster option if RSF preference is consistent and doesn't change much with availability

    //RSF_mat0 input:
        //col for each lc type,
        //first row is lc integer
        //all subsequent rows are availability combinations

    //RSF_mat input
        //col for each lc type
        //first row is lc integer
        //all subsequent rows are rsf vals

    if(pref==3){

    //Rcout << "entering if pref==3" << "\n";

    //convert rsf val col from centroids to arma vec
    arma::vec cent_rsf = arma::conv_to<arma::vec>::from(arma::colvec(acent3.col(2)));

    //subset to get rsf_vals of cells in selected set
    arma::vec rsf_vals = cent_rsf(set);

    //get unique rsf values
    arma::vec avail = unique(rsf_vals);
    //Rcout << "avail: " << avail << "\n";

    //dummy code rsf values
    //init dummy code matrix: one row for lc vals, 2nd row for dummy codes of lc avail
    arma::imat avail_mat(2,RSF_mat03.n_cols,fill::zeros);

    //first row of rsf_mat should be integers of lc types
    avail_mat.row(0) = RSF_mat03.row(0);
    
    
    //loop through avail, if hit matches with lc type row, second row gets 1
    for(std::size_t ai = 0; ai < avail.n_elem; ++ai){

        //loop through lc indices
        for(std::size_t lc = 0; lc < RSF_mat03.n_cols; ++lc){
        if(avail(ai)==RSF_mat03(0,lc)) avail_mat(1,lc)=1;
        //else avail_mat(1,lc)=0;
        } //for lc indices
    
    } //for elem in avail
    //Rcout << "avail_mat: " << avail_mat << "\n";

    //get sum of avail_mat dummy codes
    int avail_sums = sum(avail_mat.row(1));
    //Rcout << "avail_sums: " << avail_sums << "\n";

    //get sum of all rsf dummy codes
    //sum rows of rsf_mat_dummy (dim=1)
    arma::icolvec rsf_sums = sum(RSF_mat03,1);
    //Rcout << "rsf_sums: " << rsf_sums << "\n";

    //find indices of rsf_sums that match avail_sums
    //these indices will be set to further subset and match
    //think this will be faster than testing by each index right out of the gate
    arma::ivec sums_mask(rsf_sums.n_elem); //initialize mask vector
    for(std::size_t su = 0; su < rsf_sums.n_elem; ++su){

        //if sums match, go through rest of index matching
        if(avail_sums==rsf_sums(su)){
        //loop through each column of rsf_mat row su, look for match with avail
        arma::ivec match_mask(RSF_mat03.n_cols);

        for(std::size_t lc = 0; lc < RSF_mat03.n_cols; ++lc){
        if(RSF_mat03(su,lc)==avail_mat(1,lc)) match_mask(lc) = 1;
        else match_mask(lc) = 0;
        } //for cols in rsf mat

        //all function tests if all elements in object are non-zero
        if(all(match_mask)) sums_mask(su) = 1;
        else sums_mask(su) = 0;

        } else{ //else if sums don't match
        sums_mask(su) = 0;
        } //else close
        
    } //for rsf sums

    //Rcout << "sums_mask: " << sums_mask << "\n";

    //by end of loop, should have a sums_mask where exactly one element is 1
        //if all zero, no matches
        //if multiple 1, multiple matches
        //should configure stop here for these conditions

    //use sums_mask to subset RSF probabilities for each lc
    //find(sums_mask==1)

    arma::uvec indices(2);
    indices(0) = 0;
    arma::uvec sel = find(sums_mask==1);
    indices(1) = sel(0);

    //arma::uvec indices = {0,find(sums_mask==1)};
    arma::mat RSF_group = RSF_mat3.rows(indices);

    //use subset rsf probs to reclassify lc's in set to probabilities
    //rsf_vals
    //loop through cols of RSF_group
        //if avail_mat(col)==1 //is in group
            //rsf_vals(find(rsf_vals==RSF_group(0,col))) = RSF_group(1,col);
    for(std::size_t col = 0; col < RSF_group.n_cols; ++col){
        if(avail_mat(1,col)==1){
            arma::uvec reclass_ind= find(rsf_vals == RSF_group(0,col));

            //create vector of length reclass_ind.n_elem with RSF_group(1,col) value
            arma::vec RSF_probs(reclass_ind.n_elem);
            for(std::size_t r = 0; r < reclass_ind.n_elem; ++r){
                RSF_probs(r) = RSF_group(1,col);
            } //for reclass

            //use vec of repeated RSF prob values, place to matching rsf vals found
            rsf_vals(reclass_ind) = RSF_probs;

        } //if availmat =1
    } //for vals in rsf_group

    //result should be rsf_vals with probabilities, not whole number lc integers
    //add stop condition for this

    //do weighted selection using probabilities as before
    //use rsf_vals as probabilities for which cell to move to
    truemin = Rcpp::RcppArmadillo::sample(set,1,0,rsf_vals);
    
    } //if close

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

    } //close else


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
NumericMatrix parallelMovement_RSFavail(const NumericMatrix& apop,
           const IntegerMatrix& apopmat,
           const IntegerMatrix& apoplocs,
           const NumericMatrix acent,
           const IntegerMatrix& RSF_mat0,
           const NumericMatrix& RSF_mat,
           const int pref){

  // allocate the output matrix
  //essentially just defining the size of the output
  NumericMatrix outpop(apop.nrow(), 1);
  
  // x is input matrix defined above
  //outpop is output matrix defined above
    MoveLoop moveloop(apop,apopmat,apoplocs,acent,RSF_mat0,RSF_mat,pref,outpop);

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
                       IntegerMatrix& RSF_mat0,
                       NumericMatrix& RSF_mat,
                       const int pref) {

//define the main function, MovementRcpp

Rcpp::NumericMatrix popout = apop;

//set present locations to previous locations
popout(_,6)=popout(_,2);

//get new locations using parallel movement function
popout(_,2)=parallelMovement_RSFavail(apop,apopmat,apoplocs,acent,RSF_mat0,RSF_mat,pref);


return popout;

}
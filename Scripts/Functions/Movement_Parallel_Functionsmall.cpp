// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

/////////////////////////////////////////////////////
///////Set up wrapper to run with RcppParallel//////
///////////////////////////////////////////////////

//this wrapper makes the function available to RcppParallel sort of
//the function that you want made parallel goes here in void operator
struct MoveLoop : public Worker {
   // input matrices
   const RMatrix<double> apop;
   const RMatrix<int> apopabund;
   const RMatrix<int> apoplocs;
   const RMatrix<double> acent;

   // output matrix
   RMatrix<double> outpop;

   // initialize with input and output
   //kind of like calling the function with a header, but with slightly different format



MoveLoop(const NumericMatrix& apop,
           const IntegerMatrix& apopabund,
           const IntegerMatrix& apoplocs,
           const NumericMatrix& acent,
           NumericMatrix outpop) 
           : apop(apop), apopabund(apopabund), apoplocs(apoplocs), acent(acent), outpop(outpop) {}

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

arma::imat convertapopabund()
  {
    RMatrix<int> tmp_mat = apopabund;
    const arma::imat apopabund2(tmp_mat.begin(), tmp_mat.nrow(), tmp_mat.ncol(), false);
    return apopabund2;
  }

arma::imat convertapoplocs()
  {
    RMatrix<int> tmp_mat = apoplocs;
    const arma::imat apoplocs2(tmp_mat.begin(), tmp_mat.nrow(), tmp_mat.ncol(), false);
    return apoplocs2;
  }

///////////////////////////////////////////
///////Movement function begins here//////
/////////////////////////////////////////


void operator()(std::size_t begin, std::size_t end) {
    // rows we will operate on
        arma::mat apop3 = convertpop();
        arma::mat acent3 = convertcent();
        arma::imat apopabund3 = convertapopabund();
        arma::imat apoplocs3 = convertapoplocs();
        arma::mat diff(acent.nrow(),1);

///////Starting loop through piggy pop matrix//////

//loop through j rows of pop matrix
for(std::size_t j = begin; j < end; j++) {

//get pointers for centroids matrix
double* cent_x = acent3.colptr(0);
double* cent_y = acent3.colptr(1);

///////Getting distances between sounder and each cell in grid//////

//initialize integers for pig movement distance and abundance
//slightly faster than grabbing the numbers each time
double pop_j_3=apop3(j,3); //distancem pop[,4]
int pop_j_0=apop3(j,0); //abundance, pop[,1]

//loop through each element in centroids to get distance, then take the difference

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

//find diffk_0 closest to assigned movement distance, set mask to isolate those
if(diffk_0>=0 & diffk_0<=0.4){
mask[k]=1;
} else {
mask[k]=0;    
}

} //going through centroids closing bracket

//get the indices for set of locations with distance near the assigned movement distance
arma::uvec set = find(mask==1);

//get size of set
const int setsize = set.n_elem;

///////Decide which cells in set to move to, based on abundance//////

//if some possible cells to move to...
if(setsize>0){

//initialize vector for total abundance in each cell in set
arma::ivec abund(setsize);

//make mask vector 'setmask' for grabbing indices of set in loop below
arma::ivec setmask(apoplocs3.n_rows);

//loop through each cell in set, get total abundance in each cell in set
for(std::size_t s = 0; s < setsize; ++s) {

//get index of current cell in set
int set_s = set(s);

//run through each item in apoplocs, get locations of each cell in set
//need to do this before get abundances, because there is possibility that
//multiple sounders could be in a single cell
//this means that can't just find cell with min abundance first and then get that location
//need to find all sounders in each cell in set, then sum abundance
for(std::size_t p = 0; p < apoplocs3.n_rows; ++p){

//set poplocs mask to 1 if any loc matches cell in set
if(apoplocs3(p)==set_s) setmask[s]=1; 
else setmask[s]=0;

}

//get abundances of cells in set
imat abundinset_s = apopabund3.rows(find(setmask==1));
abund(s)=sum(abundinset_s.col(0));
//Rcout << "abund " << abund << "\n";
}

//find cell in set with minimum abundance
//get the index of cell(s) in set equal to the minimum value
//this is subset of set, so the index is actually of centroids.. aka cell ID

//get minimum abundance value
int minabundinset = abund.min();
arma::ivec abundmask(abund.n_rows);

//Return indices of set which are equal to the minimum abundance value
for(std::size_t p = 0; p < abund.n_rows; ++p){
  if(abund(p)==minabundinset) abundmask[p]=1; //set poplocs mask to 1 if any loc matches cell in set
  else abundmask[p]=0;
}

//then, cellindarma shouldbe wherever mask==1 (ie, wherever the value was equal to the minimum)
//imat abundinset_s = apopabund3.rows(find(setmask==1));
arma::uvec cellindarma = set.elem(find(abundmask==1));

//initialize selected minimum value vector
arma::uvec truemin; 

//if there are more than 1 cell with same minimum abundance, choose one at random
if(cellindarma.size()>1){
truemin = Rcpp::RcppArmadillo::sample(cellindarma,1,false);

} else {
truemin = cellindarma;
}

///////Assign location for sounder//////

//set location to selected cell in set with minimum abundance
outpop(j,0)=truemin[0]+1; //+1 is to get appropriate index

} else{ //this is the else to the if statement 'if any cells in set'

//there should always be a possible cell to move to (unless barriers introduced in model)
//currently, if no cells in set, should generate error
//this will output unrealistic location number that will be used in R script to generate error
outpop(j,0)=acent.nrow()+1000;


}


} else{ //if movement distance is zero, or abundance is zero

outpop(j,0)=apoplocs(j,0);


}

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
           const IntegerMatrix& apopabund,
           const IntegerMatrix& apoplocs,
           const NumericMatrix acent){

  // allocate the output matrix
  //essentially just defining the size of the output
  NumericMatrix outpop(apop.nrow(), 1);
  
  
  // x is input matrix defined above
  //outpop is output matrix defined above
    MoveLoop moveloop(apop,apopabund,apoplocs,acent,outpop);

  // call parallelFor to do the work
  // starting from 0 to length of apop, run the squareRoot function defined above in the worker
  //parallelFor(0, x.length(), squareRoot);
    parallelFor(0,apop.nrow(), moveloop);
  
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
                       IntegerMatrix& apopabund,
                       IntegerMatrix& apoplocs,
					   NumericMatrix& acent) {
                                            
Rcpp::NumericMatrix popout = apop;

//set present locations to previous locations
popout(_,6)=popout(_,2);

//get new locations using parallel movement function
popout(_,2)=parallelMovementRcpp_portion(apop,apopabund,apoplocs,acent);

return popout;

}
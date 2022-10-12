// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

//this wrapper makes the function available to RcppParallel sort of
//the function that you want made parallel goes here in void operator
//this sample function takes the square root of its input and writes it into its output
//struct SquareRoot : public Worker
struct MoveLoop : public Worker {
   // source matrix
   //const RMatrix<double> input;
   const RMatrix<double> apop;
   const RMatrix<int> apopabund;
   const RMatrix<int> apoplocs;
   const RMatrix<double> acent;

   // destination matrix
   //RMatrix<double> output;
   RMatrix<double> outpop;

   // initialize with source and destination
   //kind of like calling the function with a header, but with slightly different format
   //SquareRoot(const NumericMatrix input, NumericMatrix output) 
   //  : input(input), output(output) {}

  MoveLoop(const NumericMatrix& apop,
           const IntegerMatrix& apopabund,
           const IntegerMatrix& apoplocs,
           const NumericMatrix& acent,
           NumericMatrix outpop) 
           : apop(apop), apopabund(apopabund), apoplocs(apoplocs), acent(acent), outpop(outpop) {}

   // take the square root of the range of elements requested
   //this is where the function goes
   //need to use format for initialization-- void operator()(input stuff like in function) {function goes here}
   //tells it to go through the specified range size_t begin and end...
   //void operator()(std::size_t begin, std::size_t end) {
   //   std::transform(input.begin() + begin, 
   //                  input.begin() + end, 
   //                  output.begin() + begin, 
   //                  ::sqrt);
   //}

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

void operator()(std::size_t begin, std::size_t end) {
    // rows we will operate on
        arma::mat apop3 = convertpop();
        arma::mat acent3 = convertcent();
        arma::imat apopabund3 = convertapopabund();
        arma::imat apoplocs3 = convertapoplocs();
        arma::mat diff(acent.nrow(),1);
        //arma::uvec truemin;

//for( j = 0; j < np; ++j) {
for(std::size_t j = begin; j < end; j++) {

int pop_j_3=apop3(j,3);

if(pop_j_3 > 0){
//loop through each element in centroids to get distance, then take the difference
arma::vec mask(diff.n_rows);
double diffk_0;
for(std::size_t k = 0; k < acent.nrow(); k++) {
diffk_0=sqrt(pow((acent3(k,0)-apop3(j,4)),2)+pow((acent3(k,1)-apop3(j,5)),2))-apop3(j,3);
diff(k,0)=diffk_0;
if(diffk_0>=0 & diffk_0<=0.1){
mask[k]=1;
} else {
mask[k]=0;    
}
}

arma::uvec set = find(mask==1);

const int setsize = set.n_elem;

if(setsize>0){

//initialize vector for total abundance in each cell in set
arma::ivec abund(setsize); 
arma::ivec setmask(apoplocs3.n_rows);
//loop through each cell in set, get total abundance in each cell in set
for(std::size_t s = 0; s < setsize; ++s) {
int set_s = set(s);
//abund(s) = ArmaSubsetPoplocs(apoplocs3, apopabund3, set_s);
for(std::size_t p = 0; p < apoplocs3.n_rows; ++p){
if(apoplocs3(s)==set_s) setmask[s]=1;
else setmask[s]=0;
}

imat abundinset_s = apopabund3.rows(find(setmask==1));
abund(s)=sum(abundinset_s.col(0));

}


//get the index of cell(s) in set equal to the minimum value
//this is subset of set, so the index is actually of centroids.. aka cell ID
arma::uvec cellindarma=set.row(abund.index_min());

arma::uvec truemin;

///trying to get sample to work with arma type
if(cellindarma.size()>1){

truemin = Rcpp::RcppArmadillo::sample(cellindarma,1,false);

} else {
truemin = cellindarma;
}

outpop(j,0)=truemin[0]+1;

} else{

outpop(j,0)=apoplocs(j,0);

}


//store new locations in pop matrix
 //change cell ID of new location (ie index of centroids stored in pop)
//apop(j,4)=acent(truemin[0],0); //get x value of new location
//apop(j,5)=acent(truemin[0],1); //get y value of new location

//}

} else{

outpop(j,0)=apoplocs(j,0);

}

}

} //if any distance assigned greater than zero

};



//So, above, Movement_worker derives from RcppParallel::Worker
//this is required for function objects passed to parallelFor
//Now that the worker is described above, can call the Movement_worker worker we defined


//Here's a function that calls the SquareRoot worker defined above
//NumericMatrix parallelMatrixSqrt(NumericMatrix x) { //create a new function wrapper that calls the part you want run in parallel
// [[Rcpp::export]]
NumericMatrix parallelMovementRcpp_portion(const NumericMatrix& apop,
           const IntegerMatrix& apopabund,
           const IntegerMatrix& apoplocs,
           const NumericMatrix acent){

  // allocate the output matrix
  //essentially just defining the size of the output
  NumericMatrix outpop(apop.nrow(), 1);
  
  // SquareRoot functor (pass input and output matrixes)
  // x is input matrix defined above
  //output is output matrix defined above
  //SquareRoot squareRoot(x, output);
    MoveLoop moveloop(apop,apopabund,apoplocs,acent,outpop);
   //apop(apop), apopabund(apopabund), apoplocs(apoplocs), acent(acent), cells(cells), np(np), popoutput(popoutput) {}
  // call parallelFor to do the work
  // starting from 0 to length of x, run the squareRoot function defined above in the worker
  //parallelFor(0, x.length(), squareRoot);
    parallelFor(0,apop.nrow(), moveloop);
  
  // return the output matrix
  return outpop;
//}
}

//[[Rcpp::export]]
//define the main function, MovementRcpp
NumericMatrix MovementRcppParallel(NumericMatrix& apop,
                       IntegerMatrix& apopabund,
                       IntegerMatrix& apoplocs,
					   NumericMatrix& acent) {
                                            
Rcpp::NumericMatrix popout = apop;
//set present locations to previous locations
popout(_,6)=popout(_,2);

popout(_,2)=parallelMovementRcpp_portion(apop,apopabund,apoplocs,acent);
//popout(_,4)=acent(popout(_,2),0);
//popout(_,5)=acent(popout(_,2),1);

return popout;

}
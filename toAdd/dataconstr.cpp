// [[Rcpp::depends(RcppArmadillo)]]

#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include <RcppArmadillo.h>
using namespace Rcpp;


// function to reconstruct data s.t. constraints at interim
// [[Rcpp::export]]
List reconstruct_data(
    const arma::mat&     bData,
    const arma::rowvec&  minSSEnr_i,
    const arma::rowvec&  maxSSEnr_i,
    const arma::irowvec& active,
    const int&           K0,
    const int&           nStart,
    const int&           nMax,
    const int&           numAdd_i
) {

  // initialize result
  arma::mat datMat(K0, 2);    // stores (y, n-y)
  
  // find the number of subjects needed to meet requirements on minimum enrollment
  // without exceeding requirements on maximum enrollment;
  arma::rowvec nAcc(K0,arma::fill::zeros);
  arma::rowvec nMinCritMet(K0,arma::fill::zeros);
  arma::rowvec nMaxCritMet(K0,arma::fill::zeros);
  for ( int k = 0; k < K0; k++ ) {
    if ( minSSEnr_i(k) == 0) { 
      nMinCritMet(k) = 1; 
    }
  }
  double duration = 0.0;
  int idx2 = 0;
  for (int n = nStart; n < nMax; n++) {
    int b = bData(n, 1);
    int y = bData(n, 0);
    
    if (active(b) == 1 and nMaxCritMet(b) == 0) {
      nAcc(b) +=1;
      if ( nAcc(b) >= minSSEnr_i(b) ) { 
        nMinCritMet(b) = 1; 
      }
      if ( nAcc(b) >= maxSSEnr_i(b) ) { 
        nMaxCritMet(b) = 1; 
      }
      datMat(b,1-y)   += 1;            // if y = 1, adds to first column of datMat, otherwise adds to second column
      duration         = bData(n,4);   // = ft						
      
      if ( (sum(nAcc) >= numAdd_i) and (sum(nMinCritMet) == sum(active)) ) {
        n = nMax + 100;
      }
    }	
    idx2++;
  }
  
  // store number of added observations for each basket
  arma::vec nVec_add = arma::sum(datMat, 1);   
  
  List res = List::create(
    _["datMat"]   = datMat,     // collapsed data matrix
    _["idx2"]     = idx2,       // to update total number of subjects
    _["nVec_add"] = nVec_add,   // number of additional patients in this analysis
    _["duration"] = duration    // max duration (ft)
  );
  return res;
}
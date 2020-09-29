// [[Rcpp::depends(RcppArmadillo)]]

#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include <RcppArmadillo.h>

using namespace Rcpp;

//' Compute number of models
//' 
//' Given a basket size and maximal number of distinct
//' response rates, compute the number of possible models
//' 
//' @param K positive integer giving number of baskets
//' @param P positive integer giving maximal number of distinct rates
//' 
//' @return integer giving number of possible models
//' @noRd
// [[Rcpp::export]]
int numModels_cpp( int const& K, int const& P ) {
  double res = 0;
  for ( int p = 1; p <= P; p++ ) {
    double temp = 0;
    for ( int j = 0; j <= P; j++ ) {
      temp += pow(-1, p-j) * pow(j, K) * R::choose(p, j);
    }
    res += pow( R::gammafn(p + 1.0) , -1) * temp;
  }
  return( (int) res);
}




//' Model matrix of data given model
//' 
//' Get the (collapsed) model data specified by a particular partition (model)
//' 
//' @param datMat     matrix of \code{(y, n-y)} data
//' @param partition  vector of indices giving how to partition the data
//' 
//' @keywords internal
//' @noRd
arma::mat modelMatrix(
    arma::mat const& datMat,
    arma::vec const& partition
) {
  
  // max index = # of elements of list
  int m = partition.max() + 1;                
  
  // initialize matrix for results
  arma::mat partMat = arma::mat( m, datMat.n_cols, arma::fill::zeros );
  
  for ( int i = 0; i < m; i++ ) {
    // ith row of result is column sum of datMat corresponding partition vector being equal to i (elements in ith set of partition)
    partMat.row(i) = arma::sum( datMat.rows( find( partition == i ) ) , 0 );
  }
  return partMat;
}




//' Bayesian model averaging
//' 
//' Computes posterior model probabilities and Bayes model averaged survival function \math{P(\pi_k > x | D)}
//' 
//' @param x              scalar between 0 and 1 giving cutoff for effect size
//' @param datMat         matrix of data \code{(y, n - y)}
//' @param partitionMat   martrix giving how to partition the data for each model
//' @param mu0            scalar giving prior mean for beta prior
//' @param phi0           scalar giving prior dispersion for beta prior
//' @param logchoose      scalar giving \code{sum(lchoose(n, y))}
//' @param P              integer scalar giving number of models
//' @param logModelPriors vector of length P giving the normalized priors for each model
//' 
//' @return a list giving \code{[[posteriorSurvivalFunction], [posteriorModelProbs]]}
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
List bma_cpp (
    double    const& x,
    arma::mat const& datMat,
    arma::mat const& partitionMat,
    double    const& mu0,
    double    const& phi0,
    double    const& logchoose,
    arma::vec const& logModelPriors
) {
  
  // Initialize a vector to store posterior model probabilities and
  // a double to store the bayes average survival function
  arma::vec postModelProbs = logModelPriors + logchoose;                         // initialize vector of posterior model probs:P(M_j | D) to be the prior probs + sum(log(choose(n,y)))
  arma::vec postSurvProbs  = arma::vec(partitionMat.n_cols, arma::fill::zeros);  // initialize vector to store Pr(pi_k > x | M_j, D)
  
  // compute constants
  double a0           = phi0 * mu0;
  double b0           = phi0 * (1 - mu0);
  double logBeta_a0b0 = R::lbeta(a0, b0);
  
  // 
  // Loop through all partitions, storing the unnormalized posterior model probabilities
  // and the average Pr(\pi_k > x | D, M_j)
  // 
  for ( int j = 0; j < partitionMat.n_cols; j++ ) {
    arma::mat datPart = modelMatrix(datMat, partitionMat.col(j));     // obtain data matrix corresponding to current model
    
    // compute posterior parameters
    arma::vec a_jp   = a0 + datPart.col(0);                    // a0 + sum(y     | y in partition)
    arma::vec b_jp   = b0 + datPart.col(1);                   //  b0 + sum(n - y | (n, y) in partition)
    
    // unnormalized log posterior is 
    // sum_p [ log Beta(a_jp, b_jp) - log Beta(a0, b0) ] = sum_p [ log Beta(a_jp, b_jp) ] - P_j * log Beta (a0, b0)
    for ( int i = 0; i < datPart.n_rows; i++ ) {
      postModelProbs(j) += R::lbeta( a_jp(i), b_jp(i) );                // increment by log(Beta(a_jp(i), b_jp(i))) where i indexes row of data partiton
      postSurvProbs(j)  += R::pbeta( x, a_jp(i), b_jp(i), 0.0, 1.0 );   // increment by log survival function of Beta(a_jp(i), b_jp(i)) random variable evaluated at x
    }
    postModelProbs(j) -= datPart.n_rows * logBeta_a0b0;                 // decrement by numModels * log Beta(a0, b0)
  }
  
  // normalize posterior model probabilities
  postModelProbs = postModelProbs - log(sum(exp(postModelProbs)));
  
  // compute posterior responder probability
  double postProb = sum( exp( postModelProbs + postSurvProbs ) );
  
  // Return posterior model probabilities and survior function via BMA
  return
    List::create(
      _["postProb"]        = postProb,
      _["postModelProbs"]  = exp(postModelProbs)
    );
}


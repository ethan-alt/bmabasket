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
// [[Rcpp::export]]
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




// 
// 
// //' Model matrix of data given model
// //' 
// //' Get the (collapsed) model data specified by a particular partition (model)
// //' 
// //' @param datMat     matrix of \code{(y, n-y)} data
// //' @param partition  vector of indices giving how to partition the data
// //' 
// //' @keywords internal
// //' @noRd
// // [[Rcpp::export]]
// arma::mat modelMatrix2(
//     arma::mat const& datMat,
//     arma::vec const& partition
// ) {
//   arma::mat partMat = datMat;
//   for ( int i = 0; i < (int) partMat.n_rows ; i++ ) {
//     // ith row of result is column sum of datMat corresponding partition vector being equal to i (elements in ith set of partition)
//     partMat.row(i) = arma::sum( datMat.rows( find( partition == partition(i) ) ) , 0 );
//   }
//   return partMat;
// }
// 
// 







//' Model matrix of data given model
//' 
//' Get the (collapsed) model data specified by a particular partition (model)
//' 
//' @param datMat     matrix of \code{(y, n-y)} data
//' @param partition  vector of indices giving how to partition the data
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat modelMatrix2(
    arma::mat const& datMat,
    arma::vec const& partition
) {
  // max index = # of elements of list
  int m = partition.max() + 1;                
  
  // initialize matrix for results
  arma::mat partMat = arma::mat( datMat.n_rows, datMat.n_cols, arma::fill::zeros );
  
  for ( int i = 0; i < m; i++ ) {
    arma::uvec   rowIndx  = arma::find(partition == i);
    arma::rowvec partSum  = arma::sum( datMat.rows( rowIndx ) );
    for ( unsigned int j = 0; j < rowIndx.size(); j++ ) {
      partMat.row(rowIndx(j)) = partSum;
    }
  }
  return partMat;
}



//' Log posterior survival probability
//' 
//' Given a vector that gives a partition, computes the posterior probability 
//' that the basket proportion is larger than some specified value
//' 
//' @param pi0        vector whose elements are between 0 and 1 giving threshold for the hypothesis test. If a scalar is provided, assumes same threshold for each basket.
//' @param datMat     matrix of data. The first column gives the number of "successes" and second gives number of "failures" for each basket.
//' @param partition  vector giving a partition for a particular model
//' @param a0         beta prior shape parameter
//' @param b0         beta prior shape parameter
// [[Rcpp::export]]
List logPostSurvProb(
  arma::vec        pi0,
  arma::mat const& datMat,
  arma::vec const& partition,
  double    const& a0,
  double    const& b0,
  double    const& lbeta_a0b0
) {
  
  // Generate model matrix of the partition
  
    // max index + 1 = size of partition
    int m = partition.max() + 1;                
    
    // initialize matrix for results
    arma::mat partMat = arma::mat( m, datMat.n_cols, arma::fill::zeros );
    
    // Placeholder for posterior parameters
    // arma::vec a_jp = a0 * arma::vec(m, arma::fill::ones);
    // arma::vec b_jp = b0 * arma::vec(m, arma::fill::ones);
    // arma::vec a_jp(m);
    // arma::vec b_jp(m);
    // a_jp.fill(a0);
    // b_jp.fill(b0);
    
    
    double sumLogBeta = -m * lbeta_a0b0;   // initialize to be -m * log[Beta(a0, b0)] where m = partition size
    
    
    for ( int i = 0; i < m; i++ ) {
      
      // store index vector giving which rows in datMat pertain to partition
      arma::uvec partIndx = find( partition == i );
      
      // sum columns of datMat pertaining the i^th element of the partition: --> (y, n-y) for y in partition
      arma::rowvec partDat_i = arma::sum( datMat.rows( partIndx ) , 0 );
      
      // compute posterior parameters
      // a_jp(i) += partDat_i(0);
      // b_jp(i) += partDat_i(1);
      double a_jp = a0 + partDat_i(0);
      double b_jp = b0 + partDat_i(1);
      
      // add on log[Beta(a_jp, b_jp)] to sumLogBeta
      sumLogBeta += R::lbeta(a_jp, b_jp);
      
      
      // compute posterior surival function at pi_0k for each basket in the partition
      for ( unsigned int k = 0; k < partIndx.size(); k++ ) {
        pi0(partIndx(k)) = R::pbeta( pi0(partIndx(k)), a_jp, b_jp, 0.0, 1.0 );
      }
    }

    return List::create(
      _["logPostSurvProb"] = pi0,
      _["sumLogBeta"]      = sumLogBeta
    );
}






//' Bayesian model averaging
//'
//' Computes posterior model probabilities and Bayes model averaged survival function \math{P(\pi_k > pi0 | D)}
//'
//' @param pi0            vector whose elements are between 0 and 1 giving cutoffs for effect size for each basket
//' @param datMat         matrix of data \code{(y, n - y)}
//' @param partitionMat   martrix giving how to partition the data for each model
//' @param mu0            scalar giving prior mean for beta prior
//' @param phi0           scalar giving prior dispersion for beta prior
//' @param logModelPriors vector of length P giving the normalized priors for each model
//'
//' @return a vector giving the Bayesian model averaged posterior probability
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::vec bma_cpp(
        arma::vec const& pi0,
        arma::mat const& datMat,
        arma::mat const& partitionMat,
        double    const& mu0,
        double    const& phi0,
        arma::vec const& logModelPriors
) {
  
  // Store BMA probabilities in empty vector
  arma::vec bmaProbs = arma::vec( pi0.size(), arma::fill::zeros );
  double modelProbDenominator = 0;
  
  // compute beta prior shape parameters
  double a0         = mu0 * phi0;
  double b0         = (1 - mu0) * phi0;
  double lbeta_a0b0 = R::lbeta(a0, b0);
  
  // loop through models (partitions),
  for ( unsigned int j = 0; j < partitionMat.n_cols; j++ ) {
    
    // Obtain Pr(pi_k > pi_0k | Mj, D) and log unnormalized weight for model Mj
    List postModel      = logPostSurvProb( pi0, datMat, partitionMat.col(j), a0, b0, lbeta_a0b0 );
    
    // increment bmaProbs by unweighted posterior model probability; add to denominator of posterior model probs
    double logModelProbNumerator = postModel["sumLogBeta"];       // sum( log(Beta(a_jp, b_jp) - log(B(a0, b0) ) ) )
    arma::vec postProb           = postModel["logPostSurvProb"];  // log[ Pr(pi_k > pi_0k | Mj, D) ]
    
    bmaProbs             += arma::exp( logModelProbNumerator + postProb );
    modelProbDenominator += exp( logModelProbNumerator );
  }
  return bmaProbs / modelProbDenominator;
}














// 
// 
// //' Bayesian model averaging
// //'
// //' Computes posterior model probabilities and Bayes model averaged survival function \math{P(\pi_k > pi0 | D)}
// //'
// //' @param pi0            vector whose elements are between 0 and 1 giving cutoffs for effect size for each basket
// //' @param datMat         matrix of data \code{(y, n - y)}
// //' @param partitionMat   martrix giving how to partition the data for each model
// //' @param mu0            scalar giving prior mean for beta prior
// //' @param phi0           scalar giving prior dispersion for beta prior
// //' @param logchoose      scalar giving \code{sum(lchoose(n, y))}
// //' @param P              integer scalar giving number of models
// //' @param logModelPriors vector of length P giving the normalized priors for each model
// //'
// //' @return a list giving \code{[[posteriorSurvivalFunction], [posteriorModelProbs]]}
// //' @keywords internal
// //' @noRd
// // [[Rcpp::export]]
// List bma_cpp2 (
//     double const& pi0,
//     arma::mat const& datMat,
//     arma::mat const& partitionMat,
//     double    const& mu0,
//     double    const& phi0,
//     double    const& logchoose,
//     arma::vec const& logModelPriors
// ) {
//   
//   // Initialize a vector to store posterior model probabilities and
//   // a double to store the bayes average survival function
//   arma::vec postModelProbs = logModelPriors + logchoose;                  // initialize vector of posterior model probs:P(M_j | D) to be the prior probs + sum(log(choose(n,y)))
//   arma::vec postSurvProbs  = arma::vec(pi0.size(), arma::fill::zeros);    // initialize vector to store Pr(pi_k > pi0k | M_j, D), k = 1, ... #Baskets
//   
//   // compute constants
//   double a0           = phi0 * mu0;
//   double b0           = phi0 * (1 - mu0);
//   double logBeta_a0b0 = R::lbeta(a0, b0);
//   
//   //
//   // Loop through all partitions, storing the unnormalized posterior model probabilities
//   // and the average Pr(\pi_k > pi0k | D, M_j)
//   //
//   for ( unsigned int j = 0; j < partitionMat.n_cols; j++ ) {
//     arma::mat datPart = modelMatrix2(datMat, partitionMat.col(j));     // obtain data matrix corresponding to current model
//     
//     // compute posterior parameters
//     arma::vec a_jp   = a0 + datPart.col(0);                    // a0 + sum(y     | y in partition)
//     arma::vec b_jp   = b0 + datPart.col(1);                   //  b0 + sum(n - y | (n, y) in partition)
//     
//     // unnormalized log posterior is
//     // sum_p [ log Beta(a_jp, b_jp) - log Beta(a0, b0) ] = sum_p [ log Beta(a_jp, b_jp) ] - P_j * log Beta (a0, b0)
//     for ( unsigned int i = 0; i < datPart.n_rows; i++ ) {
//       postModelProbs(j) += R::lbeta( a_jp(i), b_jp(i) );                  // increment by log(Beta(a_jp(i), b_jp(i))) where i indexes row of data partiton
//       postSurvProbs(j)  += R::pbeta( pi0, a_jp(i), b_jp(i), 0.0, 1.0 );   // increment by log survival function of Beta(a_jp(i), b_jp(i)) random variable evaluated at pi0
//     }
//     postModelProbs(j) -= datPart.n_rows * logBeta_a0b0;                   // decrement by numModels * log Beta(a0, b0)
//   }
//   
//   // normalize posterior model probabilities
//   postModelProbs = postModelProbs - log(sum(exp(postModelProbs)));
//   
//   // compute posterior responder probability
//   double postProb = sum( exp( postModelProbs + postSurvProbs ) );
//   
//   // Return posterior model probabilities and survior function via BMA
//   return
//     List::create(
//       _["postProb"]        = postProb,
//       _["postModelProbs"]  = exp(postModelProbs)
//     );
//   
// }
// 


















// 
// 
// 
// 
// //' Bayesian model averaging
// //' 
// //' Computes posterior model probabilities and Bayes model averaged survival function \math{P(\pi_k > pi0 | D)}
// //' 
// //' @param pi0            vector whose elements are between 0 and 1 giving cutoffs for effect size for each basket
// //' @param datMat         matrix of data \code{(y, n - y)}
// //' @param partitionMat   martrix giving how to partition the data for each model
// //' @param mu0            scalar giving prior mean for beta prior
// //' @param phi0           scalar giving prior dispersion for beta prior
// //' @param logchoose      scalar giving \code{sum(lchoose(n, y))}
// //' @param P              integer scalar giving number of models
// //' @param logModelPriors vector of length P giving the normalized priors for each model
// //' 
// //' @return a list giving \code{[[posteriorSurvivalFunction], [posteriorModelProbs]]}
// //' @keywords internal
// //' @noRd
// // [[Rcpp::export]]
// List bma_cpp (
//     double const& pi0,
//     arma::mat const& datMat,
//     arma::mat const& partitionMat,
//     double    const& mu0,
//     double    const& phi0,
//     double    const& logchoose,
//     arma::vec const& logModelPriors
// ) {
//   
//   // Initialize a vector to store posterior model probabilities and
//   // a double to store the bayes average survival function
//   arma::vec postModelProbs = logModelPriors + logchoose;                         // initialize vector of posterior model probs:P(M_j | D) to be the prior probs + sum(log(choose(n,y)))
//   arma::vec postSurvProbs  = arma::vec(partitionMat.n_cols, arma::fill::zeros);  // initialize vector to store Pr(pi_k > pi0 | M_j, D)
//   
//   // compute constants
//   double a0           = phi0 * mu0;
//   double b0           = phi0 * (1 - mu0);
//   double logBeta_a0b0 = R::lbeta(a0, b0);
//   
//   // 
//   // Loop through all partitions, storing the unnormalized posterior model probabilities
//   // and the average Pr(\pi_k > pi0 | D, M_j)
//   // 
//   for ( int j = 0; j < (int) partitionMat.n_cols; j++ ) {
//     arma::mat datPart = modelMatrix(datMat, partitionMat.col(j));     // obtain data matrix corresponding to current model
//     
//     // compute posterior parameters
//     arma::vec a_jp   = a0 + datPart.col(0);                    // a0 + sum(y     | y in partition)
//     arma::vec b_jp   = b0 + datPart.col(1);                   //  b0 + sum(n - y | (n, y) in partition)
//     
//     // unnormalized log posterior is 
//     // sum_p [ log Beta(a_jp, b_jp) - log Beta(a0, b0) ] = sum_p [ log Beta(a_jp, b_jp) ] - P_j * log Beta (a0, b0)
//     for ( int i = 0; i < (int) datPart.n_rows; i++ ) {
//       postModelProbs(j) += R::lbeta( a_jp(i), b_jp(i) );                  // increment by log(Beta(a_jp(i), b_jp(i))) where i indexes row of data partiton
//       postSurvProbs(j)  += R::pbeta( pi0, a_jp(i), b_jp(i), 0.0, 1.0 );   // increment by log survival function of Beta(a_jp(i), b_jp(i)) random variable evaluated at pi0
//     }
//     postModelProbs(j) -= datPart.n_rows * logBeta_a0b0;                   // decrement by numModels * log Beta(a0, b0)
//   }
//   
//   // normalize posterior model probabilities
//   postModelProbs = postModelProbs - log(sum(exp(postModelProbs)));
//   
//   // compute posterior responder probability
//   double postProb = sum( exp( postModelProbs + postSurvProbs ) );
//   
//   // Return posterior model probabilities and survior function via BMA
//   return
//     List::create(
//       _["postProb"]        = postProb,
//       _["postModelProbs"]  = exp(postModelProbs)
//     );
// }
// 
// 
// 
// 
// 


// [[Rcpp::depends(RcppArmadillo)]] 
#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;

            
// basket trial data object;
class datObj {
public:
  
  double y;
  int b;
  double et;
  double at;
  double ft;
};

// model space data object;
class modelObj {
public:
  
  arma::imat models;
  arma::ivec nDistinctParms;
};   



// function to order patient data by time to outcome ascertainment;
bool compare_follow_up_time(const datObj &a, const datObj &b)
{
  return a.ft < b.ft;
}

int factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
} 


double numModels (int K0,int P0)
{   
  //std::vector<std::vector<int> > models;
  int nModels = 0;
  double t;
  for (int p=1;p<=P0;p++)
  {
    t=0;
    for (int j=1;j<=p;j++)
    {
      t = t + std::pow(-1,p-j)*std::pow(j,K0)*R::choose(p,j);
    }
    nModels = nModels + (int)(t/factorial(p));
  }
  //arma::mat models(K0,P0,arma::fill::zeros);
  return nModels;
}


void recursive_fill(int i, int kMin, int kMax, int & K0, int & P0, int & r, arma::irowvec & m, arma::imat & models, arma::ivec & nDistinctParms)
{
  for (int k=kMin;k<=kMax;k++)
  {
    m[i] = k;
    if (i+1==K0)
    {
      if (max(m) <= P0)
      {
        models.row(r) = m;
        nDistinctParms[r] = max(m);
        r++;
      }
    }
    else
    {
      int i2    = i + 1;
      int kMin2 = 1;
      int kMax2 = 0;
      for (int s=0;s<=i;s++){ if (m[s]>kMax2) {kMax2=m[s];}}
      kMax2++;
      
      recursive_fill(i2,kMin2,kMax2,K0,P0,r,m,models,nDistinctParms);	
    }
  }
}

modelObj modelMatrix (int K0,int P0)
{   
  
  modelObj d;		
  
  int M0 = numModels(K0,P0);
  
  d.models.set_size(M0,K0);
  d.models.zeros();
  
  d.nDistinctParms.set_size(M0);
  d.nDistinctParms.zeros();
  
  arma::irowvec m(K0,arma::fill::zeros);
  
  int i    = 0;
  int r    = 0;
  int kMin = 1;
  int kMax = 1;
  
  recursive_fill(i,kMin,kMax,K0,P0,r,m,d.models,d.nDistinctParms);	
  
  return d;
  
}    

std::vector<datObj> simData(const int & K0,const int & I0, const arma::ivec & targSSPer, const arma::vec & rRates, const arma::vec & eScales, const arma::vec & aParms)
{
  // calculate maximum possible number of data points (enough so that the trial could enroll the total N from one basket);
  int basket_nMax = 0;
  for (int i=0;i<I0;i++)
  {
    basket_nMax += targSSPer[i]*K0;
  }
  int nMax = basket_nMax*K0;
  
  
  
  // allocate a data object for one subject;
  datObj d;
  
  // allocate data object container for all possible subjects;
  std::vector<datObj> bData(nMax);
  
  // simulate all possible data
  int idx = 0;
  for (int k=0;k<K0;k++) 
  {
    d.b = k;
    
    double cumTime = 0;
    for (int n=0;n<basket_nMax;n++)
    {
      cumTime   = cumTime + R::rexp(eScales[k]);
      d.et      = cumTime;                           // simulate enrollment time;
      d.at      = R::rnorm(aParms[0],aParms[1]);     // simulate response ascertainment time;
      d.ft      = d.et + d.at;                       // calculate total follow-up time;
      d.y       = R::rbinom(1,rRates[d.b]);          // simulate response;
      
      bData[idx] = d;
      idx++;
    }
  }		
  // put the data in chronological order based on outcome ascertainment;
  std::sort(bData.begin(), bData.end(), compare_follow_up_time);
  
  return bData;
}


// [[Rcpp::export]]
arma::mat simDataR(const int & nSims, const int & K0,const int & I0, const arma::ivec & targSSPer, const arma::vec & rRates, const arma::vec & eRates, const arma::vec & aParms)
{
  
  // calculate enrollment distrution scale parameters;
  arma::vec eScales = 1.0 / eRates;	
  
  // calculate maximum possible number of data points (enough so that the trial could enroll the total N from one basket);
  int basket_nMax = 0;
  for (int i=0;i<I0;i++)
  {
    basket_nMax += targSSPer[i]*K0;
  }
  int nMax = basket_nMax*K0;	
  
  
  arma::mat d(nSims*nMax,5,arma::fill::zeros);
  int r = 0;
  for (int s=0;s<nSims;s++)
  {
    // simulate all possible data
    std::vector<datObj> bData = simData(K0,I0,targSSPer,rRates,eScales,aParms);
    
    for (int n=0;n<nMax;n++)
    {
      d(r,0)=s+1;
      d(r,1)=bData[n].b+1;
      d(r,2)=bData[n].et;
      d(r,3)=bData[n].ft;
      d(r,4)=bData[n].y;
      r++;
    }
    
  }
  return d;
  
} 


// [[Rcpp::export]]
Rcpp::List BMA_Design (int nSims,arma::vec eRates, arma::vec rRates, arma::vec aParms, arma::vec ppEffCrit,arma::vec ppFutCrit, int futOnly, arma::vec rRatesNull,arma::vec rRatesAlt, int minSSFut, int minSSEff, arma::imat minSSEnr,arma::imat maxSSEnr, arma::ivec targSSPer, 
                       int I0,double pmp0, double mu0, double phi0)
{
  // 
  // 
  // Rcout << "nSims = " << nSims << "\n";
  // Rcout << "eRates = " << eRates.t() << "\n";
  // Rcout << "rRates = " << rRates.t() << "\n";
  // Rcout << "meanTime = " << aParms(0) << "\n";
  // Rcout << "sdTime = " << aParms(1) << "\n";
  // Rcout << "ppEffCrit = " << ppEffCrit.t() << "\n";
  // Rcout << "ppFutCrit = " << ppFutCrit.t() << "\n";
  // Rcout << "futOnly =" << futOnly << "\n";
  // Rcout << "rRatesNull = " << rRatesNull.t() << "\n";
  // Rcout << "rRatesAlt = " << rRatesAlt.t() << "\n";
  // Rcout << "minSSFut = " << minSSFut << "\n";
  // Rcout << "minSSEff = " << minSSEff << "\n";
  // Rcout << "minSSEnr = \n" << minSSEnr << "\n";
  // Rcout << "maxSSEnr = \n" << maxSSEnr << "\n";
  // Rcout << "targSSPer = \n" << targSSPer.t() << "\n";
  // Rcout << "I0 = " << I0 << "\n";
  // Rcout << "pmp0 = " << pmp0 << "\n";
  // Rcout << "mu0 = " << mu0 << "\n";
  // Rcout << "phi0 = " << "phi0" << "\n";

  
  // initialize random number generator;
  RNGScope scope;    
  
  
  // calculate the number of baskets;
  int K0 = eRates.size();
  
  // calculate number of models and construct model matrix;
  int M0 = numModels(K0,K0);
  
  //Rcout << "# of Models = " << M0 << std::endl;
  
  modelObj mod = modelMatrix(K0,K0);
  
  
  arma::imat models = mod.models;   
  arma::ivec nDistinctParms = mod.nDistinctParms;
  
  // Rcout << "models = \n" << models << "\n";
  
  
  // calculate enrollment distrution scale parameters;
  arma::vec eScales = 1.0 / eRates;
  
  
  // calculate maximum possible number of data points (enough so that the trial could enroll the total N from one basket);
  int basket_nMax = 0;
  for (int i=0;i<I0;i++)
  {
    basket_nMax += targSSPer[i]*K0;
  }
  int nMax = basket_nMax*K0;
  
  double a0           = mu0*phi0;
  double b0           = (1-mu0)*phi0;
  
  
  // calculate log prior model probabilities;
  arma::rowvec log_prior_model_prob(M0,arma::fill::zeros);
  
  
  double pmax = log((double) max(nDistinctParms));
  for (int m=0;m<M0;m++)
  {
    log_prior_model_prob[m] = exp(   pmp0*( log((double) nDistinctParms[m]) - pmax )   );
  }
  log_prior_model_prob = log_prior_model_prob / sum(log_prior_model_prob);
  log_prior_model_prob = log(log_prior_model_prob);
  
  
  // Rcout << "log_prior_model_prob = " << log_prior_model_prob << "\n";
  
  
  
  
  
  // precalculate all poster normalizing constants to avoid using lgamma function repeatedly;
  // precalculate posterior probabilities for treatment efficacy using each baskets threshold;		
  double log_prior_nc = lgamma(a0+b0) - lgamma(a0) - lgamma(b0);
  arma::mat log_post_nc(nMax+1,nMax+1,arma::fill::zeros);
  
  
  std::vector<arma::mat> post_prob_cdf(K0);
  std::vector<arma::mat> post_prob_cdf_mid(K0);
  
  for (int k=0;k<K0;k++)
  {
    arma::mat pp(nMax+1,nMax+1,arma::fill::zeros);	
    arma::mat pp_mid(nMax+1,nMax+1,arma::fill::zeros);				
    for (int y0=0;y0<=nMax;y0++)
    {
      for (int y1=0;y1<=nMax;y1++)
      {
        if ((y0+y1)<=nMax)
        {
          double a1 = a0+y1;
          double b1 = b0+y0;
          
          log_post_nc(y0,y1)   = lgamma(a1+b1) - lgamma(a1) - lgamma(b1);
          pp(y0,y1)     = R::pbeta(rRatesNull[k],a1,b1,0,0);
          pp_mid(y0,y1) = R::pbeta(rRatesNull[k]*0.5+rRatesAlt[k]*0.5,a1,b1,0,0);
        }
      }		
    }
    post_prob_cdf[k]     = pp;
    post_prob_cdf_mid[k] = pp_mid;
  }
  
  
  
  arma::mat all_PostProbs(nSims,K0);
  arma::mat all_PostMeans(nSims,K0);
  arma::mat all_Efficacy(nSims,K0     ,arma::fill::zeros);
  arma::mat all_Futility(nSims,K0     ,arma::fill::zeros);
  arma::mat all_n(nSims,K0,arma::fill::zeros);
  arma::mat all_piHat(nSims,K0,arma::fill::zeros);
  arma::mat all_pInterim(nSims,I0,arma::fill::zeros);
  
  arma::mat all_bias(nSims,K0,arma::fill::zeros);
  arma::mat all_mse(nSims,K0,arma::fill::zeros);
  
  arma::vec all_durations(nSims,arma::fill::zeros);
  arma::vec all_fwr(nSims,arma::fill::zeros);
  arma::vec all_numerr(nSims,arma::fill::zeros);
  arma::mat all_cont(nSims,I0,arma::fill::zeros);
  
  
  
  
  // arma::imat yMat(K0,2,arma::fill::zeros);	
  // arma::ivec nVec(K0,  arma::fill::zeros);				
  // arma::rowvec ppModel(M0,arma::fill::zeros);	
  
  //#pragma omp parallel for num_threads(4)
  for (int s=0;s<nSims;s++)
  {
    // allocate a data object for one subject;
    //datObj d;
    
    // simulate all possible data
    std::vector<datObj> bData = simData(K0,I0,targSSPer,rRates,eScales,aParms);
    
    // put the data in chronological order based on outcome ascertainment;
    //std::sort(bData.begin(), bData.end(), compare_follow_up_time);
    
    
    // container for whether enrollment is open in a basket (1=Yes;0=No)
    arma::irowvec active(K0,arma::fill::ones);
    
    // hypothesis testing result for basket (1=Efficacy;-1=Futility;0=Indeterminate)
    arma::irowvec decision(K0,arma::fill::zeros);
    
    // variable to store at which interim the study stopped;
    int final_interim = 0;
    
    // varuabke ti store the final study duration;
    double duration = 0;
    
    // total sample size
    int nStart = 0;
    
    // number of additional outcomes to observe before next analysis
    int numAdd = 0;
    
    // containers for final study results;
    arma::rowvec basket_specific_pp(K0,  arma::fill::zeros);	
    arma::rowvec basket_specific_pp_mid(K0,  arma::fill::zeros);			
    arma::rowvec basket_specific_mn(K0,  arma::fill::zeros);			
    
    
    // containers for final study data;
    arma::imat yMat(2,K0,arma::fill::zeros);		
    arma::irowvec nVec(K0,  arma::fill::zeros);				
    
    
    // loop over interim analyzes;
    for (int i=0;i<I0;i++)
    {
      
      // determine number of additional outcomes needed (ideally);
      numAdd = sum(active)*targSSPer[i];
      
      // find the number of subjects needed to meet requirements on minimum enrollment
      // without exceeding requirements on maximum enrollment;
      arma::rowvec nAcc(K0,arma::fill::zeros);
      arma::rowvec nMinCritMet(K0,arma::fill::zeros);
      arma::rowvec nMaxCritMet(K0,arma::fill::zeros);
      
      for (int k=0;k<K0;k++)
      {
        if (minSSEnr(i,k)==0) { nMinCritMet[k]=1; }
      }
      
      int idx2=0;
      for (int n=nStart;n<nMax;n++)
      {
        int b = bData[n].b;
        int y = bData[n].y;
        
        if (active[b]==1 and nMaxCritMet[b]==0)
        {
          nAcc[b] +=1;
          if (nAcc[b] >= minSSEnr(i,b)) { nMinCritMet[b]=1; }
          if (nAcc[b] >= maxSSEnr(i,b)) { nMaxCritMet[b]=1; }
          
          yMat(y,b)   += 1;
          nVec[b]     += 1;	
          duration    = bData[n].ft;						
          
          if ( (sum(nAcc)>= numAdd) and (sum(nMinCritMet)==sum(active)) ) {n=nMax + 100;}
        }	
        idx2++;
      }		
      
      // update total number of subjects
      nStart += idx2;
      
      /*if (sum(nVec)>nMax)
       {
       Rcout << "nAcc=" << nAcc <<std::endl;
       Rcout << "nVec=" << nVec <<std::endl;
       Rcout << "nMinCritMet=" << nMinCritMet <<std::endl;
       Rcout << "nMaxCritMet=" << nMaxCritMet <<std::endl;
       }*/
      
      
      // BMA model fitting section;
      
      arma::rowvec ppModel = log_prior_model_prob; // container for posterior model probability container; 	
      arma::mat pp(M0,K0,arma::fill::zeros);       // container for model-specific posterior probabilities for efficacy;
      arma::mat pp_mid(M0,K0,arma::fill::zeros);   // container for model-specific posterior probabilities for efficacy;				
      arma::mat mn(M0,K0,arma::fill::zeros);       // container for model-specific posterior means;
      
      for (int m=0;m<M0;m++)
      {
        arma::irowvec mID = models.row(m)-1;   // basket parameter assignments;
        
        int D0            = max(mID)+1;        // number of distinct parameters;
        
        arma::imat y(D0,2,arma::fill::zeros);  // create container for compressed data;
        
        // put data into compressed buckets;
        for (int k=0;k<K0;k++)
        {
          int d    = mID[k];
          
          y(d,0)     += yMat(0,k);
          y(d,1)     += yMat(1,k);				
        }
        
        
        // compute unnormalized posterior model probabilities;
        for (int d=0;d<D0;d++)
        {
          ppModel[m] += log_prior_nc - log_post_nc(y(d,0),y(d,1));
        }
        
        
        // compute model-specific posterior probabilities;
        for (int k=0;k<K0;k++)
        {
          int d       = mID[k];
          pp(m,k)     = post_prob_cdf[k](y(d,0),y(d,1));
          pp_mid(m,k) = post_prob_cdf_mid[k](y(d,0),y(d,1));						
          mn(m,k) = (a0 + y(d,1)) / (a0 + b0 + y(d,1) +  y(d,0));
        }		
      }
      
      
      double modMax = max(ppModel),sumProb;
      ppModel = exp(ppModel - modMax);
      sumProb = sum(ppModel);
      ppModel = ppModel / sumProb;
      
      // compute model averaged posterior mean and posterior probability of efficacy;
      basket_specific_pp     = ppModel * pp;
      basket_specific_pp_mid = ppModel * pp_mid;
      basket_specific_mn     = ppModel * mn;
      
      int prev_active = sum(active);
      for (int k=0;k<K0;k++)
      {
        int eval_crit_met = (nVec[k]>=minSSEff) or ((prev_active==1) and (active[k]==1)) or (i==(I0-1));
        if ((eval_crit_met==1) and (futOnly==0 or (i==(I0-1))) and (basket_specific_pp[k]>=ppEffCrit[k])) { active[k] = 0; decision[k] =  1; }
        
        eval_crit_met = ((nVec[k]>=minSSFut) or ((prev_active==1) and (active[k]==1))) and (i<(I0-1));
        if ((eval_crit_met==1) and (basket_specific_pp_mid[k]<=ppFutCrit[k]))                             { active[k] = 0; decision[k] = -1; }
      }
      
      if (futOnly==1 and (i<(I0-1)) )
      {
        arma::irowvec active_Fut=active;
        arma::irowvec decision_Fut=decision;
        
        prev_active = sum(active);
        int poss_stop   = 0;
        for (int k=0;k<K0;k++)
        {
          if (active[k]==1 and (nVec[k]>=minSSEff) and (basket_specific_pp[k]>=ppEffCrit[k]))
          {
            poss_stop       += 1;
            active_Fut[k]    = 0;
            decision_Fut[k]  = 1;
          }
        }
        
        if (poss_stop==prev_active)
        {
          active   = active_Fut;
          decision = decision_Fut;
        }
      }
      
      
      all_cont(s,i) = sum(active);
      
      if (max(active)==0 or i==(I0-1)) 
      { 
        final_interim = i+1; 
        i=10000;
      }
      
      
    }
    
    all_PostProbs.row(s)     = basket_specific_pp;
    all_PostMeans.row(s)     = basket_specific_mn;
    
    for (int k=0;k<K0;k++)
    {
      all_bias(s,k) = basket_specific_mn[k] - rRates[k];
      all_mse(s,k)  = all_bias(s,k)*all_bias(s,k);
    }
    
    for (int k=0;k<K0;k++)
    {
      if (decision[k]==1)       {all_Efficacy(s,k) = 1.0;}
      else if (decision[k]==-1) {all_Futility(s,k) = 1.0;}
      
      if ((decision[k]==1) and (rRates[k]<=rRatesNull[k])) { all_fwr[s] = 1; all_numerr[s]+=1;}
    }
    
    arma::rowvec nVecDouble = arma::conv_to<arma::rowvec>::from(nVec);
    arma::mat yMatDouble    = arma::conv_to<arma::mat>::from(yMat);
    arma::rowvec yVecDouble = yMatDouble.row(1);
    
    
    all_n.row(s) = nVecDouble;
    all_piHat.row(s) = yVecDouble/nVecDouble;
    
    all_pInterim(s,final_interim-1) += 1;
    
    all_durations[s] = duration;
    
    
    
  }
  
  
  double NERR = 0;
  if (mean(all_fwr)>0) { NERR = mean(all_numerr)/mean(all_fwr);}
  
  
  
  Rcpp::List HT = Rcpp::List::create(  
    Rcpp::Named("rr")     = mean(all_Efficacy),
    Rcpp::Named("fw.fpr") = mean(all_fwr),
    Rcpp::Named("nerr")   = NERR,
    Rcpp::Named("fut")    = mean(all_Futility)								  
  );
  
  
  Rcpp::List PE = Rcpp::List::create(  
    Rcpp::Named("PM.ave") = mean(all_PostMeans),
    Rcpp::Named("SP.ave") = mean(all_piHat),
    Rcpp::Named("PP.ave") = mean(all_PostProbs),
    Rcpp::Named("bias")   = mean(all_bias),
    Rcpp::Named("mse")    = mean(all_mse)
  );
  
  
  Rcpp::List SS = Rcpp::List::create(  
    Rcpp::Named("basket.ave")  = mean(all_n),
    
    Rcpp::Named("basket.med")  = median(all_n),								  
    Rcpp::Named("basket.min")  = min(all_n),
    Rcpp::Named("basket.max")  = max(all_n),	
    
    Rcpp::Named("overall.ave")  = mean(sum(all_n,1),0)
    /*
     ,
     Rcpp::Named("overall.med")  = median(sum(all_n,1),0),								  
     Rcpp::Named("overall.min")  = min(sum(all_n,1),0),
     Rcpp::Named("overall.max" ) = max(sum(all_n,1),0)
     */
  );														
  
  Rcpp::List DU = Rcpp::List::create(  
    Rcpp::Named("average")      = mean(all_durations)
    /*
     ,
     Rcpp::Named("median")       = median(all_durations),								  
     Rcpp::Named("minimum")      = min(all_durations),
     Rcpp::Named("maximum")      = max(all_durations)
     */
  );	
  
  Rcpp::List SP = Rcpp::List::create(  
    Rcpp::Named("interim.stop.prob")           = mean(all_pInterim),
    Rcpp::Named("baskets.continuing.ave")      = mean(all_cont)			
  );	
  
  
  
  return Rcpp::List::create(
    Rcpp::Named("hypothesis.testing")    = HT,	
    Rcpp::Named("sample.size")           = SS,										
    Rcpp::Named("point.estimation")      = PE,
    Rcpp::Named("trial.duration")        = DU,
    Rcpp::Named("early.stopping")        = SP							  			  
  );
}   


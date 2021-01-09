options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx = as.numeric(args[1]);


if (.Platform$OS.type == "windows") { node.idx = 1 }


root.path   = '~/Projects/bmabasket'
pkg.path    = file.path(root.path, 'bmabasket')
setwd(pkg.path)
devtools::load_all()
source.path = file.path(root.path, 'internalTest')
Rcpp::sourceCpp(file.path(source.path, 'matt_bma.cpp'));

setwd(source.path);

nSims     = 250000;                     
aParms    = c(0.01,0.0000000001);
meanTime  = aParms[1]
sdTime    = aParms[2]

mu0       = 0.45;
phi0      = 1.00;
ppEffCrit = 0.985;
ppFutCrit = 0.2750;
pmp.scale = 2;
n1        = 7
n2        = 16
futOnly   =1;
parms     = expand.grid(mu0,phi0,ppEffCrit,ppFutCrit,pmp.scale,n1,n2,futOnly);

set.seed(node.idx);
start.time = proc.time();

K0  = 5
row = 0;
mss = 4



rTarg     = 0.45;
rNull     = 0.15;
rRatesMod = matrix(rNull,(K0+1)+3,K0);

for (i in 1:K0)  
{
  rRatesMod[(i+1):(K0+1),i]= rTarg     
}
rRatesMod[(K0+2),]=c(0.05,0.15,0.25,0.35,0.45);
rRatesMod[(K0+3),]=c(0.15,0.30,0.30,0.30,0.45);
rRatesMod[(K0+4),]=c(0.15,0.15,0.30,0.30,0.30);

for (idx in (1:nrow(parms)))
{
  mu0  = parms[idx,1];
  phi0 = parms[idx,2];
  
  
  ppEffCrit = rep(parms[idx,3],K0);
  ppFutCrit = rep(parms[idx,4],K0);
  
  pmp0 = parms[idx,5];
  
  targSSPer    = c(parms[idx,6],parms[idx,7]) ;
  
  futOnly      = parms[idx,8];
  
  nInterim     = 2;
  
  
  minSSFut     = mss;  ## minimum number of subjects in basket to assess futility using BMA;
  minSSEff     = mss;  ## minimum number of subjects in basket to assess activity using BMA;
  
  minSSEnr     = matrix(rep(mss,K0),nrow=nInterim,ncol=K0,byrow=T); ## minimum # of new subjects per basket before next analysis - each row is an interim;
  
  maxSSEnr     = matrix(rep(100,K0),nrow=nInterim,ncol=K0,byrow=T); ## maximum # of new subjects per basket before next analysis - each row is an interim;
  
  
  rRatesNull    = rep(rNull,K0);
  rRatesMid     = rep(rTarg,K0);
  
  for (s in c(1.0,0.5,2.0))
  {
    
    if (s==2.0) { eRatesMod = rep(1,K0); } else { eRatesMod = rep(2,K0);  }
    
    .s = 0;
    .e = K0+3;
    if (s!=1.0) {.e = K0;}
    
    
    
    for (i in .s:.e)  
    {
      if (i>=1 & i<=K0) { eRatesMod[i] = eRatesMod[i]*s } else if (i>K0) { eRatesMod = rep(2,K0); }
      # x <- BMA_Design(nSims,eRatesMod,rRatesMod[i+1,],aParms,ppEffCrit,ppFutCrit,futOnly,
      #                 rRatesNull,rRatesMid,minSSFut,minSSEff,minSSEnr,maxSSEnr,targSSPer,nInterim,
      #                 pmp0,mu0,phi0);
      
      
      x <- bma_design(
        nSims, K0, K0, eRatesMod, rRatesMod[i+1,], meanTime, sdTime, ppEffCrit, ppFutCrit, as.logical(futOnly), rRatesNull, rRatesMid, minSSFut,
        minSSEff, minSSEnr, maxSSEnr, targSSPer, nInterim, mu0, phi0, priorModelProbs = NULL, pmp0 = pmp0
      )
      
      row = row +1;
      rr = as.data.frame(rbind(rRatesMod[i+1,]));
      
      parm.sim = cbind(mss,parms[idx ,],s,sum(rRatesMod[i+1,]>rNull),rr );
      colnames(parm.sim) <- c("mss","mu0","phi0","ppEffCrit","ppFutCrit","pmp","n1_targ","n2_targ","futOnly","enrScale","numAlt",paste("trr",seq(1,K0),sep=""));
      
      if (row ==1) { y = cbind(parm.sim,rbind(unlist(x)))
      } else         y = rbind(y,cbind(parm.sim,rbind(unlist(x))))
      
    }
  }
  
}

stop.time = proc.time();
elapsed.time = (stop.time-start.time)/60


names = colnames(y)
names.new = gsub("point.estimation.PM.ave",                  "PM",names)
names.new = gsub("point.estimation.SP.ave",                  "SP",names.new)
names.new = gsub("point.estimation.PP.ave",                  "PP",names.new)
names.new = gsub("point.estimation.bias",                    "BIAS",names.new)
names.new = gsub("point.estimation.mse",                     "MSE",names.new)

names.new = gsub("hypothesis.testing.rr",                    "rr",names.new)
names.new = gsub("hypothesis.testing.fut",                   "fut",names.new)
names.new = gsub("hypothesis.testing.fw.fpr",                "FWER",names.new)
names.new = gsub("hypothesis.testing.nerr",                  "nerr",names.new)


names.new = gsub("sample.size.basket.ave",                   "aveSS",names.new)
names.new = gsub("sample.size.basket.med",                   "medSS",names.new)
names.new = gsub("sample.size.basket.min",                   "minSS",names.new)
names.new = gsub("sample.size.basket.max",                   "maxSS",names.new)
names.new = gsub("sample.size.overall.ave",                  "aveSSovr",names.new)
names.new = gsub("sample.size.overall.med",                  "medSSovr",names.new)
names.new = gsub("sample.size.overall.min",                  "minSSovr",names.new)
names.new = gsub("sample.size.overall.max",                  "maxSSovr",names.new)

names.new = gsub("trial.duration.average",                  "aveDur",names.new)
names.new = gsub("trial.duration.median" ,                  "medDur",names.new)
names.new = gsub("trial.duration.maximum",                  "maxDur",names.new)
names.new = gsub("trial.duration.minimum",                  "minDur",names.new)

names.new = gsub("early.stopping.interim.stop.prob",        "esProb",names.new)
names.new = gsub("early.stopping.baskets.continuing.ave",   "nbCont",names.new)



colnames(y) <- names.new;
rownames(y) <- NULL;

name = paste("K",K0,"-DESIGN-OPTIMAL-",formatC(node.idx, width = 4, format = "d", flag = "0"),"_PKG.CSV",sep="")
write.csv(y,file=name ,row.names=F)


stop.time = proc.time();
elapsed.time = (stop.time-start.time)/60






# nSims = 1
# Rcpp::sourceCpp(file.path(source.path, 'matt_bma.cpp'));
# 
# set.seed(1)
# matt <- BMA_Design(nSims,eRatesMod,rRatesMod[i+1,],aParms,ppEffCrit,ppFutCrit,futOnly,
#                 rRatesNull,rRatesMid,minSSFut,minSSEff,minSSEnr,maxSSEnr,targSSPer,nInterim,
#                 pmp0,mu0,phi0);
# 
# set.seed(1)
# ethan <- bma_design(
#   nSims, K0, K0, eRatesMod, rRatesMod[i+1,], meanTime, sdTime, ppEffCrit, ppFutCrit, as.logical(futOnly), rRatesNull, rRatesMid, minSSFut,
#   minSSEff, minSSEnr, maxSSEnr, targSSPer, nInterim, mu0, phi0, priorModelProbs = NULL, pmp0 = pmp0
# )
# 


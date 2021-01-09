options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
node.idx = as.numeric(args[1]);


if (.Platform$OS.type == "windows") { node.idx = 1 }
setwd('~/Projects/Psioda_GRA')


# mu0=0.45 phi0=1 ppEffCrit=0.985 ppFutCrit=0.275 pmp=2 n1_targ=7 n2_targ=16 futOnly=1

nSims     = 200000;
aParms    = c(0.01, 0.0000000001)
meanTime  = aParms[1]
sdTime    = aParms[2]
mu0       = 0.45;
phi0      = 1.0;
ppEffCrit = 0.985
ppFutCrit = 0.275
pmp.scale = 2;
n1        = 7
n2        = 16
futOnly   = TRUE;
parms     = expand.grid(mu0,phi0,ppEffCrit,ppFutCrit,pmp.scale,n1,n2,futOnly);
parms     = parms[which(parms$Var6+parms$Var7>=18),]
parms     = parms[sample(nrow(parms)),]
nInterim  = 2;
minSSFut  = 5;  ## minimum number of subjects in basket to assess futility using BMA;
minSSEff  = 5;  ## minimum number of subjects in basket to assess activity using BMA;

names(parms) = c('mu0', 'phi0', 'ppEffCrit', 'ppFutCrit', 'pmp.scale', 'n1', 'n2', 'futonly')


head(parms);

nPerNode = 250;
nrow(parms)/nPerNode;

set.start = 1  + nPerNode *(node.idx-1)
set.stop  = min(nPerNode *(node.idx),nrow(parms))


parms = parms[seq(set.start,set.stop),];
nParmSettings = nrow(parms);

set.seed(node.idx);
start.time = proc.time();

K0 = 5
eRates    = rep(2,K0);
row = 0;

rTarg     = 0.45;
rNull     = 0.15;
rRatesMod = matrix(rNull,(K0+1),K0);

for (i in 1:K0)  
{
  rRatesMod[(i+1):(K0+1),i]= rTarg     
}
for (idx in (1:nParmSettings))
{
  mu0          = parms[idx,1];
  phi0         = parms[idx,2];
  ppEffCrit    = rep(parms[idx,3],length(eRates));
  ppFutCrit    = rep(parms[idx,4],length(eRates));
  pmp0         = parms[idx,5];
  targSSPer    = c(parms[idx,6],parms[idx,7]) ;
  futOnly      = parms[idx,8];
  rRatesNull   = rep(rNull,K0);
  rRatesMid    = rep(rTarg,K0);
  
  minSSEnr     = matrix(c( rep(5,K0),rep(5,K0),
                           rep(5,K0),rep(5,K0)),nrow=nInterim,byrow=T);      ## minimum # of new subjects per basket before next analysis - each row is an interim;
  
  maxSSEnr     = matrix(c( rep(100,K0),rep(100,K0),
                           rep(100,K0),rep(100,K0)),nrow=nInterim,byrow=T); ## maximum # of new subjects per basket before next analysis - each row is an interim;
  ## make sure dimensions are correct
  minSSEnr = minSSEnr[, 1:K0]
  maxSSEnr = maxSSEnr[, 1:K0]
  
  for (s in c(1.0,0.5))
  {
    
    eRatesMod     = eRates;
    
    if (s==1.0) {.s=0;.e=K0} else {.s=1;.e=2;}
    
    
    for (i in .s:.e)  
    {
      if (i>=1 & i<=K0) { eRatesMod[i] = eRates[i]*s }
      x <- bma_design(
        nSims, nBaskets = length(eRates), maxDistinct = length(eRates), 
        eRatesMod, rRatesMod[i+1,], meanTime, sdTime, ppEffCrit, ppFutCrit, futOnly,
        rRatesNull, rRatesMid, minSSFut, minSSEff, minSSEnr, maxSSEnr, targSSPer, I0 = nInterim,
        mu0, phi0, pmp0 = pmp0
      )
      
      row = row +1;
      rr = as.data.frame(rbind(rRatesMod[i+1,]));
      
      parm.sim = cbind(parms[idx ,],s,sum(rRatesMod[i+1,]>rNull),rr );
      if (parm.sim[1,5]<0) { parm.sim[1,5] = -parm.sim[1,5]; }
      colnames(parm.sim) <- c("mu0","phi0","ppEffCrit","ppFutCrit","pmp","n1_targ","n2_targ","futOnly","enrScale","numAlt",paste("trr",seq(1,K0),sep=""));
      
      if (row ==1) { y = cbind(parm.sim,rbind(unlist(x)))
      } else         y = rbind(y,cbind(parm.sim,rbind(unlist(x))))
      
    }
  }
  
  stop.time = proc.time();
  elapsed.time = (stop.time-start.time)/60
  print(c(idx,elapsed.time[[3]],elapsed.time[[3]]/idx)); flush.console();
}

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

setwd('~/Projects/bmabasket/internalTest')
name = paste("./K",K0,"-DESIGN-",formatC(node.idx, width = 4, format = "d", flag = "0"),".CSV",sep="")
write.csv(y,file=name ,row.names=F)
  
stop.time = proc.time();
elapsed.time = (stop.time-start.time)/60
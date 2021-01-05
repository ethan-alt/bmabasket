



check_var <- function(
  var, varname, length = NULL, length.name = NULL, 
  lower = NULL, lower.strict = TRUE, lower.name = NULL,
  upper = NULL, upper.strict = TRUE, upper.name = NULL, 
  int = FALSE
  ) {
  if ( any(is.na(var)) | any(is.null(var)) ) {
    stop( paste(varname, 'must not be NA or NULL') )
  }
  if ( is.null(length.name) & !(is.null(length)) ) {
    length.name <- as.character(length)
  }
  if ( is.null(lower.name) & !(is.null(lower))) {
    lower.name <- as.character(lower)
  }
  if ( is.null(upper.name) & !(is.null(upper)) ) {
    upper.name <- as.character(upper)
  }
  
  ## check length if applicable
  if (!is.null(length)) {
    if ( length(var) != length ) {
      stop( paste(varname, 'must be of length', length.name) )
    }
  }
  ## check lower if applicable
  if(!is.null(lower)) {
    check <- ifelse(lower.strict, any(var < lower), any(var <= lower))
    if ( check ) {
      msg   <- ifelse(lower.strict, 'greater than', 'greater than or equal to')
      stop( paste( varname, 'must be', msg, lower.name ) )
    }
  }
  ## check upper if applicable
  if(!is.null(upper)) {
    check <- ifelse(upper.strict, any(var > upper), any(var >= upper))
    if ( check ) {
      msg   <- ifelse(upper.strict, 'less than', 'less than or equal to')
      stop( paste( varname, 'must be', msg, upper.name ) )
    }
  }
  ## check if integer if applicable
  if(int) {
    if(any(var %% 1 != 0)) {
      stop( paste(varname, 'must be an integer') )
    }
  }
  return(TRUE)
}


bma_design_checks <- function(
  nSims, nBaskets, maxDistinct, eRates, rRates, meanTime, sdTime, ppEffCrit, ppFutCrit, futOnly = FALSE, rRatesNull, rRatesAlt, minSSFut,
  minSSEff, minSSEnr, maxSSEnr, targSSPer, I0, mu0, phi0, priorModelProbs, pmp0, nModels
) {
  check_var(nSims, 'nSims', length = 1, lower = 1, lower.strict = TRUE, int = TRUE)
  check_var(maxDistinct, 'maxDistinct', length = 1, lower = 1, lower.strict = TRUE, upper = nBaskets, upper.strict = TRUE, upper.name = 'nBaskets', int = TRUE)
  check_var(eRates, 'eRates', length = nBaskets, length.name = 'nBaskets', lower = 0, lower.strict = FALSE)
  check_var(rRates, 'rRates', length = nBaskets, length.name = 'nBaskets', lower = 0, lower.strict = FALSE)
  check_var(meanTime, 'meanTime', length = 1)
  check_var(sdTime, 'sdTime', length = 1, lower = 0, lower.strict = TRUE)
  
  if ( !is.logical(futOnly) ) {
    stop('futOnly must be TRUE or FALSE')
  }
  if ( !futOnly ) {
    check_var(ppEffCrit, 'ppEffCrit', length = nBaskets, length.name = 'nBaskets', lower = 0, lower.strict = TRUE, upper = 1, upper.strict = TRUE)
  }
  check_var(ppFutCrit, 'ppFutCrit', length = nBaskets, length.name = 'nBaskets', lower = 0, lower.strict = TRUE, upper = 1, upper.strict = TRUE)
  check_var(rRatesNull, 'rRatesNull', length = nBaskets, length.name = 'nBaskets', lower = 0, lower.strict = TRUE, upper = 1, upper.strict = TRUE)
  check_var(rRatesAlt, 'rRatesAlt', length = nBaskets, length.name = 'nBaskets', lower = 0, lower.strict = TRUE, upper = 1, upper.strict = TRUE)
  check_var(minSSFut, 'minSSFut', length = 1, lower = 1, lower.strict = TRUE, int = TRUE)
  check_var(minSSEff, 'minSSEff', length = 1, lower = 1, lower.strict = TRUE, int = TRUE)
  check_var(minSSEnr, 'minSSEnr', lower = 1, lower.strict = TRUE, int = TRUE)
  check_var(I0, 'I0', length = 1, lower = 1, lower.strict = TRUE, int = TRUE)
  if ( ( nrow(minSSEnr) != I0 ) | (ncol(minSSEnr) != nBaskets ) ) {
    stop('minSSEnr must be a matrix with I0 rows and nBaskets columns')
  }
  check_var(maxSSEnr, 'maxSSEnr', lower = 1, lower.strict = TRUE, int = TRUE)
  if ( ( nrow(maxSSEnr) != I0 ) | (ncol(maxSSEnr) != nBaskets ) ) {
    stop('maxSSEnr must be a matrix with I0 rows and nBaskets columns')
  }
  check_var(targSSPer, 'targSSPer', length = nBaskets, length.name = 'nBaskets', lower = 1, lower.strict = TRUE, int = TRUE)
  check_var(mu0, 'mu0', length = nBaskets, length.name = 'nBaskets', lower = 0, lower.strict = FALSE, upper = 1, upper.strict = FALSE)
  check_var(phi0, 'phi0', length = nBaskets, length.name = 'nBaskets', lower = 0, lower.strict = FALSE)
  if (!is.null(pmp0)) {
    check_var(pmp0, 'pmp0', lower = 0, lower.strict = TRUE)
  }
  if (!is.null(priorModelProbs)) {
    check_var(priorModelProbs, 'priorModelProbs', length = nModels, length.name = 'nModels', lower = 0, lower.strict = FALSE, upper = 1, upper.strict = FALSE)
  }
  return(TRUE)
}



#' Simulate a BMA design
#'
#' Simulates a BMA design given hyperparameters
#'
#' @param nSims number of simulation studies to be performed
#' @param nBaskets number of baskets
#' @param maxDistinct integer between 1 and \code{nBaskets} giving number of distinct model probabilities to use. Defaults to \code{nBaskets}. It is recommended to call \code{\link{numModels}} to ensure that computation is tractable.
#' @param eRates scalar or vector of Poisson process rates for each basket
#' @param rRates scalar or vector of true response rates for each basket
#' @param meanTime mean parameter for time to outcome ascertainment
#' @param sdTime standard deviation parameter for time to outcome ascertainment
#' @param ppEffCrit scalar or vector giving basket-specific posterior probability threshold for activity (i.e., efficacy). Must be specified if \code{futOnly==FALSE}
#' @param ppFutCrit scalar or vector giving basket-specific posterior probability threshold for futility
#' @param futOnly \code{logical} giving whether design allows only for futility stopping (\code{TRUE} = futility only, \code{FALSE} = both futility and efficacy)
#' @param rRatesNull scalar or vector of basket-specific null hypothesis values (for efficacy determination)
#' @param rRatesAlt scalar or vector of basket-specific hypothesized alternative values (for futility determination)
#' @param minSSFut minimum number of subjects in basket to assess futility
#' @param minSSEff minimum number of subjects in basket to assess activity
#' @param minSSEnr matrix giving minimum number of new subjects per basket before next analysis (each row is an interim analysis, each column is a basket)
#' @param maxSSEnr matrix giving maximum number of new subjects per basket before next analysis (each row is an interim analysis, each column is a basket)
#' @param targSSPer scalar or vector giving target sample size increment for each basket
#' @param I0 maximum number of analyses
#' @param mu0 prior mean for the response probabilities
#' @param phi0 prior dispersion response probabilities
#' @param priorModelProbs vector giving prior probabilities for models. Default is prior of each model is proportional to the number of unique probabilities raised to the power \code{pmp0}
#' @param pmp0 scalar giving power for \code{priorModelProbs}. If \code{pmp0==0}, a uniform prior is used for model probabilities. Defaults to 1. Ignored if \code{priorModelProbs} is not \code{NULL}
#' 
#' @return a list giving aspects of the simulation
#' @export
bma_design <- function(
  nSims, nBaskets, maxDistinct = nBaskets, eRates, rRates, meanTime, sdTime, ppEffCrit = NULL, ppFutCrit, futOnly = FALSE, rRatesNull, rRatesAlt, minSSFut,
  minSSEff, minSSEnr, maxSSEnr, targSSPer, I0, mu0 = 0.5, phi0 = 1, priorModelProbs = NULL, pmp0 = 1
) {
  
  ## Perform initial checks
  check_var(nBaskets, 'nBaskets', length = 1, lower = 1, lower.strict = TRUE, int = TRUE)
  check_var(maxDistinct, 'nBaskets', length = 1, lower = 1, lower.strict = TRUE, upper = nBaskets, upper.strict = TRUE, upper.name = 'nBaskets')
  
  ## Compute matrix of models
  models  <- as.matrix( setparts( restrictedparts(nBaskets, maxDistinct) ) ) - 1  ## -1 for c++ indexing
  nModels <- ncol(models)                                                         ## number of total models
  
  ## For values that could be scalar or vector, replace to make vector or return error
  ## if specified variable has incompatible length.
  scalar_to_vector <- function(var, var.name, length, length.name = as.character(length)) {
    if ( length(var) == length ) {
      return(var)
    } else if ( length(var) == 1 ) {
      return( rep(var, length) )
    } else {
      stop(paste0(var.name, ' must be length 1 or ', length.name))
    }
  }
  eRates     <- scalar_to_vector(eRates, 'eRates', nBaskets, 'nBaskets')
  rRates     <- scalar_to_vector(rRates, 'rRates', nBaskets, 'nBaskets')
  ppEffCrit  <- scalar_to_vector(ppEffCrit, 'ppEffCrit', nBaskets, 'nBaskets')
  ppFutCrit  <- scalar_to_vector(ppFutCrit, 'ppFutCrit', nBaskets, 'nBaskets')
  rRatesNull <- scalar_to_vector(rRatesNull, 'rRatesNull', nBaskets, 'nBaskets')
  rRatesAlt  <- scalar_to_vector(rRatesAlt, 'rRatesAlt', nBaskets, 'nBaskets')
  targSSPer  <- scalar_to_vector(targSSPer, 'targSSPer', nBaskets, 'nBaskets')
  
  ## Default value for futOnly == TRUE--ignored in code, but must be specified for c++ function to work.
  if ( futOnly ) {
    ppEffCrit <- rep(0.5, nBaskets)
  }
  
  ## Perform checks
  bma_design_checks(
    nSims, nBaskets, maxDistinct, eRates, rRates, meanTime, sdTime, ppEffCrit, 
    ppFutCrit, futOnly = FALSE, rRatesNull, rRatesAlt, minSSFut,
    minSSEff, minSSEnr, targSSPer, I0, mu0, phi0, priorModelProbs, pmp0, nModels
  )
  
  ## construct aparms
  aParms = c(meanTime, sdTime)
  
  ## default prior model prob is # sets in partition ^ pmp0
  if ( is.null(priorModelProbs ) ) {
    priorModelProbs <- do.call('pmax', c(as.data.frame(t(models)), na.rm=TRUE) ) + 1   ## gives number distinct params for each model
    priorModelProbs <- priorModelProbs^pmp0
    priorModelProbs <- priorModelProbs / sum(priorModelProbs)
  }
  
  ## Check if prior probabilities sum to 1; normalize otherwise
  if ( sum(priorModelProbs) != 1 ) {
    warning('Normalizing priorModelProbs because elements did not sum to 1')
    priorModelProbs = priorModelProbs / sum(priorModelProbs)
  }
  priorModelProbs <- log( priorModelProbs )    ## convert to log probabilities
  
  ## Call C++ function
  bma_design_cpp(
    nSims, eRates, rRates, aParms, ppEffCrit, ppFutCrit, 
    futOnly, rRatesNull, rRatesAlt, minSSFut, minSSEff, minSSEnr, maxSSEnr, 
    targSSPer, I0, mu0, phi0, models, priorModelProbs
  )
}
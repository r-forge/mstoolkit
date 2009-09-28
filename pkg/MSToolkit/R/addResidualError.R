addResidualError <- function( 
  response,                       #@ numeric vector of response data
  covariance,                     #@ lower triangle or matrix
  errStruc = "additive",          #@ function describing how to apply residual error
#@ For multivariate normal errors i.e. full variance-covariance on residual error
#@ prefer errStruc = "None" by default (like createNormalParameters)
  seed = .deriveFromMasterSeed( ) #@ Random Seed to use
  ) { 
  ################################################################################
  # Mango Solutions, Chippenham SN14 0SQ 2006
  # addResidualError.R Tue Jun 19 16:17:20 BST 2007 @678 /Internet Time/
  #
  # Author: Romain
  ################################################################################
  # Mike K Smith, Sandwich CT13 9NJ 2009
  # Revised to include multivariate residual error structure
  #
  # Mon 28th Sept 2009 21:48 BST
  ################################################################################
  # DESCRIPTION: add residual error to a response
  # KEYWORDS: component:response
  ################################################################################

  .requiredArgs(response, "The `response` variable is required")
  .requiredArgs(covariance, "The `covariance` is required")
  set.seed(seed)

  # For multivariate residual error structure
   npar <- length(diag(covariance))
   
  covariance <- parseCovMatrix( covariance, npar)                  

  errFun <- if(is.function(errStruc)) errStruc else {
    errStruc <- initialChar( errStruc, "nap", 
      "`errStruc` should be `none`,`additive`, `proportional` or a function")
    switch( errStruc, 
	   "n" = function(x,y) y,        # none - residual error is an array
       "a" = function(x,y) x+y,      # additive
       "p" = function(x,y) exp(x+y)  # proportional
       )
  }
  
  if( npar > 1){
    error <- mvrnorm( n=length(response), mu=rep(0,npar),
	Sigma = covariance)
	}
  
  if( npar==1){
  error <- rnorm( length(response), mean = 0, 
    sd = sqrt(covariance[1,1]) )
  }
  
  if( length(formals(errFun)) <2  ){
    ectdStop("The error function should take at least two arguments")
  }
  
  out <- errFun( response, error )
  if( npar==1 & out %!of% "numeric"){
    ectdStop("The error function should return a numeric vector") 
  }
  if( npar>1 & out %!of% "matrix"){
    ectdStop("The error function should return an array for >1 covariance value") 
  }
  if( npar==1 & out %!l% response ){     
    ectdStop(
      "The error function supplied generates a vector that does not have" %.nt%
      "the same length as the response vector supplied" 
      ) 
  }
  out
  
  
}

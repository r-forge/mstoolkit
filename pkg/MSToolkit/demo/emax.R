
# set up the working path
wp <- file.path( tempdir(), "emaxExample")
dir.create(wp)

generateData( replicateN = 5, subjects = 400, treatDoses = c(0, 5, 25, 50, 100), 
conCovNames = c("wt", "age"), conCovMean = c(83, 55) , conCovVCov = c(14,10)^2 , 
  conCovDigits = 1, conCovCrit = "18 <= age <= 65", 
  genParNames = "E0,ED50,EMAX", genParMean = c(2,50,10), genParVCov = diag( c(.5,30,10) ), 
  genParBtwNames = "E0,ED50,EMAX", genParBtwMean = c(0,0,0), genParBtwVCov = diag(3), 
  respEqn = "E0 + ((DOSE * EMAX)/(DOSE + ED50))",  respVCov = 5, 
  interimSubj = ".3,.7", workingPath = wp  
  )
   
emaxCode <- function(data){
  library(DoseResponse)
  with( data, {
    uniDoses <- sort( unique(DOSE))                                                                    
    eFit <- emaxalt( RESP, DOSE )
    outDf <- data.frame( DOSE = uniDoses, 
      MEAN = eFit$dm[as.character(uniDoses)], 
      SE = eFit$dsd[as.character(uniDoses)] )
    outDf$LOWER <- outDf$MEAN - 2 * outDf$SE
    outDf$UPPER <- outDf$MEAN + 2 * outDf$SE
    outDf$N     <- table(DOSE)[ as.character(uniDoses) ]
    outDf 
  }) 
}                                                                                                                   
             
macrocode <- function(data) {
  # making up a t-test
  mu0   <- data$MEAN[ data$DOSE == 0 & data$INTERIM == 0]
  mu100 <- data$MEAN[ data$DOSE == 100 & data$INTERIM == 0]
  n0    <- data$N[ data$DOSE == 0 & data$INTERIM == 0]
  n100  <- data$N[ data$DOSE == 100 & data$INTERIM == 0]
  sd0   <- data$SE[ data$DOSE == 0 & data$INTERIM == 0]
  sd100 <- data$SE[ data$DOSE == 100 & data$INTERIM == 0]
  
  sddiff <- if( n0 == n100 ){
    sqrt( (sd0^2 + sd100^2)  / (n0 + n100) )
  } else {
    sqrt( (1/n0 + 1/n100) * ( (n0-1)*sd0^2 + (n100-1)*sd100^2  ) / (n0+n100-2)  )
  }
  tstat  <- ( mu100 - mu0 ) / sddiff 
  success <- abs(tstat) > qt( .975, n0+n100-2)
  
  data.frame( SUCCESS = success, TSTAT = tstat )
}
  
interimCode <- function( data ){
  dropdose  <- with( data , DOSE [ sign(UPPER) != sign(LOWER) & DOSE != 0] )
  outList <- list()
  if( length(dropdose) > 0 ) outList$DROP <- dropdose
  outList$STOP <- length(dropdose) == nrow(data)-1
  outList
}
   
analyzeData( 1:5, analysisCode = emaxCode, macroCode = macrocode, 
  interimCode = interimCode, workingPath = wp )


  
.onAttach <- function(libname, pkgname ){
  # Add paths from ECTD.ini to the environment
  sourceOut <- try(source( file.path( .path.package("MSToolkit"), "ECTD.ini" ))$value)
  if (class(sourceOut) != "try-error") assign("externalPaths", sourceOut, env = .ectdEnv)
  copyright <- readLines(system.file("COPYRIGHT", package = "MSToolkit" ) )
  copyright <- gsub( "\\$version", packageDescription("MSToolkit", fields = "Version"), copyright)
  cat("\n"); cat( copyright, sep = "\n" )
  
  if( "RUnit" %in% search() || "RUnit" %in% .packages(all = TRUE)) {
	  unitText <- "# Unit Tests: mstoolkitUnitTests( )"
	  cat(unitText, sep="\n")
  }
  if( .checkGridAvailable() ){
	  gridText <- "# Grid execution available: use 'grid' argument in ?analyzeData"
	  cat(gridText, sep = "\n")
  }
}

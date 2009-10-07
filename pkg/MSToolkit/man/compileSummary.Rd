\name{compileSummary}
\alias{compileSummary}

\title{Compile all Micro or Macro evaluation data} 
\description{
This function reads all of the micro evaluation or macro evaluation data and amalgamates it into a single data frame which
has an additional column indicating which replicate each row belongs to.  This data frame is then written to "m[ia]croSummary.csv".  
}
\usage{
compileSummary(dataType = c("MicroEvaluation", "MacroEvaluation"), replicates = NULL, tryPerl = FALSE, workingPath = getwd())
}

\arguments{
  \item{dataType}{ (Optional) A string that indicates which of the two data types should be read. Must be either "MicroEvaluation" or "MacroEvaluation". Partial matching is used here. See \code{\link{match.arg}}.  Default is "MicroEvaluation" }
  \item{replicates}{ (Optional) Which replicates should be compiled. By default, the function will compile all available replicate data at the time }
  \item{tryPerl}{ (Optional) Should the alternative "Perl" script be tried (faster). The perl script will ony be used if perl in available on the system, most likely if it is a unix system.  Default is TRUE }
  \item{workingPath}{ (Optional) The working directory to be used as the root for the "Micro" and "Macro" evaluation directories.  Also the directory in which the compiled summary files will be created.  The current working directory is used by default }
}
\value{                                                                                      
  CompileSummary does not return a value.
}

\author{ Mango Solutions & Pfizer \email{mstoolkit@googlemail.com} }
\examples{
  \dontrun{
    # a directory MicroEvaluation need to be there and be full 
    # of files matching the pattern micro[0-9]{4}.csv
    compileSummary(dataType = "Micro")
    
    # a directory MacroEvaluation need to be there and be full 
    # of files matching the pattern macro[0-9]{4}.csv
    compileSummary(dataType = "Macro")
    
    # the perl script
    file.show( system.file( "exec", "compile", package = "MSToolkit" ))

  }
}
\keyword{IO}

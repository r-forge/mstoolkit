\name{wrapRNMExport}
\alias{wrapRNMExport}
\title{Wraps the RNMExport functionality}
\description{
Wraps the \emph{RNMExport} library functionality into a single function call, which allows users to generate a number of new control files containing a set of new initial parameter values
}
\usage{
wrapRNMExport(sampList, obj, ctlName, path = getwd(), digits = 4, 
    fixThetas = FALSE, writeFiles = TRUE, subDirs = "nmRuns") 
}
\arguments {
        \item{sampList}{ List of new parameter samples, typically generated using \link{createNmParSamples} } 
        \item{obj}{ The NONMEM object to be used as a basis } 
        \item{ctlName}{ Output file name for new control files } 
        \item{path}{ Path containing the control file } 
        \item{digits}{ Number of digits to which to round parameter values before placing them in the control file as initial parameters } 
        \item{fixThetas}{ Logical: Should the Theta values be fixed to the given value in the new control files? } 
        \item{writeFiles}{ Logical: Should the control files be written out? } 
        \item{subDirs}{ Name of the subdirectory to which control files will be written }
}
\details{
The control file text is extracted from the NONMEM object.  For each set of new parameter samples given by \code{sampList}, the \emph{RNMExport} library functionality is used to replace the initial parameters in the control file text
If \code{writeFiles} is TRUE, each updated control file object is written out to a file contained within the "subDirs" directory
}
\value{
A list of updated control file objects will be invisibly returned
}
\author{ Mango Solutions & Pfizer \email{mstoolkit@googlemail.com} }
\keyword{ nonmem }

\name{createTvDefinition}
\alias{createTvDefinition}

\title{ Create "Typical Value" Definition } 
\description{
	Create "Typical Value" Definition for use with "SIMON" typical value simulation
}
\usage{
createTvDefinition(drugName, drugValue, doses, times = NULL, type = c("Parallel", 
    "Crossover"), sequence, testCall = TRUE)
}

\arguments{
  \item{drugName}{ ... }
  \item{drugValue}{ ...  }
  \item{doses}{ ... }
  \item{times}{ ... }
  \item{type}{ ... }
  \item{sequence}{ ... }
  \item{testCall}{ ... }
}
\details{
	Constructs "tvDefinition" object
}
\value{
	A "tvDefinition" object
}

\author{ MSToolkit Support \email{support@mango-solutions.com} } 

\examples{

createTvDefinition("A", 1, 1:3)

}

\keyword{ Typical Value}


\name{checks-class}
\docType{class}
\alias{checks-class}
\title{Class "checks"}
\description{The formal S4 class for checks. This class contains
  transformed version of the raw response data as well as summaries of
  the checks.}

\section{Objects from the Class}{
Object created by a call to function \code{\link{ConjointChecks}}. 
}
\section{Slots}{
  \describe{
    \item{\code{N}:}{matrix containing the number of respondents at each
    item/ability intersection}

    \item{\code{n}:}{matrix containing the number of correct responses
    at each item/ability intersection}

    \item{\code{Checks}:}{List containing information about each checked
      3-matrix}

    \item{\code{tab}:}{matrix containing information about the detected
    violations at each item/ability intersection}

    \item{\code{means}:}{vector containing weighted and unweighted means
    for the detected violations (where weights are the number of
    individuals at each ability level)}

    \item{\code{check.counts}:}{matrix giving the number of times a
    item/ability cell was sampled}
  }
}

\details{
  Objects of class \code{\linkS4class{checks}} contains all information returned by
  \code{\link{ConjointChecks}}.
}


\author{ Ben Domingue \email{ben.domingue@gmail.com} }

\seealso{
\code{\link{ConjointChecks}}, \code{\link{summary.checks}}, \code{\link{plot.checks}}
}

\keyword{classes}

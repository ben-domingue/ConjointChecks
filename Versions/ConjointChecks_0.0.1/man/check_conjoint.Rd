% File src/library/base/man/load.Rd
\name{ConjointChecks}
\alias{ConjointChecks}
\title{Check Single and Double Cancellation in a sample of 3-matrices}
\description{
Given two matrices, \code{n} and \code{N} (which contain the number of
correct responses and the number of total responses for each cell), a
check of single and double cancellation is performed in \code{n.3mat} matrices.
}
\usage{
ConjointChecks(N,n,n.3mat,par.options=NULL)
}
\arguments{
\item{N}{Matrix containing the total number of responses.}
\item{n}{Matrix containing the number of correct responses.}
\item{n.3mat}{Number of 3-matrices to sample.}
\item{par.options}{List of length two. First element is either "snow" or
  "SMP" specifying either the \code{snow} or \code{doSMP}/\code{foreach}
  parallel options. The first is preferred but is not available on
  Windows machines. The second element is the number of processors.}
}
\examples{
library(ConjointChecks)
n.3mat<-6 #not nearly enough, but relatively quick to run.

############################################
#simulated Rasch example
n.items<-20
n.respondents<-2000
#simulate data
rnorm(n.items)->diff
rnorm(n.respondents)->abil
matrix(abil,n.respondents,n.items,byrow=FALSE)-matrix(diff,n.respondents,n.items,byrow=TRUE)->kern
exp(kern)/(1+exp(kern))->pv
runif(n.items*n.respondents)->test
ifelse(pv>test,1,0)->resp
#now check
PrepareChecks(resp)->tmp
ConjointChecks(tmp$N,tmp$n,n.3mat=n.3mat)->out
SummarizeChecks(out)->tab
}

\keyword{iteration} %need to add keyword


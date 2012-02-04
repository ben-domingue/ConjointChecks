% File src/library/base/man/load.Rd
\name{plot.checks}
\alias{plot.checks}
\title{Plot checks produced by \code{ConjointChecks()}.}
\description{
Takes output from \code{ConjointChecks()} and produces a matrix showing
the percentage of reported violations at each cell.
}
\usage{
plot(chain) 
}
\arguments{
\item{chain}{Object returned by \code{ConjointChecks}.}
}
\examples{
library(ConjointChecks)
n.3mat<-5 #not nearly enough, but relatively quick to run.

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

plot(out)

}
\keyword{file}


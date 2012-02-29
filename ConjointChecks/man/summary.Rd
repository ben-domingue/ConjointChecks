% File src/library/base/man/load.Rd
\name{summary.checks}
\alias{summary.checks}
\title{Summarize checks produced by \link{\code{ConjointChecks}}.}
\description{
Takes output from \link{\code{ConjointChecks}} and produces a matrix showing
the percentage of reported violations at each cell.
}
\usage{
\method{summary}{checks}(object, ...)
}
\arguments{
\item{object}{Object returned by \link{\code{ConjointChecks}} of class
  \link{\code{checks}}.}

  \item{\dots}{further arguments passed to or from other methods}

}
\examples{
#simulated Rasch example
n.items<-20
n.respondents<-2000
#simulate data
rnorm(n.items)->diff
rnorm(n.respondents)->abil
matrix(abil,n.respondents,n.items,byrow=FALSE)->m1
matrix(diff,n.respondents,n.items,byrow=TRUE)->m2
m1-m2 -> kern
exp(kern)/(1+exp(kern))->pv
runif(n.items*n.respondents)->test
ifelse(pv>test,1,0)->resp
#now check
PrepareChecks(resp)->tmp
n.3mat<-1 #not nearly enough, but relatively quick to run.
ConjointChecks(tmp$N,tmp$n,n.3mat=n.3mat)->out

summary(out)

}


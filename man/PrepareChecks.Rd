\name{PrepareChecks}
\alias{PrepareChecks}
\title{Prepare raw response data for \code{\link{ConjointChecks}}.}
\description{
Takes output from \code{\link{ConjointChecks}} and produces a matrix showing
the percentage of reported violations at each cell.
}
\usage{
PrepareChecks(resp,ss.lower=10) 
}
\arguments{
\item{resp}{Raw dichotomously coded response data. Columns represent
  items and rows represent individuals.}
\item{ss.lower}{Only sum scores that have at least this many distinct
  individuals with that sum score will be used.}
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
PrepareChecks(resp)->obj
}




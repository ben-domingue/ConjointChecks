% File src/library/base/man/load.Rd
\name{plot.checks}
\alias{plot.checks}
\title{Plot checks produced by \code{ConjointChecks()}.}
\description{
Takes output from \code{ConjointChecks()} and produces a matrix showing
the percentage of reported violations at each cell.
}
\usage{
\method{plot}{checks}(x,items,item.labels,....)
}
\arguments{
  \item{x}{Object returned by \code{ConjointChecks} of class
    \code{checks}.}
  \item{items}{Vector of item numbers to include in a single plot. All
    items can be included, but this is less helpful for diagnostic
    purposes.}
  \item{item.labels}{Should item numbers be included? Defaults to \code{TRUE}. If
    length of \code{\link{items}} is unity (perhaps if the small multiple
    format of Tufte, 2001 is going to be used), then the item number
    gets printed below the x-axis. If the length of items is more than
    unity, the item number gets plotted in the figure above the largest
    proportion of violations for each item.}
  \item{\dots}{further arguments passed to or from other methods}
}
\references{
  Tufte, E. R. (2001). The visual display of quantitative information (2nd ed.). Chesire, CT: Graphics
Press.
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
m1-m2 ->kern
exp(kern)/(1+exp(kern))->pv
runif(n.items*n.respondents)->test
ifelse(pv>test,1,0)->resp
#now check
PrepareChecks(resp)->tmp
n.3mat<-1 #not nearly enough, but relatively quick to run.
ConjointChecks(tmp$N,tmp$n,n.3mat=n.3mat)->out

#this is going to be quite lackluster without a larger number of checks
plot(out)

}


\name{plot.checks}
\alias{plot.checks}
\title{Plot checks produced by \code{\link{ConjointChecks}}.}
\description{
Takes output from \code{\link{ConjointChecks}} and produces a \code{matplot} showing
the percentage of reported violations at each cell.
}
\usage{
\method{plot}{checks}(x,items=NULL,item.labels=TRUE,...)
}
\arguments{
  \item{x}{Object returned by \code{\link{ConjointChecks}} of class
    \code{\linkS4class{checks}}.}
  \item{items}{Vector of item numbers to include in a single
    plot. Defaults to all, but this is less helpful for diagnostic
    purposes.}
  \item{item.labels}{Should item numbers be included? Defaults to \code{TRUE}. If
    length of \code{items} is unity (perhaps if the small multiple
    format of Tufte, 2001 is going to be used), then the item number
    gets printed below the x-axis. If the length of items is more than
    unity, the item number gets plotted in the figure above the largest
    proportion of violations for each item.}
  \item{\dots}{further arguments passed to or from other methods}
}
\value{No return value, called for side effects}
\references{
  Tufte, E. R. (2001). The visual display of quantitative information (2nd ed.). Chesire, CT: Graphics
Press.
}
\examples{
opar <- par()
par(mfrow=c(3,2))
plot(rasch1000)
plot(rasch1000,items=c(5,10,15))
for (i in c(3,9,13,18)) plot(rasch1000,items=i)
par(opar)
}


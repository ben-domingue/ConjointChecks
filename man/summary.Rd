\name{summary.checks}
\alias{summary.checks}
\title{Summarize checks produced by \code{\link{ConjointChecks}}.}
\description{
Takes output from \code{\link{ConjointChecks}} and produces a matrix showing
the percentage of reported violations at each cell.
}
\usage{
\method{summary}{checks}(object, ...)
}
\arguments{
\item{object}{Object returned by \code{\link{ConjointChecks}} of class
  \code{\linkS4class{checks}}.}
  \item{\dots}{further arguments passed to or from other methods}

}
\value{No return value, called for side effects}
\examples{
summary(rasch1000)
}


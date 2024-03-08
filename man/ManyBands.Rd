\name{ManyBands}
\alias{ManyBands}
\title{Check What Ordering to Use}
\description{
Internal funcion. This checks to see whether the p-value (rasch difficulty) ordering should be used or if ordering should be 'as is'.
}
\usage{
ManyBands(th, se, cc.type, resp, bands=seq(10,50,by=10), uniform.bands=TRUE,
          n.workers=NULL, trim.window=NULL, pv.order=TRUE)
}
\arguments{
\item{th}{Threshold.}
\item{se}{Standard error.}
\item{cc.type}{Type of cancellation check.}
\item{resp}{resp.}
\item{bands}{Values of the bands.}
\item{uniform.bands}{Is the distribution of the bands uniform?}
\item{n.workers}{Is the distribution of the bands uniform?}
\item{trim.window}{trim.window.}
\item{pv.order}{Use the p-value ordering.}
}
\references{
Perline, R., Wright, B. D., & Wainer, H. (1979). The Rasch model as
additive conjoint measurement. Applied Psychological Measurement, 3(2), 237-255.
}
\author{Ben Domingue \email{ben.domingue@gmail.com}} 
\examples{
### INTERNAL FUNCTION ###
}


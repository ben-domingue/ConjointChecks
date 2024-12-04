\name{SingleCancel}
\alias{SingleCancel}
\title{Check Single Cancellation in a sample}
\description{
Internal function. Wrapper of the code that checks ONLY the single cancellation bits.
}
\usage{
SingleCancel(N,n,CR=c(.025,.975),single,mc.cores=1)
}
\arguments{
\item{N}{Matrix containing the total number of responses.}
\item{n}{Matrix containing the number of correct responses.}
\item{CR}{Width of the credible region taken from the
posterior. Defaults to a 95\% credible region (\code{c(.025,.975)}).}
\item{single}{Also test single cancellation.}
\item{mc.cores}{The number of cores to parallelize over.}
}
\references{
Perline, R., Wright, B. D., & Wainer, H. (1979). The Rasch model as
additive conjoint measurement. Applied Psychological Measurement, 3(2), 237-255.
}
\author{Ben Domingue \email{ben.domingue@gmail.com}} 
\examples{
### INTERNAL FUNCTION ###
}


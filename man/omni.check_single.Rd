\name{omni.check_single}
\alias{omni.check_single}
\title{Check Single Cancellation in a sample}
\description{
Internal function. Code that checks ONLY the single cancellation bits.
}
\usage{
omni.check_single(N,n,n.iter,burn=1000,thin=4,CR,single)
}
\arguments{
\item{N}{Matrix containing the total number of responses.}
\item{n}{Matrix containing the number of correct responses.}
\item{n.iter}{Total number of samples.}
\item{burn}{Number of initial samples that should be discarded.}
\item{thin}{Amount of thinning.}
\item{CR}{Width of the credible region taken from the
posterior. Defaults to a 95\% credible region (\code{c(.025,.975)}).}
\item{single}{Also test single cancellation.}
}
\references{
Perline, R., Wright, B. D., & Wainer, H. (1979). The Rasch model as
additive conjoint measurement. Applied Psychological Measurement, 3(2), 237-255.
}
\author{Ben Domingue \email{ben.domingue@gmail.com}} 
\examples{
### INTERNAL FUNCTION ###
}


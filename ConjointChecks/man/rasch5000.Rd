\name{rasch5000}
\alias{rasch5000}
\docType{data}
\title{5000 sampled 3-matrices from simulated Rasch data.}
\description{
  Object created by first generating Rasch data and then running \code{ConjointChecks}  on 5000 sampled 3 matrices
}
\usage{rasch5000}
\format{
  An object of class \code{\linkS4class{checks}}.
}
\source{
  Simulated via Rasch model.
}
\keyword{datasets}

% #simulated via
% library(ConjointChecks)
% n.items<-20
% n.respondents<-2000
% #simulate data
% rnorm(n.items)->diff
% rnorm(n.respondents)->abil
% matrix(abil,n.respondents,n.items,byrow=FALSE)->m1
% matrix(diff,n.respondents,n.items,byrow=TRUE)->m2
% m1-m2 -> kern
% exp(kern)/(1+exp(kern))->pv
% runif(n.items*n.respondents)->test
% ifelse(pv>test,1,0)->resp
% ##now check
% PrepareChecks(resp)->tmp
% detectCores()->n.workers
% ConjointChecks(tmp$N,tmp$n,n.3mat=5000,par.options=list(n.workers=n.workers,type="SOCK"),seed=1234)->rasch5000

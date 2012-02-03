library(ConjointChecks)
n.3mat<-100 #not nearly enough, but relatively quick to run.

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
prepare(resp)->tmp
#
n.3mat<-100
par.seed<-as.integer(rep(1234,6))
conjointChecks(tmp$N,tmp$n,n.3mat=n.3mat,par.options=list(c("snow","MPI"),2),seq.seed=3,par.seed=par.seed)->out1

ConjointChecks(tmp$N,tmp$n,n.3mat=n.3mat,par.options=list(c("snow","SOCK"),2),seq.seed=3,par.seed=par.seed)->out2

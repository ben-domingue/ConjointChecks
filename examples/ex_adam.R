library(ConjointChecks)
n.3mat<-20 #not nearly enough, but relatively quick to run.

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
par.seed<-as.integer(rep(1234,6))
conjointChecks(tmp$N,tmp$n,n.3mat=n.3mat,par.options=list(c("snow","MPI"),2),seq.seed=3,par.seed=par.seed)->out1
conjointChecks(tmp$N,tmp$n,n.3mat=n.3mat,par.options=list(c("snow","SOCK"),2),seq.seed=3,par.seed=par.seed)->out2
conjointChecks(tmp$N,tmp$n,n.3mat=n.3mat,par.options=list(c("foreach","doSMP"),2),seq.seed=3,par.seed=par.seed)->out3
conjointChecks(tmp$N,tmp$n,n.3mat=n.3mat,par.options=list(c("foreach","doMC"),2),seq.seed=3,par.seed=par.seed)->out4

#if you wanted
identical(out1,out2)
identical(out2,out3) #maybe false
identical(out3,out4)

#the following is kind of interesting if you turn n.3mat up to something reasonable (1k, 5k, or something like that)
#plot(out1)
#summary(out1)

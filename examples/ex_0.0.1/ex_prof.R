library(ConjointChecks)
n.3mat<-10 #not nearly enough, but relatively quick to run.

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
PrepareChecks(resp)->tmp

Rprof("/tmp/prof_v1.out",memory.profiling=TRUE,interval=.005)
#ConjointChecks(tmp$N,tmp$n,n.3mat=n.3mat,par.options=list(c("snow","MPI"),2))->out1
ConjointChecks(tmp$N,tmp$n,n.3mat=n.3mat,par.options=list(c("foreach","doMC"),2))->out1
Rprof(NULL)
Rprof("/tmp/prof_v2.out",memory.profiling=TRUE,interval=.005)
#ConjointChecks_v2(tmp$N,tmp$n,n.3mat=n.3mat,par.options=list(c("snow","MPI"),2))->out2
ConjointChecks_v2(tmp$N,tmp$n,n.3mat=n.3mat,par.options=list(c("foreach","doMC"),2))->out2
Rprof(NULL)

#processing
library(profr)
par(mfrow=c(2,1))
plot(parse_rprof("/tmp/prof_v1.out"))
plot(parse_rprof("/tmp/prof_v2.out"))

#memory
fun<-function(fn) {
  summaryRprof(fn,memory="tseries")->pro
  pro[,1:4]
}
fun("/tmp/prof_v1.out")->m1
fun("/tmp/prof_v2.out")->m2

par(mfrow=c(2,2))
for (i in 1:4) {
  plot(log(m1[,i]),type="l",xlab=names(m1)[i])
  lines(log(m2[,i]),type="l",col="red")
}
  

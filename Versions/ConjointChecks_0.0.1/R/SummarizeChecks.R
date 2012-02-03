#i should turn this into a generic function
compare<-function(dat,lim) {
  ifelse(dat<lim[[2]] & dat>lim[[1]],0,1)->tab
  tab
}

SummarizeChecks<-function(chain,weight=chain$N) {
  chain$N->N
  chain$n->n
  chain$out->out
  (n/N)->dat
  mat.num<-matrix(0,nrow(N),ncol(N))
  mat.den<-matrix(-1,nrow(N),ncol(N))
  for (i in 1:length(out)) {
    out[[i]]->x
    x[[1]]->ro
    x[[2]]->co
    compare(dat[ro,co],x[[3]])->comp
    for (i in 1:3) for (j in 1:3) {
      1+mat.den[ro[i],co[j]]->mat.den[ro[i],co[j]]
      comp[i,j]+mat.num[ ro[i],co[j] ]->mat.num[ro[i],co[j]]
    }
  }
  ifelse(mat.den<0,NA,mat.den)->mat.den
  mat.den+1->mat.den
  mat.num/mat.den->tab
  mean(tab,na.rm=TRUE)->m1
  mean(tab*weight/sum(weight),na.rm=TRUE)->m2
  list(n=chain$n,N=chain$N,tab=tab,mean=m1,wt.mean=m2,tab.checked=mat.den)
}

PrepareChecks<-function(resp,ss.lower=10,collapse.columns=FALSE) {
  if (any(is.na(resp))) stop("Checks will only work with complete data. Suggestion: remove respondents with missing responses.")
  if (ss.lower==1) {
    message("ss.lower must be greater than 1, setting to 2.")
    ss.lower<-2
  }
  ncol(resp)->n.items
  #reorder columns
  colSums(resp)->cs
  resp[,order(cs)]->resp
  #group by sum scores
  rowSums(resp)->rs
  n<-N<-list()
  table(rs)->tab
  as.numeric(names(tab))[tab>=ss.lower]->lev
  for (s in lev) {
    resp[rs==s,]->tmp
    rep(nrow(tmp),n.items)->N[[as.character(s)]]
    colSums(tmp)->n[[as.character(s)]]
  }
  do.call("rbind",N)->N
  do.call("rbind",n)->n
  if (collapse.columns) {
    colSums(n)->cs
    sort(unique(cs))->cs.index
    n2<-N2<-list()
    for (i in 1:length(cs.index)) {
      cs[i]->lev
      n[,cs==lev]->tmp
      if (is(tmp,"matrix")) rowSums(tmp)->tmp
      tmp->n2[[i]]
      N[,cs==lev]->tmp
      if (is(tmp,"matrix")) rowSums(tmp)->tmp
      tmp->N2[[i]]
    }      
    do.call("cbind",n2)->n
    do.call("cbind",N2)->N
  }
  list(N=N,n=n)
}

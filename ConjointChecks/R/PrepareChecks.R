PrepareChecks<-function(resp,ss.lower=10) {
  if (ss.lower==1) {
    message("ss.lower must be greater than 1, setting to 2.")
    ss.lower<-2
  }
  ncol(resp)->n.items
  #reorder things
  colSums(resp)->cs
  resp[,order(cs)]->resp
  #group by sum scores
  rowSums(resp)->rs
  lev[lev>=ss.lower]->lev
  n<-N<-list()
  sort(unique(rs))->lev
  for (s in lev) {
    resp[rs==s,]->tmp
    rep(nrow(tmp),n.items)->N[[as.character(s)]]
    colSums(tmp)->n[[as.character(s)]]
  }
  do.call("rbind",N)->N
  do.call("rbind",n)->n
  list(N=N,n=n)
}

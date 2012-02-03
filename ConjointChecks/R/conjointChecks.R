##main function
omni.check<-function(N,n,n.iter,burn=1000,thin=4,CR) {#this checks both single and double cancellation
  n/N->dat
  chain<-list()
  #initialize
  inits<-dat
  inits->foo1->foo2
  for (i in 1:nrow(foo1)) foo1[i,]<-i
  for (i in 1:ncol(foo2)) foo2[,i]<-i
  foo1+foo2->index
  sort(as.numeric(dat))->hold
  counter<-1
  for (i in unique(as.numeric(index))) {
    grep(paste("^",i,"$",sep=""),index)->index2
    stop.here<-length(index2)
    counter:(counter+stop.here-1)->replace
    inits[index2]<-hold[replace]
    counter<-max(replace)+1
  }
  inits->old#chain[[1]]
  like<-function(theta,N,n) {
    #choose(N,n)->x1 #this part cancels!
    1->x1
    n*log(theta)->x2
    (N-n)*log(1-theta)->x3
    sum(x2+x3)
  }
  old.ll<-inits
  for (i in 1:nrow(old.ll)) for (j in 1:ncol(old.ll)) like(inits[i,j],N[i,j],n[i,j])->old.ll[i,j]
  dc.counter<-rep(0,n.iter-1)
  #iterate
  for (I in 2:n.iter) {
    #chain[[I-1]]->old
    for (i in 1:nrow(dat)) for (j in 1:ncol(dat)) {
      #get left hand
      if (j==1) lh1<-0 else lh1<-old[i,j-1]
      if (i==1) lh2<-0 else lh2<-old[i-1,j]
      lh3<-0
      if (i==1 & j==3) lh3<-old[3,1]
      #get right hand
      if (j==ncol(dat)) rh1<-1 else rh1<-old[i,j+1]
      if (i==nrow(dat)) rh2<-1 else rh2<-old[i+1,j]
      rh3<-1
      if (i==3 & j==1) rh3<-old[1,3]
      #now make sure double cancellation needs to hold
      if (!(old[2,1]<old[1,2] & old[3,2]<old[2,3])) {
        lh3<-0
        rh3<-1
      } else dc.counter[i]<-1
      #now work everything out all nice like...
      lh<-max(lh1,lh2,lh3)
      if (rh3>lh) rh<-min(rh1,rh2,rh3) else rh<-min(rh1,rh2)
      if (rh<lh) rh<-1
      #sample new point
      runif(1,lh,rh)->draw
      #acceptance ratio
      ar<-2
      like(draw,N[i,j],n[i,j])->new.ll
      if (!(old[i,j] %in% 0:1)) exp(new.ll-old.ll[i,j])->ar
      if (ar>runif(1)) {
        draw->old[i,j]
        new.ll->old.ll[i,j]
      }
    }
    if (I>burn & I%%4==0) old->chain[[as.character(I)]]
    #if (I%%100==0) print(I)
  }
  hi<-lo<-M<-chain[[1]]
  for (i in 1:3) for (j in 1:3) {
    post<-sapply(chain,function(x) x[i,j])
    quantile(post,CR[1])->lo[i,j]
    quantile(post,CR[2])->hi[i,j]
    mean(post)->M[i,j]    
  }
  list(low=lo,high=hi,mean=M)
}


############################################################
#glorified wrapper
conjointChecks<-function(N,n,n.3mat=10,par.options=NULL,seq.seed=NULL,par.seed=NULL,CR=c(.025,.975)) {#N is the number of total tries per cell, n is the number correct
  if (!is.null(seq.seed)) set.seed(seq.seed)
  ## chain.2.ci<-function(l,burn=1000,thin=4) {
  ##   l[-(1:burn)]->l
  ##   length(l)->n
  ##   seq(1,n,thin)->index
  ##   l[index]->l
  ##   hi<-lo<-M<-l[[1]]
  ##   for (i in 1:3) for (j in 1:3) {
  ##     post<-sapply(l,function(x) x[i,j])
  ##     quantile(post,.025)->lo[i,j]
  ##     quantile(post,.975)->hi[i,j]
  ##     mean(post)->M[i,j]    
  ##   }
  ##   list(low=lo,high=hi,mean=M)
  ## }  
  #processing function
  proc.fun<-function(dummy,arg.list) {
    arg.list[[1]]->N
    arg.list[[2]]->n
    arg.list[[3]]->lof
    arg.list[[4]]->CR
    lof[[1]]->omni.check
    #lof[[2]]->chain.2.ci
    #lof[[3]]->compare
    test<-1
    nrow(N)->nrows
    ncol(N)->ncols
    while (test>0) {
      sort(sample(1:nrows,3,replace=FALSE))->rows
      sort(sample(1:ncols,3,replace=FALSE))->cols
      N[rows,cols]->nt
      n[rows,cols]->nc
      nc/nt->dat
      sum (dat==1|dat==0)->test
    }
    omni.check(nt,nc,n.iter=3000,CR=CR)->out
    #chain.2.ci(out)->out
    list(rows,cols,out)->save.dat
    save.dat
  }
  out<-list()
  n/N->dat
  ifelse(abs(dat-.5)<=.5,TRUE,FALSE)->test
  if (!all(test)) stop("There is a problem with n/N, values not between 0 and 1 (inclusive)")
  #list(omni.check,chain.2.ci,compare)->lof
  list(omni.check)->lof
  if (is.null(par.options)) {#sequential first
    for (i in 1:n.3mat) {
      list(N,n,lof,CR)->arg.list
      proc.fun(i,arg.list)->out[[i]]
    }
    #lapply(out,chain.2.ci,burn=1000,thin=4)->posterior
  } else {      #snow
    par.options[[1]][1]->par.type
    par.options[[1]][2]->bonus.option
    par.options[[2]]->n.proc
    list(N,n,lof,CR)->arg.list
    dummy<-list()
    for (i in 1:n.3mat) dummy[[i]]<-i
    if (!is.null(par.seed)) as.integer(par.seed)->par.seed
    switch(par.type,
           "snow"={
             if (!require(snow)) stop("Package 'snow' not available.")
             if (!require(rlecuyer)) stop("Package 'rlecuyer' not available.")
             cl<-makeCluster(n.proc,type=bonus.option,verbose=FALSE)
             if (is.null(par.seed)) par.seed<-round(runif(6)*1000)
             clusterSetupRNG(cl,type="RNGstream",seed=par.seed)
             clusterApply(cl,dummy,proc.fun,arg.list=arg.list)->out #could probably wrap these two things together...
             #clusterApply(cl,out,chain.2.ci,burn=1000,thin=4)->posterior
             stopCluster(cl)
           },"foreach"={ 
             if (!require(foreach)) stop("Package 'foreach' not available.")
             if (!require(doRNG)) stop("Package 'doRNG' not available.")
             switch(bonus.option,
                    "doSMP"={
                      if (!require(doSMP)) stop("Package 'doSMP' not available.")
                      w <- startWorkers(n.proc)
                      registerDoSMP(w)
                      "stopWorkers(w)"->final.cmd
                    },"doMC"={
                      if (!require(doMC)) stop("Package 'doMC' not available.")
                      registerDoMC(n.proc)
                      final.cmd<-NULL
                    }
                    )
             if (is.null(par.seed)) par.seed <- doRNGseed()
             doRNGseed(par.seed)
             (foreach(i = unlist(dummy)) %dorng% proc.fun(dummy=i,arg.list=arg.list))->out
             #(foreach(i = iter(out)) %dorng% chain.2.ci(i,burn=1000,thin=4))->posterior
             if (!is.null(final.cmd)) eval(parse(text=final.cmd))
           }
           )
  }
  #for (i in 1:length(out)) posterior[[i]]->out[[i]][[3]]
  list(N=N,n=n,Checks=out)->out
  #now do some summarizing
  compare<-function(dat,lim) {
    ifelse(dat<lim[[2]] & dat>lim[[1]],0,1)->tab
    tab
  }
  out$N->N
  out$n->n
  out$Checks->out
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
  weight<-N
  sum(tab*weight,na.rm=TRUE)/sum(weight)->m2
  list(n=n,N=N,tab=tab,mean=m1,wt.mean=m2,tab.checked=mat.den)
  new("checks", N=N,n=n,Checks=out,tab=tab,means=list(unweighted=m1,weighted=m2),check.counts=mat.den)
  #
}



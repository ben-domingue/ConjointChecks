plot.checks<-function(x, items, item.labels=TRUE, ...) {
  matplot(x@tab[,items],xlab="",xaxt="n",ylab="Proportion Violations",type="l",lty=1:length(items),col="black",...)
  mtext(side=1,line=1,paste("Increasing Ability"))
  if (item.labels) {
    if (length(items)==1) mtext(side=1,line=2,paste("Item",items)) else {
      apply(x@tab,2,which.max)->maxes.x
      apply(x@tab,2,max,na.rm=TRUE)->maxes.y
      for (i in 1:length(items)) text(maxes.x[items],maxes.y[items],items,pos=3)
    }
  }
}


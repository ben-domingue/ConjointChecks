summary.checks<-function(checks) {
  list(Means=checks@means,items=colMeans(checks@tab,na.rm=TRUE))
}

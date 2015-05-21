ManyBands<-function(th,se,cc.type,resp,bands=seq(10,50,by=10),n.workers=NULL) {
    banding.fun<-function(banding,theta,theta.se) { #banding is a vector of cutpoints (no -Inf or Inf)
        cut(theta,c(-Inf,banding,Inf))->cl
        fun<-function(x,se,lims) {
            as.character(lims)->lims
            substr(lims,2,nchar(lims))->lims
            substr(lims,1,(nchar(lims)-1))->lims
            strsplit(lims,",")[[1]]->lims
            as.numeric(lims[1])->lo
            as.numeric(lims[2])->hi
            #
            -abs((lo-x))/se->l1
            abs(x-hi)/se->l2
            pnorm(l1)->p1
            pnorm(l2)->p2
            p2-p1
        }
        Vectorize(fun)->fun
        fun(theta,theta.se,cl)->pv
        -sum(log(pv))
    }
    cc.fun<-function(th,se,banding,cc.type,resp) {
        #making matrices for ConjointChecks
        cut(th,c(-Inf,banding,Inf),ordered_result=TRUE)->cl
        N<-n<-list()
        colSums(resp)->cs
        resp[,order(cs)]->resp
        for (lev in levels(cl)) {
            cl==lev -> index
            resp[index,,drop=FALSE]->tmp
            rep(nrow(tmp),ncol(tmp))->N[[as.character(lev)]]
            colSums(tmp)->n[[as.character(lev)]]
        }
        do.call("rbind",N)->N
        do.call("rbind",n)->n
        which(N[,1]==0) -> index
        if (length(index)>0) {
            N[-index,]->N
            n[-index,]->n
        }
        if (is.null(n.workers)) detectCores()->n.workers
        ConjointChecks(N,n,
                       n.3mat=cc.type,
                       par.options=list(n.workers=n.workers,type="PSOCK")
                       )->out
        summary(out)$Means$weighted->viw
        summary(out)$Means$unweighted->viu
        c(viw,viu)
    }
    hold<-list()
    #for (len in c(10,25,50,75,100)) for (offset in c(-.005,0,.005)) {
    for (len in bands) {
        S<-0
        qu.low<-0.01
        while(S==0) {
            quantile(th,qu.low)->qu1
            sum(th<qu1) -> S
            qu.low<-qu.low+.005
        }
        S<-0
        qu.high<-0.99
        while(S==0) {
            quantile(th,qu.high)->qu2
            sum(th>qu2) -> S
            qu.high<-qu.high-.005
        }
        #print(c(qu1,qu2))
        #quantile(th,.01)->qu1
        #quantile(th,.99)->qu2
        seq(qu1,qu2,length.out=len)->banding
        banding.fun(banding,th,se)->vp
        cc.fun(th,se,banding,cc.type,resp)->cc.out
        list(len=len,banding=banding,vp=vp,cc.out=cc.out)->zz
        c(len,rev(zz$cc.out),zz$vp)->hold[[as.character(len)]]
    }
    do.call("rbind",hold)->tab
    colnames(tab)<-c("n.bands","vp.unweight","vp.weight","stringency")
    data.frame(tab)
}


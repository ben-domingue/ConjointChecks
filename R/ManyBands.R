ManyBands<-function(th,se,cc.type,resp,bands=seq(10,50,by=10),
                    uniform.bands=TRUE,
                    trim.window=NULL,
                    pv.order=TRUE, #this checks to see whether the p-value (rasch difficulty) ordering should be used or if ordering should be 'as is'
                    mc.cores=1
                    ) {
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
    cc.fun<-function(th,se,banding,cc.type,resp,trim.window) { #trim window should be given in SD of theta units
        #order columns. this happens first, prior to trimming.
        if (pv.order) {
            colSums(resp)->cs
            resp[,order(cs)]->resp
        }
        #trimming
        #plot(density(th)); for (jjj in 1:length(banding)) abline(v=banding[jjj],lty=2)
        nrow(resp)->n.ppl.orig
        if (!is.null(trim.window)) {
            matrix(th,length(th),length(banding),byrow=FALSE)->m.th
            matrix(banding,length(th),length(banding),byrow=TRUE)->m.ba
            abs(m.th-m.ba) -> del
            apply(del,1,min) -> min.dist
            trim.gap<-sd(th,na.rm=TRUE)*trim.window
            min.dist>trim.gap -> trim.test
            th[trim.test] -> th
            se[trim.test] -> se
            resp[trim.test,] -> resp
        }
        nrow(resp)->n.ppl.new
        #banding
        cut(th,c(-Inf,banding,Inf),ordered_result=TRUE)->cl
        #making matrices for ConjointChecks
        N<-n<-list()
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
        ConjointChecks(N,n,
                       n.3mat=cc.type,
                       mc.cores=mc.cores
                       )->out
        summary(out)$Means$weighted->viw
        summary(out)$Means$unweighted->viu
        list(n.ppl.new/n.ppl.orig,c(viw,viu))
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
        if (uniform.bands) {
            seq(qu1,qu2,length.out=len)->banding
        } else {
            quantile(th,seq(qu.low-0.005,qu.high+0.005,length.out=len))->banding
        }
        banding.fun(banding,th,se)->vp
        cc.fun(th,se,banding,cc.type,resp,trim.window=trim.window)->cc.out
        list(len=len,banding=banding,vp=vp,trim.ratio=cc.out[[1]],cc.out=cc.out[[2]])->zz
        c(len,zz$trim.ratio,rev(zz$cc.out),zz$vp)->hold[[as.character(len)]]
    }
    do.call("rbind",hold)->tab
    colnames(tab)<-c("n.bands","trim.ratio","vp.unweight","vp.weight","stringency")
    data.frame(tab)
}


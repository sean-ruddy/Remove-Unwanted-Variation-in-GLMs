RUVg4.poi <- function (y, x1, x2=NULL, cIdx, k, linkFun=c("logv","log1"))
 {
    linkFun=paste0(match.arg(linkFun,c("logv","log1")),"(L)")
    coFun <- function( x , y, off=NULL)
             {
                logv <- function(L)
                     {
                       ## link
                       linkfun <- function(y) log2((y+0.5)*1e6/(L+1))

                       ## inverse link
                       linkinv <- function(eta)  2^eta*(L+1)/1e6 - 0.5

                       ## derivative of invlink wrt eta
                       mu.eta <- function(eta) log(2)*2^eta*(L+1)/1e6

                       valideta <- function(eta) TRUE
                       link <- "logv"
                       structure(list(linkfun = linkfun, linkinv = linkinv,
                                      mu.eta = mu.eta, valideta = valideta,
                                      name = link),
                                 class = "link-glm")

                     }
                log1 <- function(L)
                     {
                       ## link
                       linkfun <- function(y) log((y+1)/(L+1))

                       ## inverse link
                       linkinv <- function(eta)  exp(eta)*(L+1) - 1

                       ## derivative of invlink wrt eta
                       mu.eta <- function(eta) exp(eta)*(L+1)

                       valideta <- function(eta) TRUE
                       link <- "log1"
                       structure(list(linkfun = linkfun, linkinv = linkinv,
                                      mu.eta = mu.eta, valideta = valideta,
                                      name = link),
                                 class = "link-glm")

                     }
              if( is.null(off) )
               {
                   return( sapply(1:ncol(y),FUN=function(i)
                            {
                              fit = tryCatch(glm(y[,i]~x-1,family=poisson(link=eval(parse(text=linkFun)))),error=function(e) switch(linkFun,"logv(L)"=lm(log2((y[,i]+0.5)/(L+1)*1e+06)~x-1),"log1(L)"=lm(log((y[,i]+1)/(L+1))~x-1)),silent=TRUE)
                              coefs = coef(fit)
                              return( coefs )
                            } ))
               }
              else return( sapply(1:ncol(y),FUN=function(i)
                           {
                             fit = tryCatch(glm(y[,i]~x-1,offset = off[,i],family=poisson(link=eval(parse(text=linkFun)))),error=function(e) switch(linkFun,"logv(L)"=lm(log2((y[,i]+0.5)/(L+1)*1e+06)~x-1,offset = off[,i]),"log1(L)"=lm(log((y[,i]+1)/(L+1))~x-1,offset = off[,i])),silent=TRUE)
                             coefs = coef(fit)
                             return( coefs )
                           }))
             }

    if( is.null(x2) ) x2 = matrix(1,nc=1,nr=ncol(y))
    if( is.null(ncol(x1)) ) x1 = matrix(x1,nc=1,nr=length(x1))
    if( is.null(ncol(x2)) ) x2 = matrix(x2,nc=1,nr=length(x2))

    dge <- DGEList(counts=y)
    dge <- calcNormFactors(dge)
    lib.size = dge$samples[,2]*dge$samples[,3] #TMM normalization
    L = lib.size
    Y <- switch( linkFun , "logv(L)" = log2(t(y+0.5)/(L+1)*1e+06) , "log1(L)" = log((t(y)+1)/(L+1)) ) #this is the transformation that voom uses
    Yc <- Y[, cIdx]
    ty <- t(y)
    tyc <- ty[,cIdx]
    m <- nrow(Y)
    n <- ncol(Y)
    nc <- ncol(Yc)

    D = cbind(x1 , x2)
    co = coFun(D , ty)
    Z = Y - D %*% co

    svdWa <- svd(Z)
    W_o <- svdWa$u[, (1:k), drop = FALSE]

    D = cbind(D, W_o)
    co = coFun(D , ty)
    colnames(co) <- colnames(Y)
    ac = co[(ncol(D)-ncol(W_o)+1):(ncol(D)-ncol(W_o)+k),cIdx,drop=FALSE]
    x2c = co[(ncol(x1)+1):(ncol(x1)+ncol(x2)),cIdx,drop=FALSE]

    OFF = W_o %*% ac + x2 %*% x2c
    Z = matrix(as.vector(tyc),nc=1,nr=m*nc)
    OFF = matrix(as.vector(OFF),nc=1,nr=m*nc)
    x1.expand = apply( x1 , 1 , FUN=function(vec) rep(vec,each=k))
    if( is.null(dim(x1.expand)) ) x1.expand = matrix(x1.expand,nr=k)
    D = do.call("rbind", lapply( 1:ncol(ac), FUN=function(i) t(x1.expand * ac[,i]) ) )
    L = rep(lib.size,nc)
    co = coFun(D, Z, OFF)
    if( is.null(dim(co)) ) co=matrix(co,nr=k)

    W = W_o + do.call("+",lapply(1:ncol(x1),FUN=function(i) x1[,i,drop=FALSE] %*% t(co[((i-1)*k+1):(i*k),])))

    return( W )
 }


RUVg4.voom <- function (y, x1, x2=NULL, cIdx, k)
 {
    coFun <- function( x , y, off=NULL, weights = NULL , retWeights=FALSE )
             {
                if( is.null(off) )
                 {
                    dge <- DGEList(counts=y)
                    dge <- calcNormFactors(dge)
                    v=voom(dge,x)
                    weights=t(v$weights)
                    fit=lmFit(v,x)
                    coefs=t(fit$coefficients)
                    if( retWeights )
                      return( rbind( coefs , weights ) )
                    else return( coefs )
                 }
               else
                {
                  fit = lm(y~x-1,weights=weights,offset=off)
                  return( coef(fit) )
                }
             }

    if( is.null(x2) ) x2 = matrix(1,nc=1,nr=ncol(y))
    if( is.null(ncol(x1)) ) x1 = matrix(x1,nc=1,nr=length(x1))
    if( is.null(ncol(x2)) ) x2 = matrix(x2,nc=1,nr=length(x2))

    dge <- DGEList(counts=y)
    dge <- calcNormFactors(dge)
    lib.size = dge$samples[,2]*dge$samples[,3] #TMM normalization
    Y <- log2(t(y+0.5)/(lib.size+1)*1e+06) #this is the transformation that voom uses
    Yc <- Y[, cIdx]
    m <- nrow(Y)
    n <- ncol(Y)
    nc <- ncol(Yc)

    D = cbind(x1 , x2)
    co = coFun(D , y, retWeights = TRUE)
    weights = co[(ncol(D)+1):nrow(co),]
    co = co[1:ncol(D),]
    Z = Y - D %*% co

    tmp = bwpca(Z,matw=weights,npcs=k,center=FALSE)
    W_o = apply(tmp$scores,2,FUN=function(vec) vec/sqrt(sum(vec^2)))

    D = cbind(D, W_o)
    co = coFun(D , y, retWeights=TRUE)
    ac = co[(ncol(D)-ncol(W_o)+1):(ncol(D)-ncol(W_o)+k),cIdx,drop=FALSE]
    x2c = co[(ncol(x1)+1):(ncol(x1)+ncol(x2)),cIdx,drop=FALSE]
    weights = co[(ncol(D)+1):nrow(co),]

    OFF = W_o %*% ac + x2 %*% x2c
    Z = matrix(as.vector(Yc),nc=1,nr=m*nc)
    Weights = matrix(as.vector(weights[,cIdx]),nc=1,nr=m*nc)
    OFF = matrix(as.vector(OFF),nc=1,nr=m*nc)
    x1.expand = apply( x1 , 1 , FUN=function(vec) rep(vec,each=k))
    if( is.null(dim(x1.expand)) ) x1.expand = matrix(x1.expand,nr=k)
    D = do.call("rbind", lapply( 1:ncol(ac), FUN=function(i) t(x1.expand * ac[,i]) ) )
    co = coFun( D , Z , OFF , Weights ) #give weights??? which weights??? Giving the original weight vector.
    #Can also give no weights or the weights from the second GLM fit
    if( is.null(dim(co)) ) co=matrix(co,nr=k)

    W = W_o + do.call("+",lapply(1:ncol(x1),FUN=function(i) x1[,i,drop=FALSE] %*% co[((i-1)*k+1):(i*k)]))

    return( W )
 }

RUV4.ols <- function (y, x1, x2=NULL, cIdx, k, linkFun=c("logv","log1"))
 {
    linkFun=paste0(match.arg(linkFun,c("logv","log1")),"(L)")
    coFun <- function( x , y, off=NULL)
             {
                if( is.null(off) )
                 {
                    fit=lmFit(y,x)
                    coefs=t(fit$coefficients)
                    return( coefs )
                 }
               else
                {
                  fit = lm(y~x-1,offset=off)
                  return( coef(fit) )
                }
             }

    if( is.null(x2) ) x2 = matrix(1,nc=1,nr=ncol(y))
    if( is.null(ncol(x1)) ) x1 = matrix(x1,nc=1,nr=length(x1))
    if( is.null(ncol(x2)) ) x2 = matrix(x2,nc=1,nr=length(x2))

    dge <- DGEList(counts=y)
    dge <- calcNormFactors(dge)
    lib.size = dge$samples[,2]*dge$samples[,3] #TMM normalization
    Y <- switch(linkFun, "logv(L)" = log2(t(y+0.5)/(lib.size+1)*1e+06), "log1(L)" = log(t(y+1)/(lib.size+1))) #this is the transformation that voom uses

    Yc <- Y[, cIdx]
    m <- nrow(Y)
    n <- ncol(Y)
    nc <- ncol(Yc)

    D = cbind(x1 , x2)
    co = coFun(D , t(Y))
    Z = Y - D %*% co

    svdWa <- svd(Z %*% t(Z))
    W_o <- svdWa$u[, (1:k), drop = FALSE]

    D = cbind(D, W_o)
    co = coFun(D , t(Y))
    ac = co[(ncol(D)-ncol(W_o)+1):(ncol(D)-ncol(W_o)+k),cIdx,drop=FALSE]
    x2c = co[(ncol(x1)+1):(ncol(x1)+ncol(x2)),cIdx,drop=FALSE]
    OFF = W_o %*% ac + x2 %*% x2c
    Z = matrix(as.vector(Yc),nc=1,nr=m*nc)
    OFF = matrix(as.vector(OFF),nc=1,nr=m*nc)
    x1.expand = apply( x1 , 1 , FUN=function(vec) rep(vec,each=k))
    if( is.null(dim(x1.expand)) ) x1.expand = matrix(x1.expand,nr=k)
    D = do.call("rbind", lapply( 1:ncol(ac), FUN=function(i) t(x1.expand * ac[,i]) ) )
    co = coFun( D , Z , OFF ) #give weights??? which weights??? Giving the original weight vector.
    #Can also give no weights or the weights from the second GLM fit
    if( is.null(dim(co)) ) co=matrix(co,nr=k)

    W = W_o + do.call("+",lapply(1:ncol(x1),FUN=function(i) x1[,i,drop=FALSE] %*% co[((i-1)*k+1):(i*k)]))

    return( W )
 }

RUVg2.voom <- function(y, x1, x2=NULL, cIdx, k, ret.counts=FALSE , subW = NULL )
{
    coFun <- function( x , y, off=NULL, weights = NULL , retWeights=FALSE )
             {
                if( is.null(off) )
                 {
                    dge <- DGEList(counts=y)
                    dge <- calcNormFactors(dge)
                    v=voom(dge,x)
                    weights=t(v$weights)
                    fit=lmFit(v,x)
                    coefs=t(fit$coefficients)
                    if( retWeights )
                      return( rbind( coefs , weights ) )
                    else return( coefs )
                 }
               else
                {
                  fit = lm(y~x-1,weights=weights,offset=OFF)
                  return( coef(fit) )
                }
             }

    if( is.null(x2) ) x2 = matrix(1,nc=1,nr=ncol(y))
    if( is.null(ncol(x1)) ) x1 = matrix(x1,nc=1,nr=length(x1))
    if( is.null(ncol(x2)) ) x2 = matrix(x2,nc=1,nr=length(x2))

    dge <- DGEList(counts=y)
    dge <- calcNormFactors(dge)
    lib.size = dge$samples[,2]*dge$samples[,3] #TMM normalization
    Y <- log2(t(y+0.5)/(lib.size+1)*1e+06) #this is the transformation that voom uses
    Yc <- Y[, cIdx]
    m <- nrow(Y)
    n <- ncol(Y)
    nc <- ncol(Yc)


    D=x2
    co = coFun(D,y,retWeights=TRUE)
    weights = co[(ncol(D)+1):nrow(co),]
    co = co[1:ncol(D),]
    Z = Y - D %*% co
    tmp = bwpca(Z[,cIdx],matw=weights[,cIdx],npcs=k,center=FALSE)
    W = apply(tmp$scores,2,FUN=function(vec) vec/sqrt(sum(vec^2)))

    if( ret.counts )
     {
       if( !is.null(subW) ) W = subW
       D=cbind(x1,x2,W)
       co=coFun(D, y)
       a <- co[(ncol(x1)+ncol(x2)+1):(ncol(x1)+ncol(x2)+k),,drop=FALSE]
       Ynorm <- Y - W %*% a
       ynorm <- pmax(round(2^Ynorm*(lib.size+1)/1e+06-0.5,0),0)
       return(t(ynorm))
     }
    else return( W )
}



RUV2.ols <- function(y, x1, x2=NULL, cIdx, k, linkFun=c("logv","log1"))
{

    linkFun=paste0(match.arg(linkFun,c("logv","log1")),"(L)")
    coFun <- function( x , y, off=NULL)
             {
                if( is.null(off) )
                 {
                    fit=lmFit(y,x)
                    coefs=t(fit$coefficients)
                    return( coefs )
                 }
               else
                {
                  fit = lm(y~x-1,offset=off)
                  return( coef(fit) )
                }
             }

    if( is.null(x2) ) x2 = matrix(1,nc=1,nr=ncol(y))
    if( is.null(ncol(x1)) ) x1 = matrix(x1,nc=1,nr=length(x1))
    if( is.null(ncol(x2)) ) x2 = matrix(x2,nc=1,nr=length(x2))

    dge <- DGEList(counts=y)
    dge <- calcNormFactors(dge)
    lib.size = dge$samples[,2]*dge$samples[,3] #TMM normalization
    Y <- switch(linkFun, "logv(L)" = log2(t(y+0.5)/(lib.size+1)*1e+06), "log1(L)" = log(t(y+1)/(lib.size+1))) #this is the transformation that voom uses
    Yc <- Y[, cIdx]
    m <- nrow(Y)
    n <- ncol(Y)
    nc <- ncol(Yc)


    D=x2
    co = coFun(D,t(Y))
    Z = Y - D %*% co

    svdWa <- svd(Z[,cIdx])
    W <- svdWa$u[, (1:k), drop = FALSE]
    return(W)

}

##Needs to be tested with larger x2 (at least 3 variables)
##And tested with larger x1 (at least 4 groups)
##For RUV4, do we make W orthogonal to x2 as well? Probably not. Otherwise, the X\Gamma term doesn't make sense. It should then
##be (X,Z)\Gamma

getW <- function( counts , cIdx, k , x1=NULL, x2=NULL,
                  method = c("RUV2","RUVg2v","RUV4","RUVg4poi","RUVg4v","RUVg4perp","RUVg4moderated") ,
                  linkFun=c("logv","log1"), ret.counts=FALSE , k2=5)
         {
           if( k==0 ) return( NULL )
           #For getting W:
             #RUV2: voom transform. no weights. std linear model to intercept. Same as RUVg2 except for the tranform
             #RUVg2v: voom transform. voom weights. weighted linear model to intercept. Same as RUV2 above but uses weights in the PCA step.
             #RUV4: voom transform. no weights. std linear model.
             #RUVg4poi: log-e+1 transform. no weights. poisson GLM.
             #RUVg4v: voom transform. voom weights. weighted linear model.
           method = match.arg(method,c("RUV2","RUVg2v","RUV4","RUVg4poi","RUVg4v","RUVg4perp","RUVg4moderated"))
           cIdx2 <- rownames(counts)%in%cIdx
            W = switch(method, RUV2 = RUV2.ols(y=counts, x1=x1, x2=x2, cIdx=cIdx, k=k, linkFun=linkFun),
                               RUVg2v = RUVg2.voom(y=counts, x1=x1, x2=x2, cIdx=cIdx, k=k, ret.counts=ret.counts),
                               RUV4 = RUV4.ols(y=counts, x1=x1, cIdx=cIdx, k=k, linkFun=linkFun),
                               RUVg4poi = RUVg4.poi(y=counts, x1=x1, cIdx=cIdx, k=k, linkFun=linkFun),
                               RUVg4v = RUVg4.voom(y=counts, x1=x1, cIdx=cIdx, k=k),
                               RUVg4perp = RUVg4.perp(y=counts, x1=x1, cIdx=cIdx, k=k, linkFun=linkFun),
                               RUVg4moderated = RUVg4.moderated(y=counts, x1=x1, cIdx=cIdx, k=k, k2=k2, linkFun=linkFun))
           return(W)
         }


getDE <- function(counts,x1,x2=NULL,W,groups=NULL,fdr=0.05,coef=1,pIdx=NULL, deMethod=c("edgeR","voom","empvar","lm","limma","deseq2","wald"),use.weights=FALSE)
          {

            deMethod = match.arg(deMethod,c("edgeR","voom","empvar","lm","limma","deseq2","wald"))
            eRFun <- function(counts, x1, x2, W, groups, fdr, coef)
                      {
                          D=cbind(x1,x2,W)
                          y <- DGEList(counts = counts, group = groups)
                          y <- calcNormFactors(y, method = "upperquartile")
                          y <- estimateGLMCommonDisp(y, D)
                          y <- estimateGLMTagwiseDisp(y, D)
                          fit <- glmFit(y, D)
                          lrt <- glmLRT(fit, coef=coef)
                          tbl <- topTags(lrt, n=Inf)$table
                          pvals <- tbl$PValue
                          names(pvals) <- rownames(tbl)
                          de <- rownames(tbl[tbl$FDR<=fdr,])
                          posDE=NULL
                          posAll=NULL
                          if( !is.null(pIdx) )
                           {
                            posDE <- match(pIdx,de)
                            posDE[is.na(posDE)] <- Inf
                            posDE <- sort(posDE)
                            posAll <- match(pIdx,names(pvals))
                            posAll[is.na(posAll)] <- Inf
                            posAll <- sort(posAll)
                            posAll[which(posDE==Inf)] <- -1*posAll[which(posDE==Inf)]
                           }
                          pvals <- pvals[rownames(counts)]
                          return( list( "Pvals" = pvals , "PosDE" = posDE , "PosRank" = posAll ) )
                      }
            vFun <- function(counts, x1, x2, W, fdr, coef)
                     {
                        D=cbind(x1,x2,W)
                        dge <- DGEList(counts=counts)
                        dge <- calcNormFactors(dge)
                        v=voom(dge,D)
                        fit=lmFit(v,D)
                        fit2 = eBayes(fit)
                        tbl <- topTable(fit2,coef=coef,number=nrow(counts),sort.by="P")
                        pvals <- tbl$P.Value
                        names(pvals) <- rownames(tbl)
                        de <- rownames(tbl[tbl$adj.P.Val<=fdr,])
                        posDE=NULL
                        posAll=NULL
                        if( !is.null(pIdx) )
                         {
                          posDE <- match(pIdx,de)
                          posDE[is.na(posDE)] <- Inf
                          posDE <- sort(posDE)
                          posAll <- match(pIdx,names(pvals))
                          posAll[is.na(posAll)] <- Inf
                          posAll <- sort(posAll)
                          posAll[which(posDE==Inf)] <- -1*posAll[which(posDE==Inf)]
                         }
                        pvals <- pvals[rownames(counts)]
                        return( list( "Pvals" = pvals , "PosDE" = posDE , "PosRank" = posAll ) )
                     }
         evFun <- function(counts, x1, x2, W, groups, fdr, coef, use.weights = FALSE)
                        {
                            n = nrow(counts)
                            p = ncol(x1)
                            varbetahat = pvals = tvals = matrix(0, p, n)

                            dge <- DGEList(counts=counts)
                            dge <- calcNormFactors(dge)
                            lib.size = dge$samples[,2]*dge$samples[,3] #TMM normalization
                            Y <- log2(t(counts+0.5)/(lib.size+1)*1e+06) #this is the transformation that voom uses
                            D = cbind(x1, x2, W)
                            weights=NULL
                            if( use.weights ) weights <- voom(dge,D)$weights

                            fit=lmFit(t(Y),D,weights=weights) #same as lm(transformed_counts~D, weights=weights)
                            betahat=t(fit$coefficients)[1:p,,drop=FALSE]
                            sigma2 = fit$sigma^2 #verified that this is the true sigma ie sqrt(sum(w*r^2)/df), where df = nsmp - ncol(D)

                            for (l in 1:p) {
                                varbetahat[l, ] = get_empirical_variances(sigma2,betahat[l, ], bin = 10, rescaleconst = NULL)
                                tvals[l, ] = betahat[l, ]/sqrt(varbetahat[l,])
                                 pvals[l, ] = 2 * pt(-abs(tvals[l, ]), Inf)
                            }
                            pvals <- pvals[coef,]
                            names(pvals) <- rownames(counts)
                            padj <- p.adjust(pvals,method="BH")
                            de <- names(padj[padj<=fdr])
                            posDE=NULL
                            posAll=NULL
                            if( !is.null(pIdx) )
                             {
                              posDE <- match(pIdx,de)
                              posDE[is.na(posDE)] <- Inf
                              posDE <- sort(posDE)
                              posAll <- match(pIdx,names(pvals))
                              posAll[is.na(posAll)] <- Inf
                              posAll <- sort(posAll)
                              posAll[which(posDE==Inf)] <- -1*posAll[which(posDE==Inf)]
                             }
                            pvals <- pvals[rownames(counts)]
                            return( list( "Pvals" = pvals , "PosDE" = posDE , "PosRank" = posAll ) )
                         }
            lmFun <- function(counts, x1, x2, W, fdr, coef)
                     {
                        n = nrow(counts)
                        p = ncol(x1)
                        dge <- DGEList(counts=counts)
                        dge <- calcNormFactors(dge)
                        lib.size = dge$samples[,2]*dge$samples[,3] #TMM normalization
                        Y <- log2(t(counts+0.5)/(lib.size+1)*1e+06) #this is the transformation that voom uses
                        D = cbind(x1, x2, W)
                        fit=lmFit(t(Y),D,weights=NULL) #same as lm(transformed_counts~D, weights=weights)
                        betahat=t(fit$coefficients)[coef,]
                        sigma2 = fit$sigma^2 #verified that this is the true sigma ie sqrt(sum(w*r^2)/df), where df = nsmp - ncol(D)
                        df = fit$df.residual
                        varbetahat <- diag(fit$cov.coefficients)[coef] * sigma2
                        tvals = betahat/sqrt(varbetahat)
                        pvals = 2 * pt(-abs(tvals), df)
                        p.adj <- p.adjust(pvals,method="BH")
                        de <- names(p.adj)[p.adj<=fdr]
                        posDE=NULL
                        posAll=NULL
                        if( !is.null(pIdx) )
                         {
                          posDE <- match(pIdx,de)
                          posDE[is.na(posDE)] <- Inf
                          posDE <- sort(posDE)
                          posAll <- match(pIdx,names(pvals))
                          posAll[is.na(posAll)] <- Inf
                          posAll <- sort(posAll)
                          posAll[which(posDE==Inf)] <- -1*posAll[which(posDE==Inf)]
                         }
                        pvals <- pvals[rownames(counts)]
                        return( list( "Pvals" = pvals , "PosDE" = posDE , "PosRank" = posAll ) )
                     }
            limmaFun <- function(counts, x1, x2, W, fdr, coef)
                     {
                        n = nrow(counts)
                        p = ncol(x1)
                        dge <- DGEList(counts=counts)
                        dge <- calcNormFactors(dge)
                        lib.size = dge$samples[,2]*dge$samples[,3] #TMM normalization
                        Y <- log2(t(counts+0.5)/(lib.size+1)*1e+06) #this is the transformation that voom uses
                        D = cbind(x1, x2, W)
                        fit=lmFit(t(Y),D,weights=NULL) #same as lm(transformed_counts~D, weights=weights)
                        fit2 <- ebayes(fit)
                        pvals <- fit2$p.value[,coef]
                        padj <- p.adjust(pvals,method="BH")
                        de <- names(padj[padj<=fdr])
                        posDE=NULL
                        posAll=NULL
                        if( !is.null(pIdx) )
                         {
                          posDE <- match(pIdx,de)
                          posDE[is.na(posDE)] <- Inf
                          posDE <- sort(posDE)
                          posAll <- match(pIdx,names(pvals))
                          posAll[is.na(posAll)] <- Inf
                          posAll <- sort(posAll)
                          posAll[which(posDE==Inf)] <- -1*posAll[which(posDE==Inf)]
                         }
                        pvals <- pvals[rownames(counts)]
                        return( list( "Pvals" = pvals , "PosDE" = posDE , "PosRank" = posAll ) )
                     }
            de2Fun <- function(counts, x1, x2, W, groups, fdr, coef)
                     {
                        D=cbind(x1,x2,W)
                        suppressMessages(dds <- DESeqDataSetFromMatrix(counts, DataFrame(groups), ~ groups))
                        suppressMessages(dds <- estimateSizeFactors(dds))
                        suppressMessages(dds <- estimateDispersions(dds, fitType = "parametric", quiet = TRUE, modelMatrix = D))
                        suppressMessages(dds <- nbinomWaldTest(dds, betaPrior = FALSE, quiet = TRUE, modelMatrix = D, useT=TRUE, df=ncol(counts)-ncol(D)))
                        wh.coef <- resultsNames(dds)[coef]
                        wh.stat <- paste0("WaldStatistic_",resultsNames(dds)[coef])
                        wh.SE <- paste0("SE_",resultsNames(dds)[coef])
                        b <- mcols(dds)[,wh.coef]
                        stats <- mcols(dds)[,wh.stat]
                        s2 <- mcols(dds)[,wh.SE]^2
                        varb.ev <- get_empirical_variances(s2,b,rescaleconst=0.688888)
                        wald <- stats
                        wald.ev <- b/sqrt(varb.ev)
                        df = nrow(D) - ncol(D)
                        pvalsZ <- 2*pnorm(abs(wald),lower.tail=FALSE)
                        pvalsT <- 2*pt(abs(wald),df=df,lower.tail=FALSE)
                        pvalsZ.ev <- 2*pnorm(abs(wald.ev),lower.tail=FALSE)
                        pvalsT.ev <- 2*pt(abs(wald.ev),df=df,lower.tail=FALSE)
                        names(pvalsZ) <- rownames(counts)
                        names(pvalsT) <- rownames(counts)
                        names(pvalsZ.ev) <- rownames(counts)
                        names(pvalsT.ev) <- rownames(counts)
                        return( list( "PvalsZ" = pvalsZ,  "PvalsT" = pvalsT,  "PvalsZ.ev" = pvalsZ.ev,  "PvalsT.ev" = pvalsT.ev ) )
                     }
            waldFun <- function(counts, x1, x2, W, groups, fdr, coef) #de2 does this already
                     {
                        D=cbind(x1,x2,W)
                        dds <- DESeqDataSetFromMatrix(counts, DataFrame(groups), ~ groups)
                        suppressMessages(dds <- estimateSizeFactors(dds))
                        suppressMessages(dds <- estimateDispersions(dds, fitType = "parametric", quiet = TRUE, modelMatrix = D))
                        suppressMessages(dds <- nbinomWaldTest(dds, betaPrior = FALSE, quiet = TRUE, modelMatrix = D, useT=TRUE, df=ncol(counts)-ncol(D)))
                        wh.coef <- paste0("SE_",resultsNames(dds)[coef])
                        s2 <- mcols(dds)[,wh.coef]^2
                        b <- mcols(dds)[,resultsNames(dds)[coef]]
                        varb.ev <- get_empirical_variances(s2,b,rescaleconst=0.688888)
                        wald.ev <- b/sqrt(varb.ev)
                        df = nrow(D) - ncol(D)
                        pvalsZ <- 2*pnorm(abs(wald.ev),lower.tail=FALSE)
                        pvalsT <- 2*pt(abs(wald.ev),df=df,lower.tail=FALSE)
                        names(pvalsZ) <- rownames(counts)
                        names(pvalsT) <- rownames(counts)
                        return( list( "PvalsZ" = pvalsZ , "PvalsT"=pvalsT ) )
                     }

            if( is.null(x2) ) x2 = matrix(1,nc=1,nr=ncol(counts))
            if( is.null(ncol(x1)) ) x1 = matrix(x1,nc=1,nr=length(x1))
            if( is.null(ncol(x2)) ) x2 = matrix(x2,nc=1,nr=length(x2))

            if( deMethod=="edgeR" )
             out <- eRFun(counts=counts, x1=x1, x2=x2, W=W, groups=groups, fdr=fdr, coef=coef)
            else if( deMethod=="voom" )
             out <- vFun(counts=counts, x1=x1, x2=x2, W=W, fdr=fdr, coef=coef)
            else if( deMethod=="empvar" )
             out <- evFun(counts=counts, x1=x1, x2=x2, W=W, fdr=fdr, coef=coef, use.weights=use.weights)
            else if( deMethod=="lm" )
             out <- lmFun(counts=counts, x1=x1, x2=x2, W=W, fdr=fdr, coef=coef)
            else if( deMethod=="limma" )
             out <- limmaFun(counts=counts, x1=x1, x2=x2, W=W, fdr=fdr, coef=coef)
            else if( deMethod=="deseq2" )
             out <- de2Fun(counts=counts, x1=x1, x2=x2, W=W, groups=groups, fdr=fdr, coef=coef)
            else if( deMethod=="wald" )
             out <- waldFun(counts=counts, x1=x1, x2=x2, W=W, groups=groups, fdr=fdr, coef=coef)
            else stop("not a valid 'deMethod' input")

            return(out)
          }


trcounts <- function( counts )
            {
              dge <- DGEList(counts=counts)
              dge <- calcNormFactors(dge)
              lib.size = dge$samples[,2]*dge$samples[,3] #TMM normalization
              counts <- log2(t(counts+0.5)/(lib.size+1)*1e+06) #this is the transformation that voom uses
              return(t(counts))
            }


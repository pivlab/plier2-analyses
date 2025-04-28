

mat_mult <- function(mat1, mat2) {
  is_fbm <- inherits(mat1, "FBM") 
  if (is_fbm ) {
    # For FBM objects, use the specific multiplication method
    return(big_prodMat(mat1, as.matrix(mat2)))
  } else {
    # Regular matrix multiplication
    return(mat1 %*% mat2)
  }
}

rotateSVD=function(svdres){
  upos=svdres$u
  uneg=svdres$u
  upos[upos<0]=0
  uneg[uneg>=0]=0
  uneg=-uneg
  sumposu=colSums(upos)
  sumnegu=colSums(uneg)
  
  
  for(i in 1:ncol(svdres$u)){
    if(sumnegu[i]>sumposu[i]){
      svdres$u[,i]=-svdres$u[,i] 
      svdres$v[,i]=-svdres$v[,i] 
    }
  }
  svdres
}


binarizeTop=function(Z, top, keepVals=T){
  for(i in 1:ncol(Z)){
    cutoff=sort(Z[,i],T)[top+1]
    
    
    if(cutoff==0){
      
      cutoff=min(Z[Z[,i]>0,i])
    }
    Z[Z[,i]<cutoff,i]=0
    if(!keepVals){
      Z[Z[,i]>0,i]=1
    }
  }
  Z
}


solveU=function(Z,  Chat=NULL, priorMat, penalty.factor,pathwaySelection="fast", alpha=0.9, maxPath=5,  nfolds=5, intercept=T, useSE=F, top=NULL, binary=FALSE, ...){
  if(nrow(Z)!=nrow(priorMat)){
    cm=commonRows(Z, priorMat)
    Z=Z[cm,]
    iim=match(cm, rownames(priorMat))
    priorMat=priorMat[iim,]
    show(dim(priorMat))
  }
  if(is.null(Chat)){
    Chat=pinv.ridge(crossprod(priorMat), 5)%*%(t(priorMat))
  }
  message("New solve U")
  Ur=Chat%*%Z #get U by OLS
  show(dim(Ur))
  Ur=apply(-Ur,2,rank) #rank
  Urm=apply(Ur,1,min)
  
  
  U=matrix(0,nrow=ncol(priorMat), ncol=ncol(Z))
  pathwaySelection=match.arg(pathwaySelection, c("fast", "complete"))
  
  if(!is.null(top)){
    Z=binarizeTop(Z, top, keepVals = TRUE)
  }
  if(pathwaySelection=="complete"){
    
    iip=which(Urm<=maxPath)
    message(paste("Picked", length(iip), "pathways"))
  }
  
  
  for(i in 1:ncol(Z)){
    
    if(pathwaySelection=="fast"){
      iip=which(Ur[,i]<=maxPath)
      
    }
    if(!binary){  
      
      set.seed(1); gres=cv.glmnet(y=Z[,i], x=priorMat[,iip], alpha=alpha, lower.limits=0, foldid =((1:nrow(Z)) %%nfolds)+1, 
                                  intercept=intercept, keep = T , nfolds = nfolds, standardize=F, dfmax=maxPath, ...)
      #plot(gres)
    }
    else{
      set.seed(1);gres=cv.glmnet(y=(Z[,i]>0)+1-1, x=priorMat[,iip], family="binomial", alpha=alpha, lower.limits=0, foldid =((1:nrow(Z)) %%nfolds)+1, intercept=intercept, keep = T , nfolds = nfolds,  standardize = T, type.measure="auc", dfmax=dfm, ...)
    }
    betaI=getNonZeroBetas(gres,se = useSE, index=T)
    beta=getNonZeroBetas(gres,se = useSE, index=F)
    
    U[iip[betaI],i]=as.vector(beta)
    
  }
  message(paste("Number of annotated columns is", sum(colSums(U)>0)))
  rownames(U)=substr(colnames(priorMat),1,30)
  colnames(U)=paste("LV", 1:ncol(U))
  return(U)
}

getChat=function(priorMat, method="fast"){
  if(method=="fast"){
    #this doesn't account for any covariance
    message("Presolving using dot method")
    L2n=sqrt(colSums(priorMat^2))
    Chat=t(sweep(priorMat,2,L2n, "/"))
    message("done")
  }
  else{
    message("Presolving using inverse method")
    Chat=pinv.ridge(crossprod(priorMat), 5)%*%(t(priorMat))
    message("done")
  }
  Chat
}






getMaxAUC=function(summary, verbose=F){
  max_auc_per_lv <- summary %>%
    group_by(`LV index`) %>%
    summarize(max_AUC = max(AUC, na.rm = TRUE)) %>%
    ungroup()
  if(verbose){
    message(paste("There are", sum(max_auc_per_lv$max_AUC>0.70), " LVs with AUC>0.70"))
    message(paste("There are", sum(max_auc_per_lv$max_AUC>0.90), " LVs with AUC>0.90"))
  }
  max_auc_per_lv
}

getAUCstats=function(summary){
  out=getMaxAUC(summary)
  unlist(lapply(c(0.7,0.8, 0.9), function(x){sum(out$max_AUC>x)}))
}






getMatchedPathwayMat=function(pathMat, new.genes, min.genes=10){
  cm=intersect(rownames(pathMat), new.genes)
  message("there are ", length(cm), " genes in the intersction between data and prior")
  matchPathMat=Matrix(0, nrow=length(new.genes), ncol=ncol(pathMat), sparse=T)
  rownames(matchPathMat)=new.genes
  colnames(matchPathMat)=colnames(pathMat)
  matchPathMat[cm,]=pathMat[cm,]
  
  genesInPath=colSums(matchPathMat)
  ii=which(genesInPath>=min.genes)
  message(paste("Removing", ncol(matchPathMat)-length(ii), "pathways"))
  
  return(matchPathMat[, ii])
  
}










AUC<-function(labels, values){
  posii=which(labels>0)
  negii=which(labels<=0)
  posn=length(posii)
  negn=length(negii)
  posval=values[posii]
  negval=values[negii]
  myres=list()
  if(posn>0&negn>0){
    res=wilcox.test(posval, negval, alternative="greater")
    myres$auc=(res$statistic)/(posn*negn)
    myres$pval=res$p.value
  }
  else{
    myres$auc=0.5
    myres$pval=NA
  }
  return(myres)
}


crossVal=function(plierRes,priorMat, priorMatcv){
  
  out=matrix(ncol=4, nrow=0)
  ii=which(colSums(plierRes$U)>0)
  Uauc=Matrix(0, nrow = nrow(plierRes$U), ncol=ncol(plierRes$U), sparse=T)
  Up=Matrix(0,nrow = nrow(plierRes$U), ncol=ncol(plierRes$U), sparse=T)
  
  
  for ( i in ii){ #for each column in U
    
    iipath=which(plierRes$U[,i]>0) #get the pathways
    
    if (length(iipath) > 1){ #more than one pathway
      for(j in iipath){
        iiheldout=which((rowSums(priorMat[,iipath, drop=F])==0) |(priorMat[,j]>0&priorMatcv[,j]==0))
        aucres=AUC(priorMat[iiheldout,j], plierRes$Z[iiheldout,i])
        out=rbind(out,c(colnames(priorMat)[j], i, aucres$auc, aucres$pval))
        Uauc[j,i]=aucres$auc
        Up[j,i]=-log10(aucres$pval)
        
      }}else{
        j <- iipath
        iiheldout=which((rowSums(matrix(priorMat[,iipath],ncol=1))==0) |(priorMat[,j]>0&priorMatcv[,j]==0))
        aucres=AUC(priorMat[iiheldout,j], plierRes$Z[iiheldout,i])
        out=rbind(out,c(colnames(priorMat)[j], i, aucres$auc, aucres$pval))
        Uauc[j,i]=aucres$auc
        Up[j,i]=-log10(aucres$pval)
      }#else
  }
  out=data.frame(out,stringsAsFactors = F)
  out[,3]=as.numeric(out[,3])
  out[,4]=as.numeric(out[,4])
  out[,5]=BH(out[,4])
  colnames(out)=c("pathway", "LV index", "AUC", "p-value", "FDR") 
  return(list(Uauc=Uauc, Upval=Up, summary=out))
}









getNonZeroBetas=function(gres,index=T, i=NULL, se=F,vector=F, intercept=F){
  if(is.null(i)){
    i=getBestIndex(gres, se=se)
    #  message(paste("index is",i))
  }
  
  ii=which(gres$glmnet.fit$beta[,i]!=0)
  if(index){
    ii
  }
  else if(vector){
    if(intercept){
      (c(gres$glmnet.fit$a0[i],gres$glmnet.fit$beta[,i, drop=F]))
    }
    else{
      gres$glmnet.fit$beta[,i, drop=F]
    }
  }
  else{
    
    gres$glmnet.fit$beta[ii,i, drop=F]
  }
}

getBestIndex=function(gres, se=F){
  if(!se){
    which(gres$lambda==gres$lambda.min)
  }
  else{
    which(gres$lambda==gres$lambda.1se)
  }
}


simpleDecomp=function(Y, k,svdres=NULL, L1=NULL, L2=NULL,
                      Zpos=T,max.iter=200, tol=5e-6, trace=F,
                      rseed=NULL, B=NULL, scale=1, pos.adj=3, adaptive.frac=0.05, adaptive.iter=30,  cutoff=0){
  
  
  #message("Checking type")
  # Detect matrix type
  is_fbm <- inherits(Y, "FBM")
  is_sparse <- inherits(Y, "dgCMatrix")
  
  
  ng=nrow(Y)
  ns=ncol(Y)
  
  Bdiff=Inf
  BdiffTrace=double()
  BdiffCount=0
  message("****")
  
  
  if (is.null(svdres)) {
    message("Computing SVD")
    if (is_fbm) {
      # For FBM, we need special handling for SVD
      set.seed(123)
      if (requireNamespace("bigstatsr", quietly = TRUE)) {
        # Use big_SVD from bigstatsr if available
        svdres <- bigstatsr::big_SVD(X = Y, k = k)
      } else {
        # Fallback: convert to regular matrix for SVD
        # This might be memory-intensive for large matrices
        svdres <- rsvd(to_matrix(Y), k = k)
      }
    } else if (is_sparse) {
      # For sparse matrices, use irlba or other sparse SVD methods
      if (requireNamespace("irlba", quietly = TRUE)) {
        set.seed(123)
        svdres <- irlba::irlba(Y, nv = k)
      } else {
        set.seed(123)
        svdres <- rsvd(Y, k = k)
      }
    } else {
      # Regular matrix
      set.seed(123)
      svdres <- rsvd(Y, k = k)
    }
    
    svdres <- rotateSVD(svdres)
  }
  

  
  if(is.null(L1)){
    L1=svdres$d[k]*scale
    if(!is.null(pos.adj)){
      L1=L1/pos.adj
    }
    
  }
  if(is.null(L2)){
    L2=svdres$d[k]*scale
  }
  L2k=L2*diag(k)
  #    L1=svdres$d[k]/2*scale
  print(paste0("L1 is set to ",L1))
  print(paste0("L2 is set to ",L2))
  
  if(is.null(B)){
    #initialize B with svd
    message("Init")
    B=t(svdres$v[, 1:k]%*%diag(sqrt(svdres$d[1:k])))
    #   B=t(svdres$v[1:ncol(Y), 1:k]%*%diag((svdres$d[1:k])))
    #   B=t(svdres$v[1:ncol(Y), 1:k])
    show(dim(B))
    message("Done")
  }
  else{
    message("B given")
  }
  
  
  
  
  
  if (!is.null(rseed)) {
    message("using random start")
    set.seed(rseed)
    B = t(apply(B, 1, sample))
    
  }
  
  
  round2=function(x){signif(x,4)}
  
  getT=function(x){-quantile(x[x<0], adaptive.frac)}
  
  pb <- txtProgressBar(min = 0, max = max.iter, style = 3)  # Set a high max value
 
  for ( i in 1:max.iter){
    setTxtProgressBar(pb, i) # Avoid exceeding 100
    #main loop    
    Zraw=Z=mat_mult(Y,t(B))%*%solve(tcrossprod(B)+L1*diag(k))
    
    if(i>=adaptive.iter && adaptive.frac>0){
      
      
      cutoffs=apply(Zraw,2, getT)
      
      for(j in 1:ncol(Z)){
        Z[Z[,j]<cutoffs[j],j]=0
      }
    }
    
    else if(Zpos){
      Z[Z<cutoff]=0
    }
    
    oldB=B
    
    
    if(is_fbm){
    ZYt=big_cprodMat(Y, as.matrix(Z))
    ZY=t(ZYt)
    B=solve(t(Z)%*%Z+L2k)%*%ZY
    }
    else{
    B=solve(t(Z)%*%Z+L2k)%*%mat_mult(t(Z),Y)
    }
    
    #update error
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)
 
    if(trace){
      message(paste0("iter",i, ))
    }
    
    #check for convergence
    if(i>52&&Bdiff>BdiffTrace[i-50]){
      BdiffCount=BdiffCount+1
    }
    else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }
    
    if(Bdiff<tol &&i>40){
      message(paste0("converged at  iteration ", i))
      break
    }
    if( BdiffCount>5&&i>40){
      message(paste0("stopped at  iteration ", i, " Bdiff is not decreasing"))
      break
    }
    
  }
  rownames(B)=colnames(Z)=paste("LV",1:k)
  return(list(B=B, Z=Z, Zraw=Zraw, L1=L1, L2=L2))
}


















PLIERv2=function(Y, priorMat,svdres=NULL, sdres=NULL,k=NULL, L1=NULL, L2=NULL, top=NULL, cvn=5, max.iter=350, trace=F, Chat=NULL, maxPath=5, doCrossval=T, penalty.factor=rep(1,ncol(priorMat)), glm_alpha=0.9, minGenes=10, tol=1e-6, seed=123456, allGenes=F, rseed=NULL, u.iter=20, max.U.updates=1, pathwaySelection=c("complete", "fast"), multiplier=1, adaptive.frac=0){
  
  getT=function(x){-quantile(x[x<0], adaptive.frac)}
  
  pathwaySelection=match.arg(pathwaySelection, c("complete", "fast"))
  
  
  
  
  message("**PLIER v2 **")
  
  
  # Detect matrix type
  is_fbm <- inherits(Y, "FBM")
  is_sparse <- inherits(Y, "dgCMatrix")
  
  
  if(nrow(priorMat)!=nrow(Y) || !all(rownames(priorMat)==rownames(data))){
    if(!allGenes){
      cm=commonRows(Y, priorMat)
      message(paste("Selecting common genes:", length(cm)))
      priorMat=priorMat[cm,]
      Y=Y[cm,]
    }
    else{
      extra.genes=setdiff(rownames(Y), rownames(priorMat))
      eMat=matrix(0, nrow=length(extra.genes), ncol=ncol(priorMat))
      rownames(eMat)=extra.genes
      priorMat=rbind(priorMat, eMat)
      priorMat=priorMat[rownames(Y),]
    }
    
  }
  numGenes=colSums(priorMat)
  
  heldOutGenes=list()
  iibad=which(numGenes<minGenes)
  if(length(iibad)>0){
    priorMat=priorMat[, -iibad]
    message(paste("Removed", length(iibad), "pathways with too few genes"))
  }
  if(doCrossval){
    
    
    priorMatCV=priorMat
    if(!is.null(seed))
      set.seed(seed)
    for(j in 1:ncol(priorMatCV)){
      
      iipos=which(priorMatCV[,j]>0)
      iiposs=sample(iipos, length(iipos)/5)
      priorMatCV[iiposs,j]=0
      heldOutGenes[[colnames(priorMat)[j]]]=rownames(priorMat)[iiposs]
      
    }
    C = priorMatCV
  }
  else{
    C=priorMat
  }
  
  nc=ncol(priorMat)
  ng=nrow(Y)
  ns=ncol(Y)
  
  Bdiff=-1
  BdiffTrace=double()
  BdiffCount=0
  if(is.null(Chat)){
    Cp=crossprod(C)
    Chat=pinv.ridge(crossprod(C), 5)%*%(t(C))
  }
  # YsqSum=sum(Y^2)
  #compute svd and use that as the starting point
  
  if(!is.null(svdres) && nrow(svdres$v)!=ncol(Y)){
    message("SVD V has the wrong number of columns")
    svdres=NULL
  }
  if(is.null(svdres) &&  is.null(sdres)){
    message("Computing SVD")
    if(ns>500){
      message("Using rsvd")
      set.seed(123456);svdres=rsvd(Y, k=min(ns, max(200, ns/4)), q=3)
    }
    else{
      svdres=svd(Y)
    }
    message("Done")
  }
  if(is.null(sdres)){
    if(is.null(k)){
      k=num.pc(svdres)*2
      k <- min(k, floor(ncol(Y)*0.9))
      message("k is set to ", k)
    }
    
    
    message("Simple decomp")
    if(is.null(sdres)){
      sdres=simpleDecomp(Y, k=k)
    }
  }
  else{
    message("usign SDres provided")
    
    if(nrow(Y)!=nrow(sdres$Z)){
      if(is.null(rownames(Y))|is.null(rownames(sdres$Z))){
        stop("Y and sdres$Z must have equal row numbers or row names")
      }    
      sdres$Z=sdres$Z[rownames(Y),]
    }
    u.iter=1
    k=ncol(sdres$Z)
  }
  
  Z=sdres$Z

  if(is.null(L1)){
    L1=sdres$L1
  }
  if(is.null(L2)){
    L2=sdres$L2
  }
  
  if(ncol(sdres$B)==ncol(Y)){
    B=sdres$B
  }
  
  oldB=B

  

  
  
  
  if(!is.null(rseed)){
    message("using random start")
    set.seed(rseed)
    B=t(apply(B, 1, sample))
    Z=apply(Z,2,sample)
  }
  
  
  
  
  
  
  U=matrix(0,nrow=ncol(C), ncol=k)
  
  
  round2=function(x){signif(x,4)}
  
  
  
  
  
  
  iter.full.start=iter.full=u.iter
  
  curfrac=0
  nposlast=Inf
  npos=-Inf
  num.U.updates=0
  L1k=L1*diag(k)
  L2k=L2*diag(k)
  
  if(is_fbm){
    ZYt=big_cprodMat(Y, as.matrix(Z))
    ZY=t(ZYt)
    B=solve(t(Z)%*%Z+L2k)%*%ZY
  }
  else{
    B=solve(t(Z)%*%Z+L2k)%*%mat_mult(t(Z),Y)
  }
  
  for ( iter in 1:max.iter){
    
    
    
    
    
    if(iter>=iter.full.start){
      
      
      
      
      
      if(iter==iter.full&&num.U.updates<max.U.updates ){ #update L3 to the target fraction
        #solveU=function(Z,  Chat, priorMat, penalty.factor,pathwaySelection="fast", glm_alpha=0.9, maxPath=10, nfolds=5){
        message(paste("Updating U at iteration", iter))
        U=solveU(Z, Chat, C, penalty.factor, pathwaySelection, glm_alpha, maxPath, binary=F, nfolds=cvn, top=top )
        
        num.U.updates=num.U.updates+1
        
        iter.full=iter.full+iter.full.start
        Z2=L1*C%*%U*multiplier
      }
      
      
      curfrac=(npos<-sum(apply(U,2,max)>0))/k
      #Z1=Y%*%t(B)
      Z1=mat_mult(Y, t(B))
      ratio=median((Z2/Z1)[Z2>0&Z1>0])
      Z=(Z1+Z2)%*%solve(tcrossprod(B)+L1k)
      
      
    }
    
    else{
      
      Z=mat_mult(Y,t(B))%*%solve(tcrossprod(B)+L1k)
    }
    
    
    
    if(adaptive.frac>0){
      
      
      cutoffs=apply(Z,2, getT)
      
      for(j in 1:ncol(Z)){
        Z[Z[,j]<cutoffs[j],j]=0
      }
    }
    else{
      
      Z[Z<0]=0
    }
    
    
    
    
    oldB=B

    
    if(is_fbm){
      ZYt=big_cprodMat(Y, as.matrix(Z))
      ZY=t(ZYt)
      B=solve(t(Z)%*%Z+L2k)%*%ZY
    }
    else{
      B=solve(t(Z)%*%Z+L2k)%*%mat_mult(t(Z),Y)
    }
    

    
    
    
    Bdiff=sum((B-oldB)^2)/sum(B^2)
    BdiffTrace=c(BdiffTrace, Bdiff)
    
    
    #err0=sum((Y-Z%*%B)^2)+sum((Z-C%*%U)^2)*L1+sum(B^2)*L2
    if(trace & iter >=iter.full.start){
      
      message(paste0("iter",iter,";pos. col. U=", sum(colSums(U)>0)))
    }
    else if (trace){
      message(paste0("iter",iter))
    }
    
    if(iter>52&&Bdiff>BdiffTrace[iter-50]){
      BdiffCount=BdiffCount+1
      message("Bdiff is not decreasing")
    }
    else if(BdiffCount>1){
      BdiffCount=BdiffCount-1
    }
    
    if(Bdiff<tol &iter>u.iter+10){
      message(paste0(Bdiff, "converged at  iteration ", iter))
      break
    }
    if( BdiffCount>5){
      message(paste0("converged at  iteration ", iter, " Bdiff is not decreasing"))
      break
    }
    
  }
  rownames(U)=colnames(priorMat)
  colnames(U)=rownames(B)=paste0("LV", 1:k)
  
  out=list( B=B, Z=Z, U=U, C=C, L1=L1, L2=L2, heldOutGenes=heldOutGenes)
  
  if(doCrossval){
    message("crossValidation")
    outAUC=crossVal(out, priorMat, priorMatCV)
    out$Uauc=outAUC$Uauc
    out$Up=outAUC$Upval
    out$summary=outAUC$summary
    out$priorMatCV=priorMatCV
    out$priorMat=priorMat
  }
  else{
    message("Not using cross-validation. No AUCs or p-values")
    
  }
  out$withPrior=which(colSums(out$U)>0)
  
  tt=apply(out$Uauc,2,max)
  message(paste("There are", sum(tt>0.70), " LVs with AUC>0.70"))
  message(paste("There are", sum(tt>0.90), " LVs with AUC>0.90"))
  #currently not working
  # rownames(out$B)=nameB(out)
  
  out
}


projectPLIER = function(PLIERres, newdata, scale=1) {
  stopifnot(nrow(PLIERres$Z) == nrow(newdata))
  
  # Check if newdata is a FBM/big.matrix object
  is_fbm = inherits(newdata, c("big.matrix", "FBM"))
  
  # Convert Matrix package matrices to standard R matrices for compatibility
  Z_matrix = if (inherits(PLIERres$Z, "Matrix")) as.matrix(PLIERres$Z) else PLIERres$Z
  
  if (is_fbm) {
    # FBM implementation - use big_cprodMat for efficient computation
    # But ensure Z is a standard matrix, not a Matrix package object
    ZYt = big_cprodMat(newdata, Z_matrix)
    ZY = t(ZYt)
  } else {
    # Standard matrix implementation - use regular matrix multiplication
    ZY = t(Z_matrix) %*% newdata
  }
  
  # Calculate the regularization matrix using standard matrix operations
  ZtZ = t(Z_matrix) %*% Z_matrix
  L2k = PLIERres$L2 * diag(ncol(Z_matrix))*scale
  
  # Solve the regularized system
  B = solve(ZtZ + L2k) %*% ZY
  
  return(B)
}

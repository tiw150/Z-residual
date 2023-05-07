cv_zresidual.coxph<- function( fit, data, nfolds=NULL,foldlist=NULL)
{
  # Required packages:
  if (!requireNamespace("pacman")) {
    install.packages("pacman")
  }
  pacman::p_load(
    "foreach",
    "parallel",
    "doParallel",
    "stringr"
  )
 
  mf <- model.frame(fit$formula, data)
  nc <- ncol (mf)
  fix_var<-mf[,-c(1,nc),drop=FALSE]
  
  if(is.null(nfolds))
  {
    ncores <- detectCores() 
    cl <- makeForkCluster(ncores)
    registerDoParallel(cl)
    nfolds<-10 %/% ncores*ncores
  }
  
  if(is.null(foldlist)) foldlist<-make_fold(fix_var=fix_var,y=mf[,nc],
                                            k=nfolds,censor=mf[,1][,2])
  
  res_fold <- as.list(1:nfolds)
  for(fid in 1:length (foldlist)) 
  {
    testcases <- foldlist[[fid]]
    data.train <- data[-testcases, ]
    data.test <- data[testcases, ]
    fit_traindata <- tryCatch(
      coxph(fit$formula,data=data.train),
      error = function(e) NA,
      warning = function(w) NA)

     if(any(is.na(fit_traindata))) {
       res_fold[[fid]]<-rep(NA,length(foldlist[[fid]]))
       attr(res_fold[[fid]], "SP") <-rep(NA,length(foldlist[[fid]]))
       attr(res_fold[[fid]], "hazfn") <-rep(NA,length(foldlist[[fid]]))
     }
    
    if(any(!is.na(fit_traindata))){
        res_fold[[fid]]<-zresidual.coxph(fit_coxph=fit_traindata,
                                         traindata=data.train,
                                         newdata = data.test)
      }
    
  }
  
  res_cv <- rep (0, nrow(data))
  attr(res_cv, "SP") <- rep(0, nrow(data))
  for(fid in 1:length (foldlist)) {
    testcases <- foldlist[[fid]]
    res_cv[testcases]<-res_fold[[fid]]
    attr(res_cv, "SP")[testcases] <- attr(res_fold[[fid]],"SP")
  }
  res_cv
}

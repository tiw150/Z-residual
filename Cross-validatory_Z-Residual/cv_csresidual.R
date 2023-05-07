cvcsresidual.coxph<- function( fit, data, nfolds=NULL,foldlist=NULL)
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
  
  #mm <- model.matrix(fit$formula, data)
  #mm_nc <- ncol (mm)
  #fix_var<-mm[,-c(1,mm_nc),drop=FALSE]
  
  if(is.null(nfolds))
  {
    ncores <- detectCores() 
    cl <- makeForkCluster(ncores)
    registerDoParallel(cl)
    nfolds<-10 %/% ncores*ncores
  }
  
  #if(is.null(foldlist)) foldlist<-kfold_fn(mf[,nc],k=nfolds)
  if(is.null(foldlist)) foldlist<-make_fold(fix_var=fix_var,y=mf[,nc],
                                            k=nfolds,censor=mf[,1][,2])
  
  res_fold <- as.list(1:nfolds)
  for(fid in 1:length (foldlist)) 
  {
    testcases <- foldlist[[fid]]
    data.train <- data[-testcases, ]
    data.test <- data[testcases, ]
    
    # execution_limit <- 1
    # exe_i <- 1
    # while(exe_i <= execution_limit)
    # {
    #   cat("(exe_i = ", exe_i," / ", execution_limit,  ")\n", sep="")
    #   Sys.sleep(0.1)
    # 
    #   fit_traindata <- tryCatch(
    #     coxph(fit$formula,data=data.train),
    #     error = function(e) NA,
    #     warning = function(w) NA)
    # 
    #   if(any(is.na(fit_traindata))) {
    #     exe_i <- exe_i + 1
    #   }else{
    #     exe_i<- execution_limit+1
    #   }
    # }
    
    fit_traindata <- tryCatch(
      coxph(fit$formula,data=data.train),
      error = function(e) NA,
      warning = function(w) NA)
    
    if(any(is.na(fit_traindata))) {
      res_fold[[fid]]<-rep(NA,length(foldlist[[fid]]))
      attr(res_fold[[fid]], "SP") <-rep(NA,length(foldlist[[fid]]))
      attr(res_fold[[fid]], "z_hat_new") <-rep(NA,length(foldlist[[fid]]))
    }
    
    if(any(!is.na(fit_traindata))){
      res_fold[[fid]]<-csresidual.coxph(fit_coxph=fit_traindata,
                                        traindata=data.train,
                                        newdata = data.test)
    }
    
  }
  
  res_cv <- rep (0, nrow(data))
  for(fid in 1:length (foldlist)) {
    testcases <- foldlist[[fid]]
    res_cv[testcases]<-res_fold[[fid]]
  }
  res_cv
}

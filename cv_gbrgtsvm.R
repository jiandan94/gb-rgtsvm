cv_gbrgtsvm <- function(x, y, tuning, kernel = "linear", kfolds = 5){
  if(!require(foreach) | !require(doParallel)){
    cat("Warning message:\nPlease install packages: foreach, doParallel!\n")
  }
  
  cv_main <- function(x, y, lam, a, lam1, lam2, kernel = "linear", ...){
    x_pos <- x[which(y == 1), ]
    x_neg <- x[which(y == -1), ]
    
    cvindex_pos <- sample(nrow(x_pos), nrow(x_pos))
    cvstep_pos <- floor(nrow(x_pos)/kfolds)
    cvindex_neg <- sample(nrow(x_neg), nrow(x_neg))
    cvstep_neg <- floor(nrow(x_neg)/kfolds)
    gbrgtsvm_acc <- rep(0, kfolds)
    
    for(i in 1:kfolds){
      a_pos <- 1 + cvstep_pos*(i - 1)
      b_pos <- a_pos + cvstep_pos - 1
      a_neg <- 1 + cvstep_neg*(i - 1)
      b_neg <- a_neg + cvstep_neg - 1
      
      trainpos_index <- cvindex_pos[-c(a_pos:b_pos)]
      trainneg_index <- cvindex_neg[-c(a_neg:b_neg)]
      testpos_index <- cvindex_pos[c(a_pos:b_pos)]
      testneg_index <- cvindex_neg[c(a_neg:b_neg)]
      
      trainx_pos <- x_pos[trainpos_index,]
      trainx_neg <- x_neg[trainneg_index,]
      testx_pos <- x_pos[testpos_index,]
      testx_neg <- x_neg[testneg_index,]
      
      trainx <- rbind(trainx_pos, trainx_neg)
      trainy <- c(rep(1, nrow(trainx_pos)), rep(-1, nrow(trainx_neg)))
      testx <- rbind(testx_pos, testx_neg)
      testy <- c(rep(1, nrow(testx_pos)), rep(-1, nrow(testx_neg)))
      
      gbrgtsvm_fit <- gbrgtsvm(x = trainx, y = trainy, lam = lam, a = a, 
                               lam1 = lam1, lam2 = lam2, 
                               kernel = kernel, ...)
      gbrgtsvm_wavg <- gbrgtsvm_fit$wavg
      gbrgtsvm_bavg <- gbrgtsvm_fit$bavg
      
      if(kernel == "linear"){
        fity <- sign(testx%*%gbrgtsvm_wavg + gbrgtsvm_bavg)
      }else{
        fity <- sign(gaussian_kernel(testx, t(gbrgtsvm_fit$xk), ...)%*%gbrgtsvm_wavg + gbrgtsvm_bavg)
      }

      gbrgtsvm_acc[i] <- sum(diag(table(testy, fity)))/sum(table(testy, fity))
    }
    
    cv_result <- c(mean(gbrgtsvm_acc), sd(gbrgtsvm_acc), lam, a, lam1, lam2, ...)
    return(cv_result)
  }
  
  core_num <- detectCores()
  cl <- makeCluster(core_num - 1)
  registerDoParallel(cl)
  
  if(kernel == "linear"){
    cv_all <- foreach(lam=tuning[,1], 
                      a=tuning[,2],
                      lam1=tuning[,3],
                      lam2=tuning[,4],
                      .export = "gbrgtsvm",
                      .errorhandling = "remove",
                      .combine = "rbind") %dopar% cv_main(x, y, lam, a, lam1, lam2, kernel = "linear")
  }
  else if(kernel == "gaussian"){
    cv_all <- foreach(lam=tuning[,1], 
                      a=tuning[,2], 
                      lam1=tuning[,3],
                      lam2=tuning[,4],
                      sigma=tuning[,5],
                      .export = c("gbrgtsvm", "gaussian_kernel"),
                      .errorhandling = "remove",
                      .combine = "rbind") %dopar% cv_main(x, y, lam, a, lam1, lam2, kernel = "gaussian", sigma)
  }
  else{
    cat("Warning message:\n'kenerl' should be linear or gaussian!")
  }
  
  stopImplicitCluster()
  stopCluster(cl)
  
  optimal_index <- which.max(cv_all[,1])# find the optimal parameters with maximum prediction accuracy
  optimal_tuning <- cv_all[optimal_index,]
  return(optimal_tuning)
}

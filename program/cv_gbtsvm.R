cv_gbtsvm <- function(x, y, tuning, kernel = "linear", kfolds = 5){
  if(!require(foreach) | !require(doParallel)){
    cat("Warning message:\nPlease install packages: foreach, doParallel!\n")
  }
  
  cv_main <- function(x, y, lam1, lam2, kernel = "linear", ...){
    x_pos <- x[which(y == 1), ]
    x_neg <- x[which(y == -1), ]
    
    cvindex_pos <- sample(nrow(x_pos), nrow(x_pos))
    cvstep_pos <- floor(nrow(x_pos)/kfolds)
    cvindex_neg <- sample(nrow(x_neg), nrow(x_neg))
    cvstep_neg <- floor(nrow(x_neg)/kfolds)
    gbtsvm_acc <- rep(0, kfolds)
    
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
      
      gbtsvm_fit <- gbtsvm(x = trainx, y = trainy, 
                           lam1 = lam1, lam2 = lam2, 
                           kernel = kernel, ...)
      gbtsvm_wpos <- gbtsvm_fit$wpos
      gbtsvm_bpos <- gbtsvm_fit$bpos
      gbtsvm_wneg <- gbtsvm_fit$wneg
      gbtsvm_bneg <- gbtsvm_fit$bneg
      
      if(kernel == "linear"){
        dist_pos <- abs(testx%*%gbtsvm_wpos + gbtsvm_bpos)/sqrt(sum(gbtsvm_wpos^2))
        dist_neg <- abs(testx%*%gbtsvm_wneg + gbtsvm_bneg)/sqrt(sum(gbtsvm_wneg^2))
      }else{
        dist_pos <- abs(gaussian_kernel(testx, t(gbtsvm_fit$xk), ...)%*%gbtsvm_wpos + gbtsvm_bpos)/sqrt(sum(gbtsvm_wpos^2))
        dist_neg <- abs(gaussian_kernel(testx, t(gbtsvm_fit$xk), ...)%*%gbtsvm_wneg + gbtsvm_bneg)/sqrt(sum(gbtsvm_wneg^2))
      }
      
      fity <- sign(dist_neg - dist_pos)
      gbtsvm_acc[i] <- sum(diag(table(testy, fity)))/sum(table(testy, fity))
    }
    
    cv_result <- c(mean(gbtsvm_acc), sd(gbtsvm_acc), lam1, lam2, ...)
    return(cv_result)
  }
  
  core_num <- detectCores()
  cl <- makeCluster(core_num - 1)
  registerDoParallel(cl)
  
  if(kernel == "linear"){
    cv_all <- foreach(lam1=tuning[,1], 
                      lam2=tuning[,2],
                      .export = "gbtsvm",
                      .errorhandling = "remove",
                      .combine = "rbind") %dopar% cv_main(x, y, lam1, lam2, kernel = "linear")
  }
  else if(kernel == "gaussian"){
    cv_all <- foreach(lam1=tuning[,1], 
                      lam2=tuning[,2], 
                      sigma=tuning[,3],
                      .errorhandling = "remove",
                      .export = c("gbtsvm", "gaussian_kernel"),
                      .combine = "rbind") %dopar% cv_main(x, y, lam1, lam2, kernel = "gaussian", sigma)
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

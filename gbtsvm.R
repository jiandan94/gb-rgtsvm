gbtsvm <- function(x, y, lam1, lam2, kernel = "linear", ...){
  library(MASS)
  
  ## main functions of GB generation
  gb_ball <- function(x, y){
    
    r_fun <- function(x, center){
      if(is.matrix(x)){
        n <- nrow(x)
        center_mat <- matrix(rep(center, n), byrow = T, nrow = n)
        r <- mean(sqrt(apply((x - center_mat)^2, 1, sum)))
      }else{# if a point is a GB
        r <- sqrt(sum((x - center)^2))
      }
      return(r)
    }
    
    class_fun <- function(gb_ind, y){
      ind_y <- y[gb_ind]
      return(sign(sum(ind_y)))
    }
    
    purity_fun <- function(gb_ind, y){
      ind_y <- y[gb_ind]
      n_pos <- sum(ind_y == 1)
      n_neg <- sum(ind_y == -1)
      return(max(n_pos, n_neg)/length(ind_y))
    }
    
    if(nrow(x) > 2){
      # GB generation
      gb  <- kmeans(x, centers = 2)
      gb_class <- gb$cluster
      
      gb_x1 <- x[which(gb_class == 1),]
      gb_x2 <- x[which(gb_class == 2),]
      
      if(is.matrix(gb_x1)){
        gb_c1 <- apply(gb_x1, 2, mean)
      }else{# if a point is a GB
        gb_c1 <- gb_x1
      }
      
      if(is.matrix(gb_x2)){
        gb_c2 <- apply(gb_x2, 2, mean)
      }else{# if a point is a GB
        gb_c2 <- gb_x2
      }
      
      gb_r1 <- r_fun(gb_x1, gb_c1)
      gb_r2 <- r_fun(gb_x2, gb_c2)
      
      gb_pur1 <- purity_fun(which(gb_class == 1), y)
      gb_pur2 <- purity_fun(which(gb_class == 2), y)
      
      gb_ball1 <- list(x = gb_x1, y = y[which(gb_class == 1)], purity = gb_pur1, 
                       center = gb_c1, r = gb_r1, class = class_fun(which(gb_class == 1), y))
      gb_ball2 <- list(x = gb_x2, y = y[which(gb_class == 2)], purity = gb_pur2,
                       center = gb_c2, r = gb_r2, class = class_fun(which(gb_class == 2), y))
    }else{# two points with opposite labels
      gb_x1 <- x[1,]
      gb_x2 <- x[2,]
      
      gb_c1 <- x[1,]
      gb_c2 <- x[2,]
      
      gb_r1 <- 0
      gb_r2 <- 0
      
      gb_pur1 <- 1
      gb_pur2 <- 1
      
      gb_ball1 <- list(x = gb_x1, y = y[1], purity = gb_pur1, 
                       center = gb_c1, r = gb_r1, class = y[1])
      gb_ball2 <- list(x = gb_x2, y = y[2], purity = gb_pur2,
                       center = gb_c2, r = gb_r2, class = y[2])
    }
    
    return(list(gb_ball1, gb_ball2))
  }
  
  gb_fun <- function(x, y, pur_level){
    # initialization
    gb_list <- gb_ball(x, y)
    
    pur_vec <- unlist(lapply(gb_list, function(u){return(u$purity)}))
    active_set <- which(pur_vec < pur_level)
    
    while (length(active_set) >= 1) {
      
      for(i in active_set){
        if(nrow(gb_list[[i]]$x) <= 10){# for the abnormal case
          gb_list[[i]]$purity <- 1
          if(gb_list[[i]]$class == 0){
            gb_list[[i]]$class <- 1
          }
          gb_list[[i]] <- list(gb_list[[i]])
        }else{
          x_sub <- gb_list[[i]]$x
          y_sub <- gb_list[[i]]$y
          gb_list[[i]] <- gb_ball(x_sub, y_sub)
        }
      }
      
      gb_list <- c(unlist(gb_list[active_set], recursive = FALSE), gb_list[-active_set])
      
      pur_vec <- unlist(lapply(gb_list, function(u){return(u$purity)}))
      
      active_set <- which(pur_vec < pur_level)
    }
    
    gb_center_mat <- matrix(unlist(lapply(gb_list, function(u){return(u$center)})), 
                            byrow = T, ncol = ncol(x))
    gb_r_vec <- unlist(lapply(gb_list, function(u){return(u$r)}))
    gb_class_vec <- unlist(lapply(gb_list, function(u){return(u$class)}))
    
    return(list(cmat = gb_center_mat, rvec = gb_r_vec, classvec = gb_class_vec))
  }
  
  ## this problem is solved by clipDCD algorithm (Peng, 2014)
  clipDCD <- function(H, l, u, LB, UB, eps=0.01){# the main part of clipDCD
    StopValue <- 1
    while(StopValue >= eps){
      # find the feasible index set 
      f1 <- l - H%*%u
      f2 <- diag(H)
      fvalue <- f1/f2
      index1 <- which(fvalue < 0)
      index2 <- which(u > LB)
      index3 <- which(fvalue > 0)
      index4 <- which(u < UB)
      FeasibleIndexSet <- c(intersect(index1, index2), intersect(index3, index4))
      if(length(FeasibleIndexSet) <= 0){
        break
      }
      
      # find the index L and compute Lambda
      ObjectValue <- (f1^2/f2)[FeasibleIndexSet]
      OptimalIndex <- FeasibleIndexSet[which.max(ObjectValue)]
      Lambda <- fvalue[OptimalIndex]
      
      # update u[OptimalIndex]
      u[OptimalIndex] <- min(max(LB[OptimalIndex], u[OptimalIndex] + Lambda), UB[OptimalIndex])
      
      # compute the StopVlue
      StopValue <- (f1^2/f2)[OptimalIndex]
    }
    return(u)
  }
  
  ## find GBs
  if(kernel == "gaussian"){
    if(nrow(x) > 100){
      xk <- x[sample(nrow(x), 100),]
    }else{
      xk <- x
    }
    x <- gaussian_kernel(x, t(xk), ...)
  }
  # GB generation
  gb_result <- gb_fun(x, y, pur_level = 0.95)
  # center matrix
  x <- gb_result$cmat
  y <- gb_result$classvec
  r <- gb_result$rvec
  
  r_pos <- r[which(y == 1)]
  r_neg <- r[which(y == -1)]

  x_pos <- x[which(y == 1), ]
  x_neg <- x[which(y == -1), ]
  
  n_pos <- nrow(x_pos)
  n_neg <- nrow(x_neg)
  p <- ncol(x_pos)
  
  e_pos <- matrix(1, nrow = n_pos, ncol = 1)
  e_neg <- matrix(1, nrow = n_neg, ncol = 1)
  
  H <- cbind(x_pos, e_pos)
  G <- cbind(x_neg, e_neg)
  
  ## compute w+ and b+
  Q <- G%*%ginv(t(H)%*%H)%*%t(G)
  l <- e_neg + r_neg
  
  LB <- rep(0, n_neg)
  UB <- rep(lam1, n_neg)
  
  u <- 0.5*(LB + UB)
  
  u <- clipDCD(Q, l, u, LB, UB)
  
  w_tilde_pos <- -ginv(t(H)%*%H)%*%t(G)%*%u
  w_pos <- w_tilde_pos[1:p]
  b_pos <- w_tilde_pos[p + 1]
  
  ## compute w- and b-
  Q <- H%*%ginv(t(G)%*%G)%*%t(H)
  u <- rep(0, n_pos)
  
  l <- e_pos + r_pos
  
  LB <- rep(0, n_pos)
  UB <- rep(lam2, n_pos)
  
  u <- 0.5*(LB + UB)
  
  u <- clipDCD(Q, l, u, LB, UB)
  
  w_tilde_neg <- ginv(t(G)%*%G)%*%t(H)%*%u
  w_neg <- w_tilde_neg[1:p]
  b_neg <- w_tilde_neg[p + 1]
  
  wpos_norm <- norm(w_pos, type = "2")
  wneg_norm <- norm(w_neg, type = "2")
  
  
  gbtsvmList <- list(wpos = w_pos/wpos_norm, bpos = b_pos/wpos_norm, 
                     wneg = w_neg/wneg_norm, bneg = b_neg/wneg_norm,
                     wavg = w_pos/wpos_norm + w_neg/wneg_norm, 
                     bavg = b_pos/wpos_norm + b_neg/wneg_norm)
  if(kernel == "gaussian"){
    gbtsvmList$xk <- xk
  }
  return(gbtsvmList)
}

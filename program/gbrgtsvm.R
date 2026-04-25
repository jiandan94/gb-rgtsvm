gbrgtsvm <- function(x, y, lam, a, lam1, lam2, kernel = "linear", ...){
  require(MASS)
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
  
  ## prox function
  prox <- function(tw){
    return(HH%*%tw)
  }
  
  ## d_f function
  d_f <- function(tw){
    d_lR <- function(u){
      return(ifelse(u > 0, lam*(a^2)*u*exp(-a*u), 0))
    }
    
    df1 <- colSums(diag(d_lR(r_pos - c(A%*%tw)))%*%A) 
    df2 <- colSums(diag(d_lR(r_neg - c(B%*%tw)))%*%B) 
    
    df <- -(lam1*df1 + lam2*df2)
    
    return(df)
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
  C <- gb_result$cmat
  y <- gb_result$classvec
  r <- gb_result$rvec + 1
  
  r_pos <- r[which(y == 1)]
  r_neg <- r[which(y == -1)]
  
  c_pos <- C[which(y == 1), ]
  c_neg <- C[which(y == -1), ]
  
  n_pos <- nrow(c_pos)
  n_neg <- nrow(c_neg)
  
  p <- ncol(c_pos)
  
  tc_pos <- cbind(c_pos, matrix(1, n_pos, 1))
  tc_neg <- cbind(c_neg, matrix(1, n_neg, 1))
  
  # compute H, A, B
  H11 <- t(tc_pos)%*%tc_pos
  H12 <- matrix(0, p+1, p+1)
  H22 <- t(tc_neg)%*%tc_neg
  H <- rbind(cbind(H11, H12), cbind(H12, H22))
  
  A <- cbind(tc_pos, tc_pos)
  
  B <- -cbind(tc_neg, tc_neg)
  
  # initialization
  b0 <- 0.75
  
  lip_const1 <- lam1*sum(sqrt(apply(A^2, 1, sum)))
  lip_const2 <- lam2*sum(sqrt(apply(B^2, 1, sum)))
  lip_const <- (lip_const1 + lip_const2)*(lam*a/exp(1))
  
  a0 <- 2*(1 - b0)/lip_const
  
  HH <- ginv(a0*H + (a0+1)*diag(nrow(H)))
  
  # update
  # main iteration
  w0 <- rep(1, 2*p+2)
  w1 <- w0 + 0.5
  t <- 0
  while(sum(abs(w0)) > 0 && sum(abs(w1-w0))/sum(abs(w0)) >= 0.001 && t <= 5000){
    ww <- w1
    w1 <- prox(w1 - a0*d_f(w1) + b0*(w1 - w0))
    w0 <- ww
    t <- t + 1
  }
  
  w_pos <- w1[1:p]
  b_pos <- w1[p+1]
  w_neg <- w1[(p+2):(2*p+1)]
  b_neg <- w1[2*p+2]
  
  wpos_norm <- norm(w_pos, type = "2")
  wneg_norm <- norm(w_neg, type = "2")
  
  gbrgtsvmList <- list(wpos = w_pos/wpos_norm, bpos = b_pos/wpos_norm, 
                       wneg = w_neg/wneg_norm, bneg = b_neg/wneg_norm,
                       wavg = w_pos/wpos_norm + w_neg/wneg_norm, 
                       bavg = b_pos/wpos_norm + b_neg/wneg_norm)
  if(kernel == "gaussian"){
    gbrgtsvmList$xk <- xk
  }
  return(gbrgtsvmList)
}

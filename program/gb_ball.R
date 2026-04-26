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

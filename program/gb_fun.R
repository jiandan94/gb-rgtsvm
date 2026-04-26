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

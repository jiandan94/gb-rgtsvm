gb_plot <- function(x, y, purity = 0.9){
  library(ggplot2)
  
  plt_dt <- as.data.frame(cbind(x, y))
  colnames(plt_dt) <- c("x1", "x2", "class")
  plt_dt$class <- as.factor(plt_dt$class)
  
  # gb generation
  gb_result <- gb_fun(x, y, pur_level = purity)
  
  gb_dt <- as.data.frame(cbind(gb_result$cmat, gb_result$rvec, gb_result$classvec))
  colnames(gb_dt) <- c(paste0("x", 1:ncol(x)), "r", "class")
  gb_dt$class <- as.factor(gb_dt$class)
  
  ## plot GBs
  plt <- ggplot(data = plt_dt, aes(x = x1, y = x2, colour = class)) + 
    geom_point(aes(shape = class), size = 2) + 
    geom_point(data = gb_dt, size = 3, shape = 8) +
    scale_colour_manual(values = c("blue", "red"))+
    theme_bw() + 
    coord_fixed()
  
  # plot the positive GBs
  gb_pos <- subset(gb_dt, class == 1 & r > 0)
  for (i in 1:nrow(gb_pos)) {
    circle_data <- circle_dt_fun(center = gb_pos[i, 1:2], radius = gb_pos$r[i])
    plt <- plt + geom_path(data = circle_data, linewidth = 0.4, colour = "red")
  }
  
  # plot the negative GBs
  gb_neg <- subset(gb_dt, class == -1 & r > 0)
  for (i in 1:nrow(gb_neg)) {
    circle_data <- circle_dt_fun(center = gb_neg[i, 1:2], radius = gb_neg$r[i])
    plt <- plt + geom_path(data = circle_data, linewidth = 0.4, colour = "blue")
  }
  
  return(plt)
}

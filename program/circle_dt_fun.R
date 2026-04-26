circle_dt_fun <- function(center, radius){
  
  # generate sequence
  theta <- seq(0, 2*pi, length.out = 1000)
  
  circle_data <- data.frame(
    x1 = center[1,1] + radius * cos(theta),  # x coordinate：x0 + r*cosθ
    x2 = center[1,2] + radius * sin(theta)   # y coordinate：y0 + r*sinθ
  )
  
  return(circle_data)
}
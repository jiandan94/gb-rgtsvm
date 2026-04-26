rm(list = ls())

library(ggplot2)
library(reshape2)
library(MASS)

source("program/gb_ball.R")
source("program/gb_fun.R")
source("program/gb_plot.R")
source("program/circle_dt_fun.R")

source("program/gbtsvm.R")
source("program/cv_gbtsvm.R")

source("program/gbrgtsvm.R")
source("program/cv_gbrgtsvm.R")

source("program/ParameterGrid.R")
source("program/gaussian_kernel.R")

# generate samples
set.seed(258)

n <- 400
xpos_x1 <- runif(0.5*n, -2, 2)
xpos_x2 <- runif(0.5*n, 0.6*sin(3*xpos_x1) + 0.25, 0.6*sin(3*xpos_x1) + 1.4) + 0.2*runif(0.5*n, -1, 1)
xpos <- cbind(xpos_x1, xpos_x2)
xneg_x1 <- runif(0.5*n, -2, 2)
xneg_x2 <- runif(0.5*n, 0.6*sin(2.95*xneg_x1 + 0.4) - 1.35, 0.6*sin(2.95*xneg_x1 + 0.4) - 0.25) + 0.2*runif(0.5*n, -1, 1)
xneg <- cbind(xneg_x1, xneg_x2)
x <- rbind(xpos, xneg)
y <- rep(c(1,-1), c(0.5*n, 0.5*n))

# add noise
noise_ratio <- 0.10
noise_ind <- sample(0.5*n, floor(noise_ratio*0.5*n))
y[noise_ind] <- -1

plt_dt <- as.data.frame(cbind(x, y))
colnames(plt_dt) <- c("x1", "x2", "class")
plt_dt$class <- as.factor(plt_dt$class)

ggplot(data = plt_dt, aes(x = x1, y = x2, colour = class)) + 
  geom_point(aes(shape = class), size = 2) + 
  scale_colour_manual(values = c("blue", "red"))+
  theme_bw() + 
  coord_fixed() + 
  theme(legend.position = c(0.9, 0.8))

# gb generation
plt <- gb_plot(x, y, purity = 0.95)

plt + 
  theme(legend.position = c(0.9, 0.8))

## gbtsvm
TuningSeq <- 2^seq(-8,8,2)
gbtsvm_tuning <- ParameterGrid(list(lam1 = TuningSeq, lam2 = TuningSeq, sigma = TuningSeq))

cvgbtsvm <- cv_gbtsvm(x, y, gbtsvm_tuning, kernel = "gaussian")
gbtsvm_op_tuning <- cvgbtsvm[-c(1,2)]

op_gbtsvm <- gbtsvm(x, y, lam1=gbtsvm_op_tuning[1], lam2=gbtsvm_op_tuning[2],
                    kernel = "gaussian", sigma = gbtsvm_op_tuning[3])

## gbrgtsvm
TuningSeq <- 2^seq(-8,8,2)
gbrgtsvm_tuning <- ParameterGrid(list(lam = c(1,3,5), a = 1, lam1 = TuningSeq, 
                                      lam2 = TuningSeq, sigma = TuningSeq))

cvgbrgtsvm <- cv_gbrgtsvm(x, y, gbrgtsvm_tuning, kernel = "gaussian")
gbrgtsvm_op_tuning <- cvgbrgtsvm[-c(1,2)]

op_gbrgtsvm <- gbrgtsvm(x, y, lam=gbrgtsvm_op_tuning[1], a=gbrgtsvm_op_tuning[2], 
                        lam1=gbrgtsvm_op_tuning[3], lam2=gbrgtsvm_op_tuning[4],
                        kernel = "gaussian", sigma = gbrgtsvm_op_tuning[5])
op_gbtsvm <- gbtsvm(x, y, lam1=gbtsvm_op_tuning[1], lam2=gbtsvm_op_tuning[2],
                    kernel = "gaussian", sigma = gbtsvm_op_tuning[3])

## plot
plot_x1 <- seq(-2.5, 2.5, length.out = 100)
plot_x2 <- seq(-2.5, 2.5, length.out = 100)

gbtsvm_avg <- function(s1, s2){
  s <- matrix(c(s1, s2), ncol = 2)
  fity <- gaussian_kernel(s, t(x), gbtsvm_op_tuning[3])%*%op_gbtsvm$wavg + op_gbtsvm$bavg
  return(fity)
}
z_avg <- outer(plot_x1, plot_x2, gbtsvm_avg)
colnames(z_avg) <- plot_x2
rownames(z_avg) <- plot_x1
gbtsvm_data_avg <- melt(z_avg)
colnames(gbtsvm_data_avg) <- c("x1", "x2", "value")

gbtsvm_data_avg$class <- as.factor(rep(1, nrow(gbtsvm_data_avg)))

gbrgtsvm_avg <- function(s1, s2){
  s <- matrix(c(s1, s2), ncol = 2)
  fity <- gaussian_kernel(s, t(op_gbrgtsvm$xk), gbrgtsvm_op_tuning[5])%*%op_gbrgtsvm$wavg + op_gbrgtsvm$bavg
  return(fity)
}
z_avg <- outer(plot_x1, plot_x2, gbrgtsvm_avg)
colnames(z_avg) <- plot_x2
rownames(z_avg) <- plot_x1
gbrgtsvm_data_avg <- melt(z_avg)
colnames(gbrgtsvm_data_avg) <- c("x1", "x2", "value")

gbrgtsvm_data_avg$class <- as.factor(rep(1, nrow(gbrgtsvm_data_avg)))

plt + 
  stat_contour(data = gbtsvm_data_avg, aes(z = value), 
               breaks = 0, linewidth = 1, 
               colour = "black", linetype = "dashed") +
  stat_contour(data = gbrgtsvm_data_avg, aes(z = value), 
               breaks = 0, linewidth = 1, 
               colour = "black") +
  scale_x_continuous(limits = c(-2.1, 2.1)) + 
  scale_y_continuous(limits = c(-2, 2)) +
  theme(legend.position = c(0.9, 0.8))


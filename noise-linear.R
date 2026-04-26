rm(list = ls())

library(ggplot2)
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

# generate samples
n <- 400

set.seed(425)

sig <- diag(c(0.02, 0.04))
xpos1 <- mvrnorm(n/4, mu = c(0.4, 0.9), Sigma = sig)
xpos2 <- mvrnorm(n/4, mu = c(-0.3, 0.9), Sigma = sig)
xneg1 <- mvrnorm(n/4, mu = c(-0.7, 0.2), Sigma = sig)
xneg2 <- mvrnorm(n/4, mu = c(0.1, 0.2), Sigma = sig)
xpos <- rbind(xpos1, xpos2)
xneg <- rbind(xneg1, xneg2)
x <- rbind(xpos1, xpos2, xneg1, xneg2)
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
  theme(legend.position = c(0.1, 0.8))
ggsave("noise10-rawdata.eps", width=5)

# gb generation
plt <- gb_plot(x, y, purity = 0.95)

plt + 
  theme(legend.position = c(0.1, 0.8))
ggsave("noise10-granulation.eps", width=5)

## gbtsvm
TuningSeq <- 2^seq(-8,8,2)
gbtsvm_tuning <- ParameterGrid(list(lam1 = TuningSeq, lam2 = TuningSeq))

cvgbtsvm <- cv_gbtsvm(x, y, gbtsvm_tuning, kernel = "linear")
gbtsvm_op_tuning <- cvgbtsvm[-c(1,2)]

op_gbtsvm <- gbtsvm(x, y, lam1=gbtsvm_op_tuning[1], lam2=gbtsvm_op_tuning[2])

## gbrgtsvm
TuningSeq <- 2^seq(-8,8,2)
gbrgtsvm_tuning <- ParameterGrid(list(lam = c(1,3,5), a = 1, lam1 = TuningSeq, lam2 = TuningSeq))

cvgbrgtsvm <- cv_gbrgtsvm(x, y, gbrgtsvm_tuning, kernel = "linear")
gbrgtsvm_op_tuning <- cvgbrgtsvm[-c(1,2)]

op_gbrgtsvm <- gbrgtsvm(x, y, lam=gbrgtsvm_op_tuning[1], a=gbrgtsvm_op_tuning[2], 
                        lam1=gbrgtsvm_op_tuning[3], lam2=gbrgtsvm_op_tuning[4])

## plot
plt + 
  geom_abline(slope = -(op_gbrgtsvm$wavg[1]/(op_gbrgtsvm$wavg[2])), 
              intercept = -(op_gbrgtsvm$bavg/(op_gbrgtsvm$wavg[2])), linewidth = 1)+
  geom_abline(slope = -(op_gbtsvm$wavg[1]/(op_gbtsvm$wavg[2])), 
              intercept = -(op_gbtsvm$bavg/(op_gbtsvm$wavg[2])), 
              linewidth = 1, linetype = "dashed")+
  theme(legend.position = c(0.1, 0.8))
ggsave("noise10-comparison.eps", width=5)


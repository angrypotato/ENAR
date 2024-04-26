##################################
## EVA for peer-reviewed paper
## Authors: Zoe Haskell-Craig
## Date Created: April 25th 2024
## Last Modified: April 25th 2024
##################################
#set wd based on who's running the code
setwd("/Users/zoehaskellcraig/Documents/NYU GPH/2023-2024/ENAR 2024/ENAR")

# packages --------
library(readr) #read csv files
library(tidyverse) #data manipulation 
library(POT) #peaks over thresholds methods
library(evmix) #POT mixture model
library(ggpubr) #arranging plots


# data ----------------
data_raw <- read_csv("data/data_weighted.csv")


## remove individuals with missing bp data
data <- data_raw %>% dplyr::filter(!is.na(bp_sys_mean)) #45576 observations removed

## filter by individuals who are not pregnant
data <- data %>% dplyr::filter(demo_pregnant == "No") #15434 observations removed


# split into training/testing sets 70/30 --------------
n <- length(data$id)
train_n <- round(0.7 * n)

#sample row numbers for split
row_train <- sample(c(1:n), train_n)
row_test <- c(1:n)[-row_train]

#split data 
data_train <- data[row_train,]
data_test <- data[row_test,]


# Data transformation: adding noise to break ties ---------------
noise <- rnorm(n = dim(data_train)[1], mean = 0, sd = 0.1)

data_train$bp_sys_mean_noise <- data_train$bp_sys_mean + 
  rnorm(n = dim(data_train)[1], mean = 0, sd = 0.01*data_train$bp_sys_mean)


# EDA ----------------------------
ggplot(data = data_train) +
  geom_histogram(aes(x = bp_sys_mean), bins = 60)

View(table(data_train$bp_sys_mean))
View(table(data_train$bp_sys_mean_noise))

# POT traditional method -------------
## Analysis step 1: determine threshold -----------------

### Method a: threshold range plots ----------
#par(mfrow=c(1,2))
#POT::tcplot(data_train$bp_sys_mean, u.range = c(120, 200))
# ERROR: observed information matrix is singular; passing std.err.type to ``expected''

par(mfrow=c(1,2))
POT::tcplot(data_train$bp_sys_mean_noise[1:400000], u.range = c(120, 240))
#400,000 at sd = 30% works
#400,000 at sd = 10% works
#400,000 at sd = 5% works at 500,000 at sd = 5% works
#500,000 at sd = 1% breaks, 400,000 at sd = 1% breaks less

### Method b: mean excess plot -----------------
par(mfrow=c(1,1))
POT::mrlplot(data_train$bp_sys_mean_noise[1:400000], u.range = c(120, 240)) 
# works at noise = 1%, data = 400,000 observations but quite slow

## Estimate parameters for GPD --------------
gpd_m1 <- fitgpd(data_train$bp_sys_mean_noise, 
                  threshold = 160, "mle")


gpd_m1$param #scale = 15.67 > 0 so a "heavy tail" according to Frigessi, Haug & Rue

# Mixture Modeling approach -------------

## normal + GPD -------------------

#following the recommendations of Hu & , 2018, we carry out a grid search over
#potential threshold listed in useq to find the value which maximizes the likelihood
#then this value is used as the initial value for the threshold in the 
# MLE with complete likelihood (option 5)
threshold_q <- seq(from = 140, to = 200, by = 2)
#3 FIX:: need to double check bulk vs parametric explanation and choose phiu better
normgpd_1 <- fnormgpd(x = data_train$bp_sys_mean_noise[1:100000],
                      useq = threshold_q,
                      fixedu = TRUE,
                      phiu = TRUE) #allowing discontinuity at threshold

# all parameters
round(normgpd_1$mle, 3)
normgpd_1$se

# parameters by name
normgpd_1$u #threshold
normgpd_1$nmean #mean of normal component
normgpd_1$nsd #sd of normal component
normgpd_1$sigmau #GPD scale parameter
normgpd_1$xi #GPD shape parameter

normgpd_1$phiu #tail fraction
normgpd_1$se.phiu #se of tail fraction estimate

normgpd_2 <- fnormgpdcon(x = data_train$bp_sys_mean_noise[1:100000],
                         useq = threshold_q,
                         fixedu = TRUE) #no discontinuity at threshold 

# parameters by name
normgpd_2$u #threshold
normgpd_2$nmean #mean of normal component
normgpd_2$nsd #sd of normal component
normgpd_2$sigmau #GPD scale parameter
normgpd_2$xi #GPD shape parameter

normgpd_2$phiu #tail fraction
normgpd_2$se.phiu #se of tail fraction estimate

### model diagnostics -------------------
rlplot(normgpd_1)
rlplot(normgpd_2)


##################################
## EVA
## Authors: Zoe Haskell-Craig
## Date Created: Dec 13th 2023
## Last Modified: Feb 26th 2024
##################################

# packages --------
library(readr) #read csv files
library(tidyverse) #data manipulation 
library(POT) #peaks over thresholds methods
library(ggpubr) #arranging plots

# data -----------
setwd("/Users/zoehaskellcraig/Documents/NYU GPH/2023-2024/ENAR 2024/ENAR")
data_raw <- read_csv("data/data_imputed.csv")


## remove individuals with missing bp data
data <- data_raw %>% dplyr::filter(!is.na(bp_sys_mean))

## filter by individuals diagnosed with hypertension
data_htn <- filter(data, htn_jnc7 == 'Yes')


# Step 1: EDA ----------------

## See EDA.Rmd for more exploratory data analysis



# Analysis 1 - Expected extreme BP given data 1999-2012 -----------------


#split data by period (1999-2012)
data_1 <- data_htn %>% filter(svy_year %in% c("1999-2000","2001-2002","2003-2004",
                                          "2005-2006","2007-2008","2009-2010",
                                          "2011-2012"))
#split data by period (2013-2020)
data_2 <- data_htn %>% filter(svy_year %in% c("2013-2014", "2015-2016","2017-2020"))

## Method a: threshold range plots ----------
par(mfrow=c(1,2))
POT::tcplot(data_1$bp_sys_mean, u.range = c(120, 200)) #from this looks like mu = 140- 160



## Method b: mean excess plot/mean residual life plot -------
par(mfrow=c(1,1))
#first four waves
POT::mrlplot(data_1$bp_sys_mean, u.range = c(120, 220)) 
#abline(v = 158, lty = 2)
abline(v = 140, lty = 2)
abline(v = 172, lty = 2)
#abline(v = 170)

## Method d: 90% quantile ---------
print(quantile(data_1$bp_sys_mean, 
               probs = c(0.5, 0.75, 0.83, 0.85, 0.9, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99)))




## estimate parameters for GPD --------------
gpd_m1a <- fitgpd(data_1$bp_sys_mean, 
                 threshold = 140, "mle")

gpd_m1b <- fitgpd(data_1$bp_sys_mean, 
                  threshold = 160, "mle")


#model fit statistics
gpd_m1a$logLik
gpd_m1b$logLik

#AIC m1a = 49182.04
#AIC m1b = 15423.87

# looks like b is the better model

#parameter estimates
gpd_m1b$param

#parameter confidence intervals
gpd.fiscale(gpd_m1b)
gpd.fishape(gpd_m1b)




## model diagnostics -----------
n_per_survey_1 <- data_1 %>% group_by(svy_year) %>% count()
npy_1 = mean(n_per_survey_1$n)

par(mfrow=c(2,2))
plot(gpd_m1b, npy = npy_1)



## calculate return period ----------
par(mfrow=c(1,1))
retlev(gpd_m1b, npy = npy_1, points = T)

return_level_1 <- retlev(gpd_m1b, npy = npy_1, points = T)
return_level_1(npy_1)

gpd_m1b$pat #proportion of exceedances

#probability that x > 180
1 - pgpd(180, loc = 160, scale = 18.0808248, shape = -0.1136259)
#probability that x > 200
1 - pgpd(200, loc = 160, scale = 18.0808248, shape = -0.1136259)

# the 99th percentile value for BP
qgpd(p = 0.99, loc = 160, scale = 18.0808248, shape = -0.1136259)


## observed vs expected extreme values 2013-2020 ---------------

# number of individuals (w/o missing bp data) surveyed 2013-2020
n_samp <- length(data_2$id)
n_over_160 <- length(filter(data_2, bp_sys_mean > 160)$id)

# expected sys > 160
e_sys_160 <- gpd_m1b$pat * n_samp
#observed sys > 160
o_sys_160 <- n_over_160

# contigency table for sys > 160
t_160 <- c(length(filter(data_1, bp_sys_mean < 160)$id), 
           length(filter(data_1, bp_sys_mean > 160)$id),
           length(data_2$id),
           n_over_160)

chisq.test(t_160)

# among observations of extreme blood pressure
# expected sys > 180 
e_sys_180 <- 0.3066362 * n_over_160
o_sys_180 <- length(filter(data_2, bp_sys_mean > 180)$id)

# expected sys > 200 
e_sys_200 <- 0.07824075 * n_over_160
o_sys_200 <- length(filter(data_2, bp_sys_mean > 200)$id)

# expected vs observed
chisq.test(x = c(e_sys_180, e_sys_200), y = c(o_sys_180, o_sys_200))


# Analysis 2 - extreme BP threshold quantile for QR ------------------
# what values constitute extreme tail of distribution, given individual
# is diagnosed with high blood pressure?


data_h <- filter(data, htn_jnc7 == "Yes", #only individuals with hypertension
                 demo_pregnant == "No") ## remove individuals who are pregnant


## Method a: threshold range plots ----------
# first four waves
POT::tcplot(data_h$bp_sys_mean, u.range = c(120, 200))  




## Method b: mean excess plot/mean residual life plot -------
#first four waves
POT::mrlplot(data_h$bp_sys_mean, u.range = c(120, 220)) #mu = 140, 160, 170
abline(v = 160)
abline(v = 140)
abline(v = 170)

## Method c: L-moments plot -----------
#first four waves
POT::lmomplot(data_h$bp_sys_mean, u.range = c(120, 220))

## Method d: 90% quantile ---------
print(quantile(data_h$bp_sys_mean, 
               probs = c(0.5, 0.75, 0.85, 0.9,0.91, 0.92, 0.95,
                         0.96, 0.97, 0.98, 0.99)))




## Fit model ------------
gpd_m2a <- fitgpd(data_h$bp_sys_mean, 
                  threshold = 160, "mle")
gpd_m2b <- fitgpd(data_h$bp_sys_mean, 
                  threshold = 170, "mle") #lower AIC


##  Check model fit ------------

## QQplot, density plot, and return level plot
# need to calculate estimated number of events (participants) per year (survey wave)
svy_num <- data_h %>% group_by(svy_year) %>% count()
avg_svy_num <- mean(svy_num$n)

par(mfrow=c(2,2))
plot(gpd_m2a, npy = avg_svy_num)


# Fig 1 plot --------------

data_1$hypertension <- as.factor(data_1$htn_jnc7)

g1 <- ggplot(data_1) +
  geom_histogram(aes(x = bp_sys_mean, 
                     fill = bp_sys_mean < 160), bins = 38) +
  scale_fill_manual(values = c( "#F4D06F", "#644375"),
                    name = "Blood Pressure", labels = c(">160 mm Hg", "<160 mm Hg")) +
  xlim(75, 225) +
  ylim(0, 4000) +
  labs(title = "Survey years from 1999 - 2012",
       x = "Systolic BP (mm Hg)", y = "Num. observations") +
  #geom_vline(xintercept = 130, color = "grey", lty = 2) + 
  #geom_vline(xintercept = 140, color = "grey", lty = 2) +
  theme(axis.text = element_text(size = 12),  
        axis.title = element_text(size = 14)) #+


g2 <- ggplot(data_2) +
  geom_histogram(aes(x = bp_sys_mean, 
                     fill = bp_sys_mean < 160), bins = 38) +
  scale_fill_manual(values = c( "#F4D06F", "#644375"),
                    name = "Blood Pressure", labels = c(">160 mm Hg", "<160 mm Hg")) +
  xlim(75, 225) +
  ylim(0, 4000) + 
  labs(title = "Survey years from 2013 - 2020",
       x = "Systolic BP (mm Hg)", y = "Num. observations") +
  #geom_vline(xintercept = 130, color = "grey", lty = 2) + 
  #geom_vline(xintercept = 140, color = "grey", lty = 2) +
  theme(axis.text = element_text(size = 12),  
        axis.title = element_text(size = 14)) #+


data_np <- filter(data, demo_pregnant == "No") ## remove individuals who are pregnant


g3 <- ggplot(data_np) +
  geom_histogram(aes(x = bp_sys_mean, 
                     fill = as.factor(htn_jnc7)), bins = 30) +
  scale_fill_manual(values = c( "#ABA9BF", "#11ADD0"),
                    name = "Hypertension", labels = c("No", "Yes")) +
  xlim(75, 225) +
  #ylim(0, 4000) + 
  labs(title = "Individuals with hypertension, 1999 - 2020",
       x = "Systolic BP (mm Hg)", y = "Num. observations") +
  geom_vline(xintercept = 160, color = "black", lty = 2) + 
  geom_vline(xintercept = 170, color = "black", lty = 2) +
  theme(axis.text = element_text(size = 12),  
        axis.title = element_text(size = 14)) #+
#guides(fill = "none")  # Remove legend

ggarrange(g1, g2, ncol = 2,
          common.legend = TRUE, legend = "right")

## participants with hypertension surveyed 1999-2013
bin_n <- round(max(data_1$bp_sys_mean) - min(data_1$bp_sys_mean))
ggplot(data_1) +
  geom_histogram(aes(x = bp_sys_mean, 
                     fill = bp_sys_mean > 160), 
                 breaks = c(70, 75, 80, 85,  90, 95,  100, 105, 110, 115,
                            120, 125, 130, 135, 140, 145, 150, 155, 160, 165,
                            170, 175, 180, 185, 190, 195, 200,
                            205, 210, 215, 220, 225, 230)) +
  scale_fill_manual(values = c( "#F4D06F", "#644375"),
                    name = "Blood Pressure", 
                    labels = c("<160 mm Hg", ">160 mm Hg")) +
  xlim(76, 225) +
  ylim(0, 2000) + 
  labs(title = "Individuals with hypertension", 
       subtitle = "Survey years from 1999 - 2012",
       x = "Systolic BP (mm Hg)", y = "Num. observations") +
  #geom_vline(xintercept = 130, color = "grey", lty = 2) + 
  #geom_vline(xintercept = 140, color = "grey", lty = 2) +
  theme(axis.text = element_text(size = 12),  
        axis.title = element_text(size = 14)) #+
  


##################################
## EVA
## Authors: Zoe Haskell-Craig
## Date Created: Dec 13th 2023
## Last Modified: Dec 15th 2023
##################################

# packages --------
library(tidyverse) #data manipulation 
library(POT) #peaks over thresholds methods

# data -----------
load("ENAR/data/nhanes_data.rda")

#split by survey year (train on 4 waves, predict prevalence on next 4 waves)
data_train1 <- nhanes_data %>%
  filter(svy_year %in% c("1999-2000", "2001-2002", "2003-2004", "2005-2006"))


## data cleaning ------------
# filter by observations without missing blood pressure data (sys_mean)
data_train1_clean <- data_train1 %>%
  select(svy_id, svy_year, bp_sys_mean) %>%
  filter(!is.na(bp_sys_mean)) #21205 - 20095 = 1110 observations removed

# create dataset for yearly means and sd
## mean blood pressure
bp_yearly_data <- nhanes_data %>%
  group_by(svy_year, demo_race) %>%
  transmute(survey_year = svy_year,
            race = demo_race,
            BP_sys = mean(bp_sys_mean, na.rm = TRUE),
            BP_sys_sd = sd(bp_sys_mean, na.rm = TRUE),
            BP_dia = mean(bp_dia_mean, na.rm = TRUE),
            BP_dia_sd = sd(bp_dia_mean, na.rm = TRUE),
            n = n()) %>%
  unique() %>%
  mutate(lower.sys = BP_sys - 1.96*BP_sys_sd/sqrt(n),
         upper.sys = BP_sys + 1.96*BP_sys_sd/sqrt(n),
         lower.dia = BP_dia - 1.96*BP_dia_sd/sqrt(n),
         upper.dia = BP_dia + 1.96*BP_dia_sd/sqrt(n))

# Step 1: EDA ----------------

# yearly mean BP and 95% CI
ggplot(data = bp_yearly_data) +
  geom_point(aes(x = survey_year, y = BP_sys, color = race)) +
  geom_line(aes(x = survey_year, y = BP_sys, color = race, group = race)) +
  geom_linerange(aes(x = survey_year, y = BP_sys, color = race,
                    ymin = lower.sys, ymax = upper.sys)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## distribution of outcome/outliers for each year of data, by race
ggplot(data = nhanes_data) +
  geom_boxplot(aes(x = svy_year, y = bp_sys_mean,  
                   fill = demo_race)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## distribution of outcome - for all observations 1999-2006
ggplot(data = data_train1_clean) + geom_histogram(aes(x = bp_sys_mean)) +
  labs(title = "Distribution of Mean Systolic BP",
       x = "BP", y = "Num. observations") +
  geom_vline(xintercept = 130, color = "blue") + #
  geom_vline(xintercept = 140, color = "blue")


# Step 2: Find the threshold value -------------

## Method a: threshold range plots ----------
# first four waves
POT::tcplot(data_train1_clean$bp_sys_mean, u.range = c(140, 200)) #looks like mu = 171?


## Method b: mean excess plot/mean residual life plot -------
#first four waves
POT::mrlplot(data_train1_clean$bp_sys_mean, u.range = c(140, 220)) #also looks like mu = 172?
abline(v = 172)

## Method c: L-moments plot -----------
#first four waves
POT::lmomplot(data_train1_clean$bp_sys_mean, u.range = c(140, 220))

## Method d: 90% quantile ---------
print(quantile(data_train1_clean$bp_sys_mean, 
               probs = c(0.5, 0.75, 0.85, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99)))


# Step 3: Estimate the parameters for Grand Pareto distribution ------
gpd_m0 <- fitgpd(data_train1_clean$bp_sys_mean, 
                 threshold = 172, "mle")

#threshold 
gpd_m0$threshold

#parameter estimates
gpd_m0$param

#parameter confidence intervals
gpd.fiscale(gpd_m0)
gpd.fishape(gpd_m0)


# Step 4: Check model fit ------------

## QQplot, density plot, and return level plot
# need to calculate estimated number of events (participants) per year (survey wave)
svy_num <- data_train1_clean %>% group_by(svy_year) %>% count()
avg_svy_num <- mean(svy_num$n)

par(mfrow=c(2,2))
plot(gpd_m0, npy = avg_svy_num)



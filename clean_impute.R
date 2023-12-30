####################################
## Data Cleaning and Missing Data Imputation
## Date Created: Dec 29th 2023
## Last Modified: Dec 29th 2023
## Author: Zoe Haskell-Craig and Jianan Zhu
####################################

# packages -----------
library(tidyverse)


# data ----------------

# data dictionary
load("data/nhanes_key.rda")
load("data/key_guide.rda")

dataset.complete.raw <- read.csv("data/nhanes_complete_raw.csv")


# select covariates

# 1. `bp_sys_mean`: Systolic blood pressure (SBP), mm Hg
# 2. `bp_dia_mean`: Diastolic blood pressure (DBP), mm Hg
# 3. `demo_age_years`: Age, years
# 4. `demo_gender`: Gender 
# 5. `cc_bmi`: Body mass index, kg/m2
# 6. `chol_total`: Total cholesterol, mg/dL
# 7. `chol_hdl`: HDL cholesterol, mg/dL
# 8. `chol_ldl`: LDL cholesterol, mg/dL
# 9. `INDFMPIR`: Family PIR -- continuous, 0 - 4.99, 5 = PIR value greater than or equal to 5.00
# 10. `DMDEDUC2`: Education level - Adults 20+ -- 7 = refused, 9 = don't know
# 11. `DMDMARTL`: Marital status -- 77 = refused, 99 = don't know
# 12. `ALQ130`: Avg # alcoholic drinks/day - past 12 mos -- 15 = 15 drinks or more, 
# 777 = refused, 999 = don't know
# 13. `SMQ020`: Smoked at least 100 cigarettes in life --- 7 = refused, 9 = don't know
# 14. `demo_pregnant`: Pregnant - note we will exclude individuals who are pregnant
# 15. `htn_jnc7`: Hypertension defined by the JNC7 guideline



dataset <- dataset.complete.raw[,c("svy_year", "demo_pregnant" , "htn_jnc7" ,"bp_sys_mean","bp_dia_mean","demo_gender","demo_age_years","cc_bmi","chol_total","chol_hdl","chol_ldl","INDFMPIR","DMDEDUC2","DMDMARTL","ALQ130","SMQ020")]


# data cleaning ------------

## replace 'missing', 'don't know' or 'refused' with NA ------
dataset_clean <- dataset %>% 
  mutate(DMDEDUC2 = if_else(DMDEDUC2 == 7 | DMDEDUC2 == 9,
                            NA_real_, DMDEDUC2),
         DMDMARTL = if_else(DMDMARTL == 77 | DMDMARTL == 99, 
                              NA_real_, DMDMARTL),
         ALQ130 = if_else(ALQ130 > 15, 
                            NA_real_, ALQ130),
         SMQ020 = if_else(SMQ020 == 7 | SMQ020 == 9,
                            NA_real_, SMQ020),
         chol_ldl = if_else(chol_ldl < 0, 
                            NA_real_, chol_ldl))





## change categorical variables to factors -----------

dataset_clean$svy_year<- factor(dataset_clean$svy_year, levels = unique(dataset_clean$svy_year))
dataset_clean$demo_pregnant<- factor(dataset_clean$demo_pregnant, 
                                     levels = unique(dataset_clean$demo_pregnant))
dataset_clean$htn_jnc7<- factor(dataset_clean$htn_jnc7, levels = unique(dataset_clean$htn_jnc7))
dataset_clean$demo_gender<- factor(dataset_clean$demo_gender, levels = unique(dataset_clean$demo_gender))
dataset_clean$cc_bmi<- factor(dataset_clean$cc_bmi, levels = unique(dataset_clean$cc_bmi))
dataset_clean$DMDEDUC2<- factor(dataset_clean$DMDEDUC2, levels = unique(dataset_clean$DMDEDUC2))
dataset_clean$DMDMARTL<- factor(dataset_clean$DMDMARTL, levels = unique(dataset_clean$DMDMARTL))
dataset_clean$SMQ020<- factor(dataset_clean$SMQ020, levels = unique(dataset_clean$SMQ020))
summary(dataset_clean)

readr::write_csv(dataset_clean, file = "data/dataset_clean.csv")

# Missing value imputation -------------

## train/test split ----------------
## check missing data imputation performance for individuals with higher and lower bp
test_id_HBP <- sample(filter(dataset_clean))





imp <- mice(dataset, m=2, maxit=2, 
            method = c("rf", "pmm", "pmm", "rf","pmm","rf","pmm",
                       "pmm","pmm","pmm","rf","rf","pmm","rf")) #rf for categorical, pmm for continuous


dataset_imp <- complete(imp)
summary(dataset_imp)




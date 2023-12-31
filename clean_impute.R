####################################
## Data Cleaning and Missing Data Imputation
## Date Created: Dec 29th 2023
## Last Modified: Dec 31st 2023
## Author: Zoe Haskell-Craig and Jianan Zhu
####################################

# packages -----------
library(tidyverse) #data cleaning/manipulation
library(readr) #read/write csv files
library(mice) #missing data imputation


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
# Dropping NAs in cc_bmi results in 55,121 observations
# Dropping NAs in chol_total results in 55,121 
# Dropping NAs in INDFMPIR results in 51,136
# Dropping NAs in DMDEDUC2 results in 52,384
# Dropping NAs in DMDMARTL results in 54,010
# Dropping NAs in ALQ130 results in 33,216
# Dropping NAs in SMQ020 results in 53,336 
# Dropping NAs in chol_hdl results in 24,043
# Dropping NAs in chol_ldl results in 24,020 

# Dropping NAs in all above results in 12,892 
# Dropping NAs in cc_bmi,chol_total,INDFMPIR,DMDEDUC2,DMDMARTL,ALQ130,and SMQ020 results in 12,903 

# So ALQ130, SMQ020, chol_hdl, chol_ldl are covariates for whom imputation is most likely
# to change the distribution of values in the sample

## train/test split ----------------
## create IDs/rownumbers
dataset_clean <- dataset_clean %>% mutate(id=row_number())

## check missing data imputation performance for individuals with higher and lower bp
test_id_HBP <- sample(filter(dataset_clean, bp_sys_mean > 140)$id, 700)
test_id_LBP <- sample(filter(dataset_clean, bp_sys_mean < 140)$id, 700)

# select 40% of data randomly from rows in test_id_HBP and test_id_LBP from 
# ALQ130, SMQ020, chol_hdl, chol_ldl for removal
ALQ130_rm <- sample(c(test_id_HBP, test_id_LBP), 560)
SMQ020_rm <- sample(c(test_id_HBP, test_id_LBP), 560)
chol_hdl_rm <- sample(c(test_id_HBP, test_id_LBP), 560)
chol_ldl_rm <- sample(c(test_id_HBP, test_id_LBP), 560)

# replace values with NA
dataset_train <- dataset_clean %>%
  mutate(ALQ130 = if_else(id %in% ALQ130_rm, NA_real_, ALQ130),
         SMQ020 = if_else(id %in% SMQ020_rm, NA, SMQ020),
         chol_hdl = if_else(id %in% chol_hdl_rm, NA_real_, chol_hdl),
         chol_ldl = if_else(id %in% chol_ldl_rm, NA_real_, chol_ldl))



## impute missing data ------------------

# remove col with row numbers/ids, blood pressure
ids <- dataset_train$id #save ids
bp_sys_mean <- dataset_train$bp_sys_mean #save bp
bp_dia_mean <- dataset_train$bp_dia_mean #save bp

dataset_train <- dataset_train %>% select(-id, -bp_dia_mean, -bp_sys_mean)

imp <- mice(dataset_train, m=5, maxit=5, 
            method = c("rf", "rf", "rf", "rf", #svy year - gender
                       "pmm", "rf", #age, bmi
                       "pmm","pmm", "pmm", #chol
                       "pmm", "rf", "rf", #demographics
                       "pmm", "rf" #alq, smoking
                       )) #rf for categorical, pmm for continuous


dataset_imp <- complete(imp)
summary(dataset_imp)
dataset_imp$id <- ids

## check imputed data performance WRT holdout observations ----------

### ALQ130 MSE
true_ALQ130 <- filter(dataset_clean, id %in% ALQ130_rm)$ALQ130
imp_ALQ130 <- filter(dataset_imp, id %in% ALQ130_rm)$ALQ130
sqerr_ALQ130 <- (true_ALQ130 - imp_ALQ130)^2

err_ALQ_LBP <- sqerr_ALQ130[which((filter(dataset_clean, id %in% ALQ130_rm)$id) %in% test_id_LBP)]
err_ALQ_HBP <- sqerr_ALQ130[which((filter(dataset_clean, id %in% ALQ130_rm)$id) %in% test_id_HBP)]
print(mean(err_ALQ_LBP, na.rm = TRUE))
print(mean(err_ALQ_HBP, na.rm = TRUE))

t.test(err_ALQ_LBP, err_ALQ_HBP, na.rm = TRUE) #not statistically significant difference in mean MSE


### SMQ020 two by two table
true_SMQ020 <- filter(dataset_clean, id %in% SMQ020_rm)$SMQ020
imp_SMQ020 <- filter(dataset_imp, id %in% SMQ020_rm)$SMQ020

table(true_SMQ020, imp_SMQ020)

true_SMQ020_LBP <- true_SMQ020[which((filter(dataset_clean, id %in% SMQ020_rm)$id) %in% test_id_LBP)]
imp_SMQ020_LBP <- imp_SMQ020[which((filter(dataset_clean, id %in% SMQ020_rm)$id) %in% test_id_LBP)]

table(true_SMQ020_LBP, imp_SMQ020_LBP)

true_SMQ020_HBP <- true_SMQ020[which((filter(dataset_clean, id %in% SMQ020_rm)$id) %in% test_id_HBP)]
imp_SMQ020_HBP <- imp_SMQ020[which((filter(dataset_clean, id %in% SMQ020_rm)$id) %in% test_id_HBP)]

table(true_SMQ020_HBP, imp_SMQ020_HBP)

### chol_hdl MSE
true_chol_hdl <- filter(dataset_clean, id %in% chol_hdl_rm)$chol_hdl
imp_chol_hdl <- filter(dataset_imp, id %in% chol_hdl_rm)$chol_hdl
sqerr_chol_hdl <- (true_chol_hdl - imp_chol_hdl)^2

err_hdl_LBP <- sqerr_chol_hdl[which((filter(dataset_clean, id %in% chol_hdl_rm)$id) %in% test_id_LBP)]
err_hdl_HBP <- sqerr_chol_hdl[which((filter(dataset_clean, id %in% chol_hdl_rm)$id) %in% test_id_HBP)]
print(mean(err_hdl_LBP, na.rm = TRUE))
print(mean(err_hdl_HBP, na.rm = TRUE))

t.test(err_hdl_LBP, err_hdl_HBP, na.rm = TRUE) #not statistically sig. different

### chol_ldl MSE
true_chol_ldl <- filter(dataset_clean, id %in% chol_ldl_rm)$chol_ldl
imp_chol_ldl <- filter(dataset_imp, id %in% chol_ldl_rm)$chol_ldl
sqerr_chol_ldl <- (true_chol_ldl - imp_chol_ldl)^2

err_ldl_LBP <- sqerr_chol_ldl[which((filter(dataset_clean, id %in% chol_ldl_rm)$id) %in% test_id_LBP)]
err_ldl_HBP <- sqerr_chol_ldl[which((filter(dataset_clean, id %in% chol_ldl_rm)$id) %in% test_id_HBP)]
print(mean(err_ldl_LBP, na.rm = TRUE))
print(mean(err_ldl_HBP, na.rm = TRUE))

t.test(err_ldl_LBP, err_ldl_HBP, na.rm = TRUE) #not stat sig difference


## return true values of hold out observations ---------

dataset_final <- dataset_imp 
# add original columns from holdout variables
dataset_final$ALQ_orig <- dataset_clean$ALQ130 #original col
dataset_final$SMQ020_orig <- dataset_clean$SMQ020
dataset_final$chol_hdl_orig <- dataset_clean$chol_hdl
dataset_final$chol_ldl_orig <- dataset_clean$chol_ldl
# re-add IDs, blood pressure
dataset_final$id <- dataset_clean$id
dataset_final$bp_sys_mean <- dataset_clean$bp_sys_mean
dataset_final$bp_dia_mean <- dataset_clean$bp_dia_mean

# if orginal cols for holdout variables did not have NA, update col with true val
dataset_final <- dataset_final %>% 
  mutate(ALQ130 = if_else(id %in% ALQ130_rm & !is.na(ALQ_orig), ALQ_orig, ALQ130),
         SMQ020 = if_else(id %in% SMQ020_rm & !is.na(SMQ020_orig), SMQ020_orig, SMQ020),
         chol_ldl = if_else(id %in% chol_ldl_rm & !is.na(chol_ldl_orig), chol_ldl_orig, chol_ldl),
         chol_hdl = if_else(id %in% chol_hdl_rm & !is.na(chol_hdl_orig), chol_hdl_orig, chol_hdl)) %>%
  select(-ALQ_orig, -SMQ020_orig, -chol_hdl_orig, -chol_ldl_orig) #remove un-imputed col

summary(dataset_final)

## save imputed dataset ---------
write_csv(dataset_final, file = "data/data_imputed.csv")


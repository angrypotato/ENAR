This repository contains the code for the ENAR datafest submission from Zoe Haskell-Craig, Jianan Zhu, Iris Zhang and Xiaoting Chen.

All data is housed in the `data` folder.

|`extract_add_covariates.Rmd` contains code to created `nhanes_complete_raw.csv` 
|
|--> `clean_impute.R` contains code to clean data and impute missing values. Output for cleaned data is `dataset_clean.csv`. Output for |imputed data is `data_imputed.csv`. All subsequent files use `data_imputed.csv`
|
|
|---> `EDA.Rmd` is the exploratory data analysis code
|
|
|----> `EVA.R` contains code for the extreme values analysis
|
|
|-----> `QR.Rmd` contains code for the quantile regression analysis
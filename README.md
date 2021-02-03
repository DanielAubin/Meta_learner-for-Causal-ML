# Meta-learners for Causal-ML

This repository collects Meta-learners for causal inference estimation. 

The datastructure has to be the following: 

* y = outcome variable
* d = treatment variable
* covariates = the covariates to map on d or y
* learners = the machine learning methods to use in the [SuperLearner package](https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html).
* df_aux and df_main = the train and test-set. df_aux is used for training the nuisance functions while df_main is used to estimate the CATE. 

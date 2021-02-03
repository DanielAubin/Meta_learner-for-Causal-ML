# Meta_learner-for-Causal-ML

This repository collects Meta-Learners for causal inference estimation. 
The datastructure has to be the following:

y = outcome variable
d = treatment variable
covariates = the covariates to map on d or y
learners = the machine learning methods to use in the SuperLearner package (see: https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html for an introduction to the SuperLearner)

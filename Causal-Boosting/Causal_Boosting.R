remotes::install_github("saberpowers/causalLearning")
library(causalLearning)

Causal_Boost <- function(df_aux,df_main){
fit_cb = causalBoosting(y=df_aux$y,tx=df_aux$d,x=df_aux[,covariates], num.trees = 500)

score_CBoost = predict(fit_cb, newx =df_main[,covariates] , num.trees = 500)

return(score_CBoost)
}

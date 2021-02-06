install.packages("xgboost", repos=c("http://dmlc.ml/drat/", getOption("repos")), type="source")
vec.pac= c("SuperLearner", "gbm", "glmnet","ranger")

lapply(vec.pac, require, character.only = TRUE) 




#Learner Library:
learners <- c( "SL.glmnet","SL.xgboost", "SL.ranger","SL.lm","SL.mean")

#CV Control for the SuperLearner
control <- SuperLearner.CV.control(V=5)



R_learner <- function(df_aux,df_main,covariates,learners){

# Train a classification model to get the propensity scores
p_mod <- SuperLearner(Y = df_aux$d, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                      verbose = FALSE, method = "method.NNLS", family = binomial(),cvControl = control)

p_hat <- p_mod$SL.predict
p_hat = ifelse(p_hat<0.025, 0.025, ifelse(p_hat>.975,.975, p_hat)) # Overlap bounding

# Train a regression model 
m_mod <- SuperLearner(Y = df_aux$y, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                      verbose = FALSE, method = "method.NNLS",cvControl = control)

m_hat <- m_mod$SL.predict

# Apply the R-learner (residual-on-residual approach)
y_tilde = df_main$y - m_hat
w_tilde = df_main$d - p_hat
pseudo_outcome = y_tilde/w_tilde

weights = w_tilde^2


a  <- tryCatch({
  R_mod <- SuperLearner(Y = pseudo_outcome, X = df_main[,covariates], newX = df_main[,covariates], SL.library = learners,
                        verbose = FALSE, method = "method.NNLS",obsWeights = weights[,1],cvControl = control)
  score_R <- R_mod$SL.predict
  a <- score_R
  
  
},error=function(e){
  
  mean_score <- weighted.mean(pseudo_outcome, w = weights)
  score_R <- rep.int(mean_score, times = nrow(df_main))
  a <- score_R
  return(a)
})

score_R <- a

return(score_R)
}

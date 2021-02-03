T_learner <- function(df_aux,df_main,covariates){
  
  
  aux_1 <- df_aux[which(df_aux$d==1),]
  aux_0 <- df_aux[which(df_aux$d==0),]
  
  m1_mod <- SuperLearner(Y = aux_1$y, X = aux_1[,covariates], newX = df_main[,covariates], SL.library = learners,
                         verbose = FALSE, method = "method.NNLS", family = binomial())
  
  m1_hat <- m1_mod$SL.predict
  
  m0_mod <- SuperLearner(Y = aux_0$y, X = aux_0[,covariates], newX = df_main[,covariates], SL.library = learners,
                         verbose = FALSE, method = "method.NNLS", family = binomial())
  
  m0_hat <- m0_mod$SL.predict
  
  score_t <- m1_hat - m0_hat
  
  return(score_t)
}

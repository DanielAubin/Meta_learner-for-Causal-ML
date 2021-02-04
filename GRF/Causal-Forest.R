
Causal_Forest <- function(df_aux,df_main,covariates,learners){
  

  p_mod <- SuperLearner(Y = df_aux$d, X = df_aux[,covariates], newX = df_aux[,covariates], SL.library = learners,
                        verbose = FALSE, method = "method.NNLS", family = binomial(),cvControl = control)
  
  p_hat <- p_mod$SL.predict
  
  p_hat <- ifelse(p_hat<0.025, 0.025, ifelse(p_hat>.975,.975, p_hat)) # Overlap bounding  
  

  m_mod <- SuperLearner(Y = df_aux$y, X = df_aux[,covariates], newX = df_aux[,covariates], SL.library = learners,
                         verbose = FALSE, method = "method.NNLS",cvControl = control)

  m_hat <- m_mod$SL.predict


tau.forest <- causal_forest(df_aux[,covariates], df_aux$y, df_aux$d,
                           W.hat = p_hat, Y.hat = m_hat,
                           tune.parameters = "all")


tau.hat <- predict(tau.forest, df_main[,covariates])

return(tau.hat$predictions)

}

test <- Causal_Forest(df_aux,df_main,covariates,learners)

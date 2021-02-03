TOM <- function(df_aux, df_main,covariates,learners){
  
  # Propensity score
  p <- rep((nrow(df_aux[df_aux[,"d"]==1,]) + nrow(df_main[df_main[,"d"]==1,]))/(nrow(df_main) + nrow(df_main)), nrow(df_aux)) 
  
  # Transformed outcome
  Y_TO <- df_aux$y * df_aux$d/p - df_aux$y * (1-df_aux$d)/(1-p)
  
  m_TO <- SuperLearner(Y =Y_TO, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                           verbose = FALSE, method = "method.NNLS")
  
 
  
  tau_TO <- m_TO$SL.predict
 
  
  return(tau_TO)
  }

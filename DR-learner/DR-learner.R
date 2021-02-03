DR_learner <- function(df_aux,df_main,covariates,learners){


  p_mod <- SuperLearner(Y = df_aux$d, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                        verbose = FALSE, method = "method.NNLS", family = binomial())
  
  p_hat <- p_mod$SL.predict
  p_hat = ifelse(p_hat<0.025, 0.025, ifelse(p_hat>.975,.975, p_hat)) # Overlap bounding
  
  aux_1 <- df_aux[which(df_aux$d==1),]
  aux_0 <- df_aux[which(df_aux$d==0),]
  
  m1_mod <- SuperLearner(Y = aux_1$y, X = aux_1[,covariates], newX = df_main[,covariates], SL.library = learners,
                         verbose = FALSE, method = "method.NNLS", family = binomial())
  
  m1_hat <- m1_mod$SL.predict
  
  m0_mod <- SuperLearner(Y = aux_0$y, X = aux_0[,covariates], newX = df_main[,covariates], SL.library = learners,
                         verbose = FALSE, method = "method.NNLS", family = binomial())
  
  m0_hat <- m0_mod$SL.predict




y_mo <- (m1_hat - m0_hat) + ((df_main$d*(df_main$y -m1_hat))/p_hat) - ((1-df_main$d)*(df_main$y - m0_hat)/(1-p_hat))





a  <- tryCatch({
  dr_mod <- SuperLearner(Y = y_mo, X = df_main[,covariates], newX = df_main[,covariates], SL.library = learners,
                         verbose = FALSE, method = "method.NNLS2")
  
  score_dr <- dr_mod$SL.predict
  a <- score_dr
  
  
},error=function(e){
  
  mean_score <- mean(y_mo)
  score_dr <- rep.int(mean_score, times = nrow(test_data))
  a <- score_dr
  return(a)
})

score_dr <- a

return(score_dr)
}


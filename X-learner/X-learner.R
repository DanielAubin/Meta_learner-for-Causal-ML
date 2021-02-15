install.packages("xgboost", repos=c("http://dmlc.ml/drat/", getOption("repos")), type="source")
vec.pac= c("SuperLearner", "gbm", "glmnet","ranger")

lapply(vec.pac, require, character.only = TRUE) 




#Learner Library:
learners <- c( "SL.glmnet","SL.xgboost", "SL.ranger","SL.lm","SL.mean")

#CV Control for the SuperLearner
control <- SuperLearner.CV.control(V=5)


X_learner <- function(df_aux,df_main,covariates,learners){

  
  p_mod <- SuperLearner(Y = df_aux$d, X = df_aux[,covariates], newX = df_aux[,covariates], SL.library = learners,
                        verbose = FALSE, method = "method.NNLS", family = binomial(),cvControl = control)
  
  # Prop-Score for df_main
  p_hat_main <- predict(p_mod,df_main[,covariates])$pred
  p_hat_main = ifelse(p_hat<0.025, 0.025, ifelse(p_hat>.975,.975, p_hat)) # Overlap bounding
  
  aux_1 <- df_aux[which(df_aux$d==1),]
  aux_0 <- df_aux[which(df_aux$d==0),]
  
  m1_mod <- SuperLearner(Y = aux_1$y, X = aux_1[,covariates], newX = df_aux[,covariates], SL.library = learners,
                         verbose = FALSE, method = "method.NNLS",cvControl = control)
  
  m1_hat <- m1_mod$SL.predict
  
  m0_mod <- SuperLearner(Y = aux_0$y, X = aux_0[,covariates], newX = df_aux[,covariates], SL.library = learners,
                         verbose = FALSE, method = "method.NNLS",cvControl = control)
  
  m0_hat <- m0_mod$SL.predict
  
  
  
  
  tau1 <- df_aux[which(df_aux$d==1),"y"] - m0_hat[which(df_aux$d==1),]
  
  tau0 <- m1_hat[which(df_aux$d==0),]  -  df_aux[which(df_aux$d==0),"y"]
  
  
  a1  <- tryCatch({
    tau1_mod <- SuperLearner(Y = tau1, X = df_aux[which(df_aux$d==1),covariates], newX = df_main[,covariates], SL.library = learners,
                             verbose = FALSE, method = "method.NNLS",cvControl = control)
    
    score_tau1 <- tau1_mod$SL.predict
    a1 <- score_tau1
    
    
  },error=function(e){
    
    mean_score <- mean(tau1)
    score_tau1 <- rep.int(mean_score, times = nrow(df_aux))
    a1 <- score_tau1
    return(a1)
  })
  
  score_tau1 <- a1
  
  a0  <- tryCatch({
    tau0_mod <- SuperLearner(Y =tau0, X = df_aux[which(df_aux$d==0),covariates], newX = df_main[,covariates], SL.library = learners,
                             verbose = FALSE, method = "method.NNLS",cvControl = control)
    
    score_tau0 <- tau0_mod$SL.predict
    a0 <- score_tau0
    
    
  },error=function(e){
    
    mean_score <- mean(tau0)
    score_tau0 <- rep.int(mean_score, times = nrow(df_aux))
    a0 <- score_tau0
    return(a0)
  })
  
  score_tau0 <- a0
  
  
  
  tau_weight <- p_hat*score_tau0 + (1-p_hat)*score_tau1
  
  
  
  return(tau_weight)
  }

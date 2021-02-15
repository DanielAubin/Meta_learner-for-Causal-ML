R <- 20

est_cate_DR <- matrix(0,nrow(df_main),R)
est_cate_R <- matrix(0,nrow(df_main),R)
est_cate_T <- matrix(0,nrow(df_main),R)
est_cate_X <- matrix(0,nrow(df_main),R)

df_aux <- df_boot

for(r in 1:R){
  set.seed(101+R)

## DR-learner 
# Train a classification model to get the propensity scores
p_mod <- SuperLearner(Y = df_aux$d, X = df_aux[,covariates], newX = df_aux[,covariates], SL.library = learners,
                      verbose = FALSE, method = "method.NNLS", family = binomial(),cvControl = control)

p_hat <- p_mod$SL.predict
p_hat = ifelse(p_hat<0.025, 0.025, ifelse(p_hat>.975,.975, p_hat)) # Overlap bounding

# Prop-Score for df_main
p_hat_main <- predict(p_mod,df_main[,covariates])$pred

# Split the training data into treatment and control observations
aux_1 <- df_aux[which(df_aux$d==1),]
aux_0 <- df_aux[which(df_aux$d==0),]

# Train a regression model for the treatment observations
m1_mod <- SuperLearner(Y = aux_1$y, X = aux_1[,covariates], newX = df_aux[,covariates], SL.library = learners,
                       verbose = FALSE, method = "method.NNLS",cvControl = control)

m1_hat <- m1_mod$SL.predict

# Train a regression model for the control observations
m0_mod <- SuperLearner(Y = aux_0$y, X = aux_0[,covariates], newX = df_aux[,covariates], SL.library = learners,
                       verbose = FALSE, method = "method.NNLS",cvControl = control)

m0_hat <- m0_mod$SL.predict



# Apply the doubly-robust estimator 
y_mo <- (m1_hat - m0_hat) + ((df_aux$d*(df_aux$y -m1_hat))/p_hat) - ((1-df_aux$d)*(df_aux$y - m0_hat)/(1-p_hat))





a  <- tryCatch({
  dr_mod <- SuperLearner(Y = y_mo, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
                         verbose = FALSE, method = "method.NNLS",cvControl = control)
  
  score_dr <- dr_mod$SL.predict
  a <- score_dr
  
  
},error=function(e){
  
  mean_score <- mean(y_mo)
  score_dr <- rep.int(mean_score, times = nrow(df_main))
  a <- score_dr
  return(a)
})

score_dr <- a

######################

est_cate_DR[,r] <- score_dr 
est_cate_T[,r] <- (predict(m1_mod,df_main[,covariates])$pred - predict(m0_mod,df_main[,covariates])$pred)

### R-learner 

# Train a regression model 
m_mod <- SuperLearner(Y = df_aux$y, X = df_aux[,covariates], newX = df_aux[,covariates], SL.library = learners,
                      verbose = FALSE, method = "method.NNLS",cvControl = control)

m_hat <- m_mod$SL.predict

# Apply the R-learner (residual-on-residual approach)
y_tilde = df_aux$y - m_hat
w_tilde = df_aux$d - p_hat
pseudo_outcome = y_tilde/w_tilde

weights = w_tilde^2


a  <- tryCatch({
  R_mod <- SuperLearner(Y = pseudo_outcome, X = df_aux[,covariates], newX = df_main[,covariates], SL.library = learners,
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



###########

est_cate_R[,r] <- score_R

###  X-learner

tau1 <- df_aux[which(df_aux$d==1),"y"] - m0_hat[which(df_aux$d==1),]

tau0 <- m1_hat[which(df_aux$d==0),]  -  df_aux[which(df_aux$d==0),"y"]


a1  <- tryCatch({
  tau1_mod <- SuperLearner(Y = tau1, X = df_aux[which(df_aux$d==1),covariates], newX = df_main[,covariates], SL.library = learners,
                           verbose = FALSE, method = "method.NNLS",cvControl = control)
  
  score_tau1 <- tau1_mod$SL.predict
  a1 <- score_tau1
  
  
},error=function(e){
  
  mean_score <- mean(tau1)
  score_tau1 <- rep.int(mean_score, times = nrow(df_main))
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
  score_tau0 <- rep.int(mean_score, times = nrow(df_main))
  a0 <- score_tau0
  return(a0)
})

score_tau0 <- a0



score_X <- p_hat_main*score_tau0 + (1-p_hat_main)*score_tau1

est_cate_X[,r] <- score_X

cat("This is Iteration: ", r, "out of", R,"\n")

}



cate_DR <- bootCI(pred_B=results_cate_DR,est_R=est_cate_DR)
cate_R <- bootCI(pred_B=results_cate_R,est_R=est_cate_R)
cate_T <- bootCI(pred_B=results_cate_T,est_R=est_cate_T)
cate_X <- bootCI(pred_B=results_cate_X,est_R=est_cate_X)

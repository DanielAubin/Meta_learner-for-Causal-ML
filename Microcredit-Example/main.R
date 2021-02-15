vec.pac= c("foreign", "quantreg", "gbm", "glmnet",
           "MASS", "rpart", "nnet", "matrixStats",
           "xtable", "readstata13","grf","remotes",
           "caret",  "multcomp","cowplot","SuperLearner",
           "ranger","reshape2","gridExtra","bartCause")
#install.packages(vec.pac)
#remotes::install_github("vdorie/bartCause")


lapply(vec.pac, require, character.only = TRUE) 

# Read the microcredit data
data        <- read.dta13("H:/data_rep.dta")


####################################### Inputs  #######################################

data$y <- data$loansamt_total
data$loansamt_total <- NULL

data$d <- data$treatment
data$treatment <- NULL
# "loansamt_total"     vector of outcome variables
# "treatment"   vector of treatment variables



# create a vector of control variables
covariates     <- c("members_resid_bl", "nadults_resid_bl", "head_age_bl", "act_livestock_bl", "act_business_bl", 
                    "borrowed_total_bl", "members_resid_d_bl", "nadults_resid_d_bl", "head_age_d_bl", "act_livestock_d_bl", 
                    "act_business_d_bl", "borrowed_total_d_bl", "ccm_resp_activ", "other_resp_activ", "ccm_resp_activ_d", 
                    "other_resp_activ_d", "head_educ_1", "nmember_age6_16")




######################################################################################################


data <- data[,c("y", "d",covariates)]
# erase some missing values 
data <- data[complete.cases(data),]
data$ID <- c(1:nrow(data))

B <- 100 # Number of Bootstrap repetitions

#Learner Library:

SL.ranger_td = create.Learner("SL.ranger", params = list(num.trees = 1000, min.node.size = 10))
learners <- c( "SL.glmnet","SL.ranger_1")

#CV Control for the SuperLearner
control <- SuperLearner.CV.control(V=10)


# Here we apply our estimates only on a subset of the data and use the training data to generate bootstrap confidence intervals.
set.seed(1011)
folds <- createFolds(data$d,k=5) # Here we use 5-fold sample splitting.

df_boot <- data[c(folds[[5]],folds[[2]],folds[[3]],folds[[4]]),]
df_main <- data[folds[[1]],]



# Create a matrix to store the CATE results from each method 
results_cate_DR <- matrix(0,nrow(df_main),B)
results_cate_R <- matrix(0,nrow(df_main),B)
results_cate_T <- matrix(0,nrow(df_main),B)
results_cate_X <- matrix(0,nrow(df_main),B)






createbootstrappedData <- function(df_boot) {
  
  smpl_0 <- sample((1:nrow(df_boot))[df_boot$d == 0],
                   replace = TRUE,
                   size = sum(1 - df_boot$d))
  smpl_1 <- sample((1:nrow(df_boot))[df_boot$d == 1],
                   replace = TRUE,
                   size = sum(df_boot$d))
  smpl <- sample(c(smpl_0, smpl_1))
  
  return(df_boot[smpl,])
}



for(b in 1:B){
  
  set.seed(1011+b)
  
  ### Apply the meta-learners with Bootstrapping  
  
  df_aux <- createbootstrappedData(df_boot)
  
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
  
  results_cate_DR[,b] <- score_dr 
  results_cate_T[,b] <- (predict(m1_mod,df_main[,covariates])$pred - predict(m0_mod,df_main[,covariates])$pred)
  
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
  
  results_cate_R[,b] <- score_R
  
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
  
  results_cate_X[,b] <- score_X
  
  
  
  
  cat("This is Iteration: ", b, "out of", B,"\n")
}



df_aux <- df_boot
### Causal BART
fit_bart <- bartc(df_aux$y,df_aux$d,df_aux[,covariates],keepTrees = TRUE,ntree=4000)
score_bart <- predict(fit_bart,newdata = df_main,type="icate")
score_bart_m <- apply(score_bart,2,mean)
ite.sd <- apply(score_bart, 2, sd)
ite.lb <- score_bart_m - 1.96 * ite.sd
ite.ub <- score_bart_m + 1.96 * ite.sd

cate_CBART <- as.data.frame(cbind(score_bart_m,ite.lb,ite.ub))
colnames(cate_CBART) <- c("pred","X5.","X95.")

### Causal Forest

forest.D <- regression_forest(df_aux[,covariates], df_aux$d, tune.parameters = "all")
W.hat <- predict(forest.D)$predictions

forest.Y <- regression_forest(df_aux[,covariates], df_aux$y, tune.parameters = "all")
Y.hat <- predict(forest.Y)$predictions

tau.forest <- causal_forest(df_aux[,covariates], df_aux$y, df_aux$d,tune.parameters = "all",num.trees = 8000,
                            min.node.size=5,W.hat = W.hat, Y.hat = Y.hat)
tau.hat <- predict(tau.forest, df_main[,covariates],estimate.variance = TRUE)
summary(tau.hat$predictions)
sigma.hat <- sqrt(tau.hat$variance.estimates)

pred = tau.hat$predictions
X5. =  tau.hat$predictions - 1.96 * sigma.hat
X95. = tau.hat$predictions + 1.96 * sigma.hat

cate_CForest <- as.data.frame(cbind(pred,X5.,X95.))



# Calculate confidence intervals based on bootstrapping 

cate_DR <- bootCI(results_cate_DR)
cate_R <- bootCI(results_cate_R)
cate_T <- bootCI(results_cate_T)
cate_X <- bootCI(results_cate_X)


# Prepare data for plot
prep_plot <- function(cate){
  cate1 <- melt(cate[order(cate$pred),])
  cate1$ID <- rep(1:nrow(cate),times=3)
  return(cate1)
}

# Define some colors
cbp <- c("#000000", "#E69F00", "#0072B2", "#009E73",
         "#F0E442", "#0072B2", "#D55E00", "#CC79A7")




mm_DR <- prep_plot(cate_DR)
mm_R <- prep_plot(cate_R)
mm_T <- prep_plot(cate_T)
mm_X <- prep_plot(cate_X)
mm_CBART <- prep_plot(cate_CBART)
mm_CF <- prep_plot(cate_CForest)

mm <- rbind(mm_DR,mm_R,mm_T,mm_X,mm_CBART,mm_CF)
mm$Method <- rep(c("DR-learner","R-learner","T-learner","X-learner","Causal-BART","Causal-Forest"),each=nrow(mm_DR))

# Plot the CATE estimates + CI for 5% and 95%
# Plot the CATE estimates + CI for 5% and 95%
ggplot(mm, aes(x=ID,y=value,group=variable))+
  geom_line(aes(color=variable,alpha=variable),size=1.1)+
  scale_alpha_manual(values=c(1.0,0.1,0.1),labels = c("CATE","CI lower 5%", "CI upper 95%")) +
  theme_cowplot() +
  facet_wrap( ~ Method, scales="free", nrow=4) + # Facet wrap with common scales 
  labs(y = "Treatment Effect ", x = "Ordered Observation") +
  theme(legend.position="none", legend.justification = 'center')+ 
  scale_color_manual(values = cbp,labels = c("CATE","CI lower 5%", "CI upper 95%")) +
  geom_smooth(data=mm[mm["variable"] == "X5." | mm["variable"] == "X95.",],linetype="dashed", size=0.5) +
  ylim(-1000,5000)


### Estimate 20% least, ATE and 80% most affected


quantiles <- function(CATE){
S2        <- CATE+runif(length(CATE), 0, 0.00001) # Include white noise to guarantee that the score (S) differs from the baseline effect
S_20 <- as.numeric(quantile(S2, c(.20)))
S_ATE <- mean(S2)
S_80 <- as.numeric(quantile(S2,c(.80)))

return(data.frame(S_20,S_ATE,S_80))

}



cate_list <- list(cate_DR$pred,cate_R$pred,cate_T$pred,cate_X$pred,cate_CBART$pred,cate_CForest$pred)
res <- lapply(cate_list, quantiles)
res


### Classification Analysis




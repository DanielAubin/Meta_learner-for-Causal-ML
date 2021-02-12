install.packages("remotes")
remotes::install_github("vdorie/bartCause")

# Fits a collection of treatment and response models using the Bayesian Additive
# Regression Trees (BART) algorithm, producing estimates of treatment effects.


Causal_BART <- function(df_aux,df_main,covariates){

fit_bart <- bartc(df_aux$y,df_aux$d,df_aux[,covariates],keepTrees = TRUE)
score_bart <- predict(fit_bart,newdata = df_main,type="icate")
score_bart_m <- apply(score_bart,2,mean)
ite.sd <- apply(score_bart, 2, sd)
ite.lb <- score_bart_m - 2 * ite.sd
ite.ub <- score_bart_m + 2 * ite.sd

cate_CBART <- as.data.frame(cbind(score_bart_m,ite.lb,ite.ub))
colnames(cate_CBART) <- c("pred","X5.","X95.")

return(cate_CBART)

}


install.packages("remotes")
remotes::install_github("vdorie/bartCause")

# Fits a collection of treatment and response models using the Bayesian Additive
# Regression Trees (BART) algorithm, producing estimates of treatment effects.

Causal_BART <- function(y,d,x){

fit_bart <- bartc(y,d,x, )
score_bart <- extract(fit_bart,type="icate")
score_bart_m <- apply(score_bart,2,mean)

      #ite.sd <- apply(score_bart, 2, sd)
      #ite.lb <- score_bart_m - 1.96 * ite.sd
      #ite.ub <- score_bart_m + 1.96 * ite.sd

      #cate_CBART <- as.data.frame(cbind(score_bart_m,ite.lb,ite.ub))
      #colnames(cate_CBART) <- c("pred","X5.","X95.")

     return(score_bart_m)
 
}

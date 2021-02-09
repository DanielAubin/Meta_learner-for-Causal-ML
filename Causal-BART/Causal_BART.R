install.packages("remotes")
remotes::install_github("vdorie/bartCause")

# Fits a collection of treatment and response models using the Bayesian Additive
# Regression Trees (BART) algorithm, producing estimates of treatment effects.


bartc(response, treatment, confounders, data, subset, weights,
      method.rsp = c("bart", "tmle", "p.weight"),
      method.trt = c("bart", "glm", "none"),
      estimand   = c("ate", "att", "atc"),
      group.by = NULL,
      commonSup.rule = c("none", "sd", "chisq"),
      commonSup.cut  = c(NA_real_, 1, 0.05),
      args.rsp = list(), args.trt = list(),
      p.scoreAsCovariate = TRUE, use.ranef = TRUE, group.effects = FALSE,
      crossvalidate = FALSE,
      keepCall = TRUE, verbose = TRUE,
      \dots)
}
# response:
# A vector of the outcome variable, or a reference to such in the \code{data}
# argument. Can be continuous or binary.

# treatment:
# A vector of the binary treatment variable, or a reference to \code{data}.

# confounders:
# A matrix or data frame of covariates to be used in estimating the treatment
# and response model. Can also be the right-hand-side of a formula (e.g.
# \code{x1 + x2 + ...}). The \code{data} argument will be searched if
# supplied.
  

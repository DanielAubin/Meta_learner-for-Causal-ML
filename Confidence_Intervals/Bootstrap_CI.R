

bootCI <- function(pred_B){
  

  # the the 5% and 95% CI from the bootstrapped procedure
  CI_b <- data.frame(
    X5. =  apply(pred_B, 1, function(x)
      quantile(x, c(.025))),
    X95. = apply(pred_B, 1, function(x)
      quantile(x, c(.975))),
    sd = apply(pred_B, 1, function(x) sd(x))
  )
  
  return(data.frame(
    pred = apply(pred_B,1,mean),
    X5. =  apply(pred_B,1,mean) - 1.96 * CI_b$sd,
    X95. = apply(pred_B,1,mean) + 1.96 * CI_b$sd
  ))
}


n <- 1000 # Amount of observations
X1 <- rnorm(n,0,1) # Distribution of characteristic
D <- rbinom(n,1,0.5) # Random treatment assigment
tau <- ifelse(X<0,4,10) # Treatment effect: If X<0, tau = 4; if X>0, tau = 10
Y <- D*tau + rnorm(n,0,0.5) # Definde observed outcome

tau_hat_smallX <- mean(Y[D==1 & X<0]) - mean(Y[D==0 & X<0]) # ATE given X<0
tau_hat_bigX <- mean(Y[D==1 & X>0]) - mean(Y[D==0 & X>0]) # ATE given X<0

print(c(tau_hat_smallX,tau_hat_bigX)) # Print CATE

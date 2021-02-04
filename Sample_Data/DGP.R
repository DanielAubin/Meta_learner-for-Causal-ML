# Generate data.
n = 2000; p = 20
X = matrix(rnorm(n * p), n, p)
theta = 1 / (1 + exp(-X[, 3]))
d = rbinom(n ,1, 1 / (1 + exp(-X[, 1] - X[, 2])))
y = pmax(X[, 2] + X[, 3], 0) + rowMeans(X[, 4:6]) / 2 + d * theta + rnorm(n)


data <- cbind(y,theta,d,X)
colnames(data) <- c("y", "theta", "d", c(paste0("V", 1:p)))
head(data)

# split data in df_aux and df_main
G <- 2 # Amout of folds 
split             <- runif(nrow(data))
cvgroup           <- as.numeric(cut(split,quantile(split,probs = seq(0, 1, 1/G)),include.lowest = TRUE))  

df_aux      <- as.data.frame(data[cvgroup == 1,])
df_main       <- as.data.frame(data[cvgroup != 1,])  



k <- ncol(data)-3
covariates <- c(paste0("V", 1:k))

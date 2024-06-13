# confidence interval for lambda from a stochastic population matrix?
library(popdemo)
library(Hmisc)
data(Pbear)
# project over 50 years with uniform matrix selection
Pbearvec <- c(0.106, 0.068, 0.106, 0.461, 0.151, 0.108)
p2 <- project(Pbear, Pbearvec, time = 50, Aseq = "unif")

# list of matrices used for projection years 
lmat_stoch <- popdemo::mat(p2)[popdemo::Aseq(p2)]
# lambda values
stoch_lambda_vals <- unlist(lapply(lmat_stoch, popdemo::eigs, what = "lambda"))
# summaries
stoch_lambda_mean <- mean(stoch_lambda_vals)
# stoch_lambda_med <- median(stoch_lambda_vals) # No. Suggests pop is growing but is not.
stoch_lambda_quant <- quantile(stoch_lambda_vals, probs = c(0.025, 0.975))
stoch_lambda_boot <- Hmisc::smean.cl.boot(stoch_lambda_vals)
stoch_lambda_boot_mean <- stoch_lambda_boot["Mean"]
stoch_lambda_boot_lcl <- stoch_lambda_boot["Lower"]
stoch_lambda_boot_ucl <- stoch_lambda_boot["Upper"]

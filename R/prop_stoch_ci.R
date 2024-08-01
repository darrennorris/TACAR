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

### now same for Rage - generation time etc
# split into survival and reproductive
testm <- lmat_stoch[[1]]
split_mat <- function(x){ 
  mat_u <- x
  mat_u[1, 6] <- 0
  mat_f <- x 
  mat_f[mat_f > 0] <- 0
  mat_f[1, 6] <- x[1, 6]
  mats <- list(mat_u = mat_u, mat_f = mat_f)
  return(mats)
    }
matlist_u_f <- lapply(lmat_stoch, split_mat)
# wrapper around Rage::gen_time function
my_gen <- function(x){ 
  gen_age_diff <- Rage::gen_time(matU = x$mat_u, x$mat_f, method = "age_diff") 
  return(gen_age_diff) 
  }
gen_vals <- unlist(lapply(matlist_u_f, my_gen))


# Now with real data
# make population matrix
x <- dt[3, ]
stage_names <- c("a1", "a2", "a3", "a4",
                 "b1", "b2", "b3", "b4",
                 "c1", "c2", "c3", "c4",
                 "d1", "d2", "d3", "d4")
stoch_stage_names <- "b1"
year_val <- x[ , stoch_stage_names]
#bad_year_thresh <- 0.1
#good_year_thresh <- 0.2
din <-  x[1 , stage_names] 
# make vector of stochastic values
vb1 <- c(0.0001, 0.1, rep(year_val, 8))
row_n <- length(vb1)
dfstoch <- din[rep(seq_len(nrow(din)), row_n), ]
dfstoch$b1 <- vb1

#list of matrices for projection
makemat <- function(x){
  vpop <- unlist(x[ , stage_names])
  pop_mat <- matrix(vpop, byrow = TRUE, ncol=4)
  dimnames(pop_mat) <- list(c("a", "b", "c", "d"),
                            c(1,2,3,4))
  mout <- pop_mat
  return(mout)
}

# Make list 
matlist <- plyr::alply(dfstoch , .margins = 1, .fun = makemat)
pop_n <-  x$adultF_n * c(11.1, 4, 2, 1) 
ttime <- 100  
spr_01 <- popdemo::project(matlist, vector = pop_n, Aseq = "unif", time = ttime)
spr_01_lmat <- popdemo::mat(spr_01)[popdemo::Aseq(spr_01)]
# Now split into survival and reproduction
split_mat <- function(x){ 
  mat_u <- x
  mat_u[1, 4] <- 0
  mat_f <- x 
  mat_f[mat_f > 0] <- 0
  mat_f[1, 4] <- x[1, 4]
  mats <- list(mat_u = mat_u, mat_f = mat_f)
  return(mats)
}
spr_01_lmat_u_f <- lapply(unname(spr_01_lmat), split_mat)
my_gen <- function(x){ 
  gen_age_diff <- Rage::gen_time(matU = x$mat_u, x$mat_f, method = "age_diff") 
  return(gen_age_diff) 
}
spr_01_gen_age_diff <- unlist(lapply(spr_01_lmat_u_f, my_gen))
mean(as.numeric(spr_01_gen_age_diff))
median(as.numeric(spr_01_gen_age_diff))
quantile(spr_01_gen_age_diff, probs = seq(0, 1, 0.25), 
         na.rm = TRUE)
my_rage <- function(x){ 
  gen_age_diff <- Rage::gen_time(matU = x$mat_u, x$mat_f, method = "age_diff") 
  # Life expectancy (time to death)
  life_exp_stage <- Rage::life_expect_mean(matU = x$mat_u, start = NULL)
  if(life_exp_stage[1] > 1){life_exp_stage[1] <- 1}
  life_exp <- sum(life_exp_stage) 
  life_exp_adult <- life_exp_stage[4]
  dfout <- data.frame(gen_age_diff = gen_age_diff, life_exp = life_exp, 
             life_exp_adult = life_exp_adult) 
  dfout
}
# Apply to each matrix
spr_01_rage <- plyr::ldply(spr_01_lmat_u_f, my_rage)

names(spr_01_lmat_u_f)

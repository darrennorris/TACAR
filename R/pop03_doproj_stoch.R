#' @title Create data.frame with projection results.
#' 
#' @description Creates data.frame with projection results for use in subsequent
#' population projection function ("PopProj.R"). 
#'
#' @param x data.frame created by "PopPrep.R".
#'
#' @details This function applies stochastic population projection moodels. 
#' 
#' @return Creates a data.frame with projection results.
#' @export
#'
#' @examples
#' \dontrun{
#' dfproj_vals <- pop03_doproj_stoch()
#' }
#' 

pop03_doproj_stoch <- function(x) {
  # make population matrix
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
  
   # Projection time and population vector
   pop_n <-  x$adultF_n * c(11.1, 4, 2, 1) 
   ttime <- 100  
   
   # Scenarios
   # Default, the projection selects each matrix with an equal 
   # probability at each time interval.
   # Bad years at 2/5 = 0.4 chance of being chosen
   #set.seed(c(1,2,3,4))
   # Uniform - random
   spr_01 <- popdemo::project(matlist, vector = pop_n, Aseq = "unif", time = ttime)
   # Bad and good years have equal probability
   p2 <- 0.5
   m2 <- m2 <- matrix(rep(c(p2/2, p2/2, (1-p2)/8, (1-p2)/8, (1-p2)/8, 
                            (1-p2)/8, (1-p2)/8, (1-p2)/8, (1-p2)/8, (1-p2)/8), 10), 10, 10)
   spr_02 <- popdemo::project(matlist, vector = pop_n, Aseq = m2, time = ttime)
   # Bad years occur twice as often as good years.
   p3 <- 0.667
   m3 <- matrix(rep(c(p3/2, p3/2, (1-p3)/8, (1-p3)/8, (1-p3)/8, 
                      (1-p3)/8, (1-p3)/8, (1-p3)/8, (1-p3)/8, (1-p3)/8), 10), 10, 10)
   spr_03 <- popdemo::project(matlist, vector= pop_n, Aseq = m3, time = ttime)
   # Bad years now occur four times as often as good years.
   p4 <- 0.8
   m4 <- matrix(rep(c(p4/2, p4/2, (1-p4)/8, (1-p4)/8, (1-p4)/8, 
                      (1-p4)/8, (1-p4)/8, (1-p4)/8, (1-p4)/8, (1-p4)/8), 10), 10, 10)
   spr_04 <- popdemo::project(matlist, vector= pop_n, Aseq = m4, time = ttime)
   
   # Stochastic projection matrices
   # List of matrices used in stochastic projection years 
   spr_01_lmat <- popdemo::mat(spr_01)[popdemo::Aseq(spr_01)]
   spr_02_lmat <- popdemo::mat(spr_02)[popdemo::Aseq(spr_02)]
   spr_03_lmat <- popdemo::mat(spr_03)[popdemo::Aseq(spr_03)]
   spr_04_lmat <- popdemo::mat(spr_04)[popdemo::Aseq(spr_04)]
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
   spr_02_lmat_u_f <- lapply(unname(spr_02_lmat), split_mat)
   spr_03_lmat_u_f <- lapply(unname(spr_03_lmat), split_mat)
   spr_04_lmat_u_f <- lapply(unname(spr_04_lmat), split_mat)
   
   # Demographic and life histroy traits from matrices
   # lambda values
   spr_01_lambda_vals <- unlist(lapply(spr_01_lmat, popdemo::eigs, what = "lambda"))
   spr_02_lambda_vals <- unlist(lapply(spr_02_lmat, popdemo::eigs, what = "lambda"))
   spr_03_lambda_vals <- unlist(lapply(spr_03_lmat, popdemo::eigs, what = "lambda"))
   spr_04_lambda_vals <- unlist(lapply(spr_04_lmat, popdemo::eigs, what = "lambda"))
   # generation times
   spr_01_gen_vals <- unlist(lapply(spr_01_lmat, popbio::generation.time)) 
   spr_02_gen_vals <- unlist(lapply(spr_02_lmat, popbio::generation.time)) 
   spr_03_gen_vals <- unlist(lapply(spr_03_lmat, popbio::generation.time)) 
   spr_04_gen_vals <- unlist(lapply(spr_04_lmat, popbio::generation.time)) 
   
   # Wrapper around Rage functions.
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
   spr_02_rage <- plyr::ldply(spr_02_lmat_u_f, my_rage)
   spr_03_rage <- plyr::ldply(spr_03_lmat_u_f, my_rage)
   spr_04_rage <- plyr::ldply(spr_04_lmat_u_f, my_rage)

   # summaries
   #as.numeric to remove names which generates warning
   spr_01_boot <- Hmisc::smean.cl.boot(spr_01_lambda_vals)
   spr_01_boot_mean <- as.numeric(spr_01_boot["Mean"])
   spr_01_boot_lcl <- as.numeric(spr_01_boot["Lower"])
   spr_01_boot_ucl <- as.numeric(spr_01_boot["Upper"])
   spr_01_quant <- quantile(spr_01_lambda_vals, probs = seq(0, 1, 0.25), 
                            na.rm = TRUE)
   spr_01_q25 <- as.numeric(spr_01_quant["25%"])
   spr_01_q75 <- as.numeric(spr_01_quant["75%"])
   spr_01_sd <- sd(spr_01_lambda_vals, na.rm = TRUE)
   spr_02_boot <- Hmisc::smean.cl.boot(spr_02_lambda_vals)
   spr_02_boot_mean <- as.numeric(spr_02_boot["Mean"])
   spr_02_boot_lcl <- as.numeric(spr_02_boot["Lower"])
   spr_02_boot_ucl <- as.numeric(spr_02_boot["Upper"]) 
   spr_02_quant <- quantile(spr_02_lambda_vals, probs = seq(0, 1, 0.25), 
                            na.rm = TRUE)
   spr_02_q25 <- as.numeric(spr_02_quant["25%"])
   spr_02_q75 <- as.numeric(spr_02_quant["75%"])
   spr_02_sd <- sd(spr_02_lambda_vals, na.rm = TRUE)
   spr_03_boot <- Hmisc::smean.cl.boot(spr_03_lambda_vals)
   spr_03_boot_mean <- as.numeric(spr_03_boot["Mean"])
   spr_03_boot_lcl <- as.numeric(spr_03_boot["Lower"])
   spr_03_boot_ucl <- as.numeric(spr_03_boot["Upper"]) 
   spr_03_quant <- quantile(spr_03_lambda_vals, probs = seq(0, 1, 0.25), 
                            na.rm = TRUE)
   spr_03_q25 <- as.numeric(spr_03_quant["25%"])
   spr_03_q75 <- as.numeric(spr_03_quant["75%"])
   spr_03_sd <- sd(spr_03_lambda_vals, na.rm = TRUE)
   spr_04_boot <- Hmisc::smean.cl.boot(spr_04_lambda_vals)
   spr_04_boot_mean <- as.numeric(spr_04_boot["Mean"])
   spr_04_boot_lcl <- as.numeric(spr_04_boot["Lower"])
   spr_04_boot_ucl <- as.numeric(spr_04_boot["Upper"]) 
   spr_04_quant <- quantile(spr_04_lambda_vals, probs = seq(0, 1, 0.25), 
                            na.rm = TRUE)
   spr_04_q25 <- as.numeric(spr_04_quant["25%"])
   spr_04_q75 <- as.numeric(spr_04_quant["75%"])
   spr_04_sd <- sd(spr_04_lambda_vals, na.rm = TRUE)
   # generation times
   spr_01_gen_boot <- Hmisc::smean.cl.boot(spr_01_gen_vals)
   spr_01_gen_boot_mean <- as.numeric(spr_01_gen_boot["Mean"])
   spr_01_gen_boot_lcl <- as.numeric(spr_01_gen_boot["Lower"])
   spr_01_gen_boot_ucl <- as.numeric(spr_01_gen_boot["Upper"])
   spr_01_gen_sd <- sd(spr_01_gen_vals, na.rm = TRUE) 
   spr_01_gen_quant <- quantile(spr_01_gen_vals, probs = seq(0, 1, 0.25), 
                            na.rm = TRUE)
   spr_01_gen_q25 <- as.numeric(spr_01_gen_quant["25%"])
   spr_01_gen_q75 <- as.numeric(spr_01_gen_quant["75%"])
   spr_02_gen_boot <- Hmisc::smean.cl.boot(spr_02_gen_vals)
   spr_02_gen_boot_mean <- as.numeric(spr_02_gen_boot["Mean"])
   spr_02_gen_boot_lcl <- as.numeric(spr_02_gen_boot["Lower"])
   spr_02_gen_boot_ucl <- as.numeric(spr_02_gen_boot["Upper"]) 
   spr_02_gen_sd <- sd(spr_02_gen_vals, na.rm = TRUE) 
   spr_02_gen_quant <- quantile(spr_02_gen_vals, probs = seq(0, 1, 0.25), 
                                na.rm = TRUE)
   spr_02_gen_q25 <- as.numeric(spr_02_gen_quant["25%"])
   spr_02_gen_q75 <- as.numeric(spr_02_gen_quant["75%"])
   spr_03_gen_boot <- Hmisc::smean.cl.boot(spr_03_gen_vals)
   spr_03_gen_boot_mean <- as.numeric(spr_03_gen_boot["Mean"])
   spr_03_gen_boot_lcl <- as.numeric(spr_03_gen_boot["Lower"])
   spr_03_gen_boot_ucl <- as.numeric(spr_03_gen_boot["Upper"])
   spr_03_gen_sd <- sd(spr_03_gen_vals, na.rm = TRUE) 
   spr_03_gen_quant <- quantile(spr_03_gen_vals, probs = seq(0, 1, 0.25), 
                                na.rm = TRUE)
   spr_03_gen_q25 <- as.numeric(spr_03_gen_quant["25%"])
   spr_03_gen_q75 <- as.numeric(spr_03_gen_quant["75%"])
   spr_04_gen_boot <- Hmisc::smean.cl.boot(spr_04_gen_vals)
   spr_04_gen_boot_mean <- as.numeric(spr_04_gen_boot["Mean"])
   spr_04_gen_boot_lcl <- as.numeric(spr_04_gen_boot["Lower"])
   spr_04_gen_boot_ucl <- as.numeric(spr_04_gen_boot["Upper"]) 
   spr_04_gen_sd <- sd(spr_04_gen_vals, na.rm = TRUE) 
   spr_04_gen_quant <- quantile(spr_04_gen_vals, probs = seq(0, 1, 0.25), 
                                na.rm = TRUE)
   spr_04_gen_q25 <- as.numeric(spr_04_gen_quant["25%"])
   spr_04_gen_q75 <- as.numeric(spr_04_gen_quant["75%"])
   
   # Generation time based on age difference. Get some huge values for
   # 2% cases when population is strongly declining.
   spr_01_gen_age_diff_med <- median(as.numeric(spr_01_rage$gen_age_diff))
   spr_01_gen_age_quant <- quantile(as.numeric(spr_01_rage$gen_age_diff), probs = seq(0, 1, 0.25), 
                                na.rm = TRUE)
   spr_01_gen_age_q25 <- as.numeric(spr_01_gen_age_quant["25%"])
   spr_01_gen_age_q75 <- as.numeric(spr_01_gen_age_quant["75%"])
   spr_02_gen_age_diff_med <- median(as.numeric(spr_02_rage$gen_age_diff))
   spr_02_gen_age_quant <- quantile(as.numeric(spr_02_rage$gen_age_diff), probs = seq(0, 1, 0.25), 
                                    na.rm = TRUE)
   spr_02_gen_age_q25 <- as.numeric(spr_02_gen_age_quant["25%"])
   spr_02_gen_age_q75 <- as.numeric(spr_02_gen_age_quant["75%"])
   spr_03_gen_age_diff_med <- median(as.numeric(spr_03_rage$gen_age_diff))
   spr_03_gen_age_quant <- quantile(as.numeric(spr_03_rage$gen_age_diff), probs = seq(0, 1, 0.25), 
                                    na.rm = TRUE)
   spr_03_gen_age_q25 <- as.numeric(spr_03_gen_age_quant["25%"])
   spr_03_gen_age_q75 <- as.numeric(spr_03_gen_age_quant["75%"])
   spr_04_gen_age_diff_med <- median(as.numeric(spr_01_rage$gen_age_diff))
   spr_04_gen_age_quant <- quantile(as.numeric(spr_01_rage$gen_age_diff), probs = seq(0, 1, 0.25), 
                                    na.rm = TRUE)
   spr_04_gen_age_q25 <- as.numeric(spr_04_gen_age_quant["25%"])
   spr_04_gen_age_q75 <- as.numeric(spr_04_gen_age_quant["75%"])
   
   # Life expectancy
   spr_01_life_exp_med <- median(as.numeric(spr_01_rage$life_exp))
   spr_01_life_exp_quant <- quantile(as.numeric(spr_01_rage$life_exp), probs = seq(0, 1, 0.25), 
                                    na.rm = TRUE)
   spr_01_life_exp_q25 <- as.numeric(spr_01_life_exp_quant["25%"])
   spr_01_life_exp_q75 <- as.numeric(spr_01_life_exp_quant["75%"])
   spr_01_life_exp_adult_med <- median(as.numeric(spr_01_rage$life_exp_adult))
   spr_01_life_exp_adult_quant <- quantile(as.numeric(spr_01_rage$life_exp_adult), probs = seq(0, 1, 0.25), 
                                     na.rm = TRUE)
   spr_01_life_exp_adult_q25 <- as.numeric(spr_01_life_exp_adult_quant["25%"])
   spr_01_life_exp_adult_q75 <- as.numeric(spr_01_life_exp_adult_quant["75%"])
   spr_02_life_exp_med <- median(as.numeric(spr_02_rage$life_exp))
   spr_02_life_exp_quant <- quantile(as.numeric(spr_02_rage$life_exp), probs = seq(0, 1, 0.25), 
                                     na.rm = TRUE)
   spr_02_life_exp_q25 <- as.numeric(spr_02_life_exp_quant["25%"])
   spr_02_life_exp_q75 <- as.numeric(spr_02_life_exp_quant["75%"])
   spr_02_life_exp_adult_med <- median(as.numeric(spr_02_rage$life_exp_adult))
   spr_02_life_exp_adult_quant <- quantile(as.numeric(spr_02_rage$life_exp_adult), probs = seq(0, 1, 0.25), 
                                           na.rm = TRUE)
   spr_02_life_exp_adult_q25 <- as.numeric(spr_02_life_exp_adult_quant["25%"])
   spr_02_life_exp_adult_q75 <- as.numeric(spr_02_life_exp_adult_quant["75%"]) 
   spr_03_life_exp_med <- median(as.numeric(spr_03_rage$life_exp))
   spr_03_life_exp_quant <- quantile(as.numeric(spr_03_rage$life_exp), probs = seq(0, 1, 0.25), 
                                     na.rm = TRUE)
   spr_03_life_exp_q25 <- as.numeric(spr_03_life_exp_quant["25%"])
   spr_03_life_exp_q75 <- as.numeric(spr_03_life_exp_quant["75%"])
   spr_03_life_exp_adult_med <- median(as.numeric(spr_03_rage$life_exp_adult))
   spr_03_life_exp_adult_quant <- quantile(as.numeric(spr_03_rage$life_exp_adult), probs = seq(0, 1, 0.25), 
                                           na.rm = TRUE)
   spr_03_life_exp_adult_q25 <- as.numeric(spr_03_life_exp_adult_quant["25%"])
   spr_03_life_exp_adult_q75 <- as.numeric(spr_03_life_exp_adult_quant["75%"])
   
   spr_04_life_exp_med <- median(as.numeric(spr_04_rage$life_exp))
   spr_04_life_exp_quant <- quantile(as.numeric(spr_04_rage$life_exp), probs = seq(0, 1, 0.25), 
                                     na.rm = TRUE)
   spr_04_life_exp_q25 <- as.numeric(spr_04_life_exp_quant["25%"])
   spr_04_life_exp_q75 <- as.numeric(spr_04_life_exp_quant["75%"])
   spr_04_life_exp_adult_med <- median(as.numeric(spr_04_rage$life_exp_adult))
   spr_04_life_exp_adult_quant <- quantile(as.numeric(spr_04_rage$life_exp_adult), probs = seq(0, 1, 0.25), 
                                           na.rm = TRUE)
   spr_04_life_exp_adult_q25 <- as.numeric(spr_04_life_exp_adult_quant["25%"])
   spr_04_life_exp_adult_q75 <- as.numeric(spr_04_life_exp_adult_quant["75%"])
   
   # data.frame to return
   dfpop <- rbind(data.frame(model = "Stochastic uniform",  
                             lambda = spr_01_boot_mean,
                             lambda_lcl = spr_01_boot_lcl, 
                             lambda_ucl = spr_01_boot_ucl, 
                             lambda_sd = spr_01_sd, 
                             lambda_q25 = spr_01_q25, 
                             lambda_q75 = spr_01_q75,
                             gen_time = spr_01_gen_boot_mean, 
                             gen_sd = spr_01_gen_sd, 
                             gen_q25 = spr_01_gen_q25, 
                             gen_q75 = spr_01_gen_q75, 
                             gen_age_diff_med = spr_01_gen_age_diff_med, 
                             gen_age_q25 = spr_01_gen_age_q25, 
                             gen_age_q75 = spr_01_gen_age_q75, 
                             life_exp_med = spr_01_life_exp_med, 
                             life_exp_q25 = spr_01_life_exp_q25, 
                             life_exp_q75 = spr_01_life_exp_q75, 
                             life_exp_adult_med = spr_01_life_exp_adult_med, 
                             life_exp_adult_q25 = spr_01_life_exp_adult_q25, 
                             life_exp_adult_q75 = spr_01_life_exp_adult_q75,
                             ayear = 0:ttime,
                       egghatch = popdemo::vec(spr_01)[,1], 
                       j_early = popdemo::vec(spr_01)[,2], 
                       j_late = popdemo::vec(spr_01)[,3], 
                       fem = popdemo::vec(spr_01)[,4], 
                       fem_t0 = as.numeric(popdemo::vec(spr_01)[1, 4]),
                pop = spr_01), 
                data.frame(model = "Stochastic equal", 
                           lambda = spr_02_boot_mean,
                           lambda_lcl = spr_02_boot_lcl, 
                           lambda_ucl = spr_02_boot_ucl, 
                           lambda_sd = spr_02_sd, 
                           lambda_q25 = spr_02_q25, 
                           lambda_q75 = spr_02_q75,
                           gen_time = spr_02_gen_boot_mean, 
                           gen_sd = spr_02_gen_sd, 
                           gen_q25 = spr_02_gen_q25, 
                           gen_q75 = spr_02_gen_q75, 
                           gen_age_diff_med = spr_02_gen_age_diff_med, 
                           gen_age_q25 = spr_02_gen_age_q25, 
                           gen_age_q75 = spr_02_gen_age_q75, 
                           life_exp_med = spr_02_life_exp_med, 
                           life_exp_q25 = spr_02_life_exp_q25, 
                           life_exp_q75 = spr_02_life_exp_q75, 
                           life_exp_adult_med = spr_02_life_exp_adult_med, 
                           life_exp_adult_q25 = spr_02_life_exp_adult_q25, 
                           life_exp_adult_q75 = spr_02_life_exp_adult_q75,
                           ayear = 0:ttime, 
                           egghatch = popdemo::vec(spr_02)[,1], 
                           j_early = popdemo::vec(spr_02)[,2], 
                           j_late = popdemo::vec(spr_02)[,3], 
                           fem = popdemo::vec(spr_02)[,4], 
                           fem_t0 = as.numeric(popdemo::vec(spr_02)[1, 4]),
                           pop = spr_02),
                data.frame(model = "Stochastic bad x2", 
                           lambda = spr_03_boot_mean,
                           lambda_lcl = spr_03_boot_lcl, 
                           lambda_ucl = spr_03_boot_ucl, 
                           lambda_sd = spr_03_sd, 
                           lambda_q25 = spr_03_q25, 
                           lambda_q75 = spr_03_q75,
                           gen_time = spr_03_gen_boot_mean, 
                           gen_sd = spr_03_gen_sd, 
                           gen_q25 = spr_03_gen_q25, 
                           gen_q75 = spr_03_gen_q75, 
                           gen_age_diff_med = spr_03_gen_age_diff_med, 
                           gen_age_q25 = spr_03_gen_age_q25, 
                           gen_age_q75 = spr_03_gen_age_q75, 
                           life_exp_med = spr_03_life_exp_med, 
                           life_exp_q25 = spr_03_life_exp_q25, 
                           life_exp_q75 = spr_03_life_exp_q75, 
                           life_exp_adult_med = spr_03_life_exp_adult_med, 
                           life_exp_adult_q25 = spr_03_life_exp_adult_q25, 
                           life_exp_adult_q75 = spr_03_life_exp_adult_q75,
                           ayear = 0:ttime, 
                           egghatch = popdemo::vec(spr_03)[,1], 
                           j_early = popdemo::vec(spr_03)[,2], 
                           j_late = popdemo::vec(spr_03)[,3], 
                           fem = popdemo::vec(spr_03)[,4], 
                           fem_t0 = as.numeric(popdemo::vec(spr_03)[1, 4]),
                           pop = spr_03),
                data.frame(model = "Stochastic bad x4", 
                           lambda = spr_04_boot_mean,
                           lambda_lcl = spr_04_boot_lcl, 
                           lambda_ucl = spr_04_boot_ucl, 
                           lambda_sd = spr_04_sd, 
                           lambda_q25 = spr_04_q25, 
                           lambda_q75 = spr_04_q75,
                           gen_time = spr_04_gen_boot_mean, 
                           gen_sd = spr_04_gen_sd, 
                           gen_q25 = spr_04_gen_q25, 
                           gen_q75 = spr_04_gen_q75, 
                           gen_age_diff_med = spr_04_gen_age_diff_med, 
                           gen_age_q25 = spr_04_gen_age_q25, 
                           gen_age_q75 = spr_04_gen_age_q75, 
                           life_exp_med = spr_04_life_exp_med, 
                           life_exp_q25 = spr_04_life_exp_q25, 
                           life_exp_q75 = spr_04_life_exp_q75, 
                           life_exp_adult_med = spr_04_life_exp_adult_med, 
                           life_exp_adult_q25 = spr_04_life_exp_adult_q25, 
                           life_exp_adult_q75 = spr_04_life_exp_adult_q75,
                           ayear = 0:ttime, 
                           egghatch = popdemo::vec(spr_04)[,1], 
                           j_early = popdemo::vec(spr_04)[,2], 
                           j_late = popdemo::vec(spr_04)[,3], 
                           fem = popdemo::vec(spr_04)[,4], 
                           fem_t0 = as.numeric(popdemo::vec(spr_04)[1, 4]),
                           pop = spr_04) 
   )
   dfpop$fem_diff <- round(((dfpop$fem - dfpop$fem_t0) / dfpop$fem_t0), 3)
   dfpop$change50_flag <- as.integer(ifelse(abs(dfpop$fem_diff) >= 0.5, 1, 0))
   dfpop$change30_flag <- as.integer(ifelse(abs(dfpop$fem_diff) >= 0.3, 1, 0))
   dfpop$double_flag <- as.integer(ifelse(dfpop$fem_diff >= 1.0, 1, 0))         
   return(dfpop)
}
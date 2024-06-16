#' @title Create data.frame with projection results.
#' 
#' @description Creates data.frame with projection results for use in subsequent
#' population projection function ("PopProj.R"). 
#'
#' @param x data.frame created by "PopPrep.R".
#'
#' @details This function applies a deterministic population projection moodel. 
#' 
#' @return Creates data.frame with projection results.
#' @export
#'
#' @examples
#' \dontrun{
#' dfproj_vals <- pop03_doproj()
#' }
#' 

pop03_doproj <- function(x) {
  pop_n <-  x$adultF_n * c(11.1, 4, 2, 1) 
  stage_names <- c("a1", "a2", "a3", "a4",
                   "b1", "b2", "b3", "b4",
                   "c1", "c2", "c3", "c4",
                   "d1", "d2", "d3", "d4")
  vpop <- unlist(x[ , stage_names])
  pop_mat <- matrix(vpop, byrow = TRUE, ncol=4)
  dimnames(pop_mat) <- list(c("a", "b", "c", "d"),
                            c(1,2,3,4))
  # project PPM 
  ttime <- 100
  pr_pop <- popdemo::project(pop_mat, vector = pop_n, time = ttime)
  
  # data for plotting
  len <- length(pr_pop)
  Time.intervals <- 0:(len - 1)
  eggs <- as.integer(trunc(pr_pop * (popbio::stable.stage(pop_mat)[1])))
  eju <- as.integer(trunc(pr_pop * (popbio::stable.stage(pop_mat)[2])))
  lju <- as.integer(trunc(pr_pop * (popbio::stable.stage(pop_mat)[3])))
  ad.fe <- as.integer(trunc(pr_pop * (popbio::stable.stage(pop_mat)[4])))
  # below avoids integer overflow on large numbers
  all_inds <- as.numeric(eggs) + as.numeric(eju) + as.numeric(lju) + as.numeric(ad.fe)
  plambda <- popbio::lambda(pop_mat)
  flag_ergo <- popdemo::isErgodic(pop_mat)
  flag_irre <- popdemo::isIrreducible(pop_mat)
  gen_time <- popbio::generation.time(pop_mat)
  # rep_val_f <- popdemo::eigs(pop_mat, "rv")[4] # reproductive value
  # make dataframe 
  dfout <- data.frame(model = "Deterministic",
                      lambda = plambda, 
                      gen_time = gen_time,
                      ergodic = flag_ergo,
                      irred = flag_irre,
                      ayear = Time.intervals, 
                      individuals = all_inds,
                      ss_egghatchling = round(as.numeric(popbio::stable.stage(pop_mat)[1]),3),
                      ss_earlyjuven = round(as.numeric(popbio::stable.stage(pop_mat)[2]),3),
                      ss_latejuven = round(as.numeric(popbio::stable.stage(pop_mat)[3]),3),
                      ss_adultfemale = round(as.numeric(popbio::stable.stage(pop_mat)[4]),3),
                      egghatch = eggs,
                      j_early = eju,
                      j_late = lju,
                      fem = ad.fe
  )
  fem0 <- dfout[(dfout$ayear == 0), 'fem']
  dft <- data.frame(dfout, fem_t0 = fem0)
  dft$fem_diff <- round(((dft$fem - dft$fem_t0) / dft$fem_t0), 3)
  dft$change50_flag <- as.integer(ifelse(abs(dft$fem_diff) >= 0.5, 1, 0))
  dft$change30_flag <- as.integer(ifelse(abs(dft$fem_diff) >= 0.3, 1, 0))
  dft$double_flag <- as.integer(ifelse(dft$fem_diff >= 1.0, 1, 0))
 return(dft)
}
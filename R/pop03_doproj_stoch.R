#' @title Create data.frame with projection results.
#' 
#' @description Creates data.frame with projection results for use in subsequent
#' population projection function ("PopProj.R"). 
#'
#' @param x data.frame created by "PopPrep.R".
#'
#' @details This function applies stochastic population projection moodels. 
#' 
#' @return Creates data.frame with projection results.
#' @export
#'
#' @examples
#' \dontrun{
#' dfproj_vals <- pop03_doproj_stoch()
#' }
#' 

pop03_doproj_stoch <- function(x) {
  pop_n <-  x$adultF_n * c(11.1, 4, 2, 1) 
  stage_names <- c("a1", "a2", "a3", "a4",
                   "b1", "b2", "b3", "b4",
                   "c1", "c2", "c3", "c4",
                   "d1", "d2", "d3", "d4")
  stoch_stage_names <- "b1"
  year_val <- x[ , stoch_stage_names]
  bad_year_thresh <- 0.1
  good_year_thresh <- 0.2
 din <-  x[1 , stage_names] 
 # make vector of stochastic values
 vb1 <- c(0.0001, 0.1, year_val, (year_val + (0.1 * year_val)), 
                                  (year_val + (0.2 * year_val)))
 vb1[vb1 > 1] <- 1
 row_n <- 5
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
 
   # Scenarios
   # Default, the projection selects each matrix with an equal 
   # probability at each time interval.
   # Bad years at 2/5 = 0.4 chance of being chosen
   #set.seed(c(1,2,3,4))
   ttime <- 100
   spr_01 <- popdemo::project(matlist, vector = pop_n, time = ttime)
   # Bad and good years have equal probability
   p2 <- 0.5
   m2 <- matrix(rep(c(p2/2, p2/2, (1-p2)/3, (1-p2)/3, (1-p2)/3), 5), 5, 5)
   spr_02 <- popdemo::project(matlist, vector = pop_n, Aseq = m2, time = ttime)
   # Bad years occur twice as often as good years.
   p3 <- 0.667
   m3 <- matrix(rep(c(p3/2, p3/2, (1-p3)/3, (1-p3)/3, (1-p3)/3), 5), 5, 5)
   spr_03 <- popdemo::project(matlist, vector= pop_n, Aseq = m3, time = ttime)
   # Bad years now occur four times as often as good years.
   p4 <- 0.8
   m4 <- matrix(rep(c(p4/2, p4/2, (1-p4)/3, (1-p4)/3, (1-p4)/3), 5), 5, 5)
   spr_04 <- popdemo::project(matlist, vector= pop_n, Aseq = m4, time = ttime)
   
   dfpop <- rbind(data.frame(model = "Stochastic uniform", ayear = 0:ttime, 
                       egghatch = popdemo::vec(spr_01)[,1], 
                       j_early = popdemo::vec(spr_01)[,2], 
                       j_late = popdemo::vec(spr_01)[,3], 
                       fem = popdemo::vec(spr_01)[,4], 
                pop = spr_01), 
                data.frame(model = "Stochastic equal", ayear = 0:ttime, 
                           egghatch = popdemo::vec(spr_02)[,1], 
                           j_early = popdemo::vec(spr_02)[,2], 
                           j_late = popdemo::vec(spr_02)[,3], 
                           fem = popdemo::vec(spr_02)[,4], 
                           pop = spr_02),
                data.frame(model = "Stochastic bad x2", ayear = 0:ttime, 
                           egghatch = popdemo::vec(spr_03)[,1], 
                           j_early = popdemo::vec(spr_03)[,2], 
                           j_late = popdemo::vec(spr_03)[,3], 
                           fem = popdemo::vec(spr_03)[,4], 
                           pop = spr_03),
                data.frame(model = "Stochastic bad x4", ayear = 0:ttime, 
                           egghatch = popdemo::vec(spr_04)[,1], 
                           j_early = popdemo::vec(spr_04)[,2], 
                           j_late = popdemo::vec(spr_04)[,3], 
                           fem = popdemo::vec(spr_04)[,4], 
                           pop = spr_04) 
   )
                
   return(dfpop)
}
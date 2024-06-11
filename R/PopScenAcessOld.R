#' Title
#' @title Extract and combine data for different scenarios.
#' @description Obtains results from "PopProj.R" and combines to 
#' enable subsequent processing. Adapted to evaluate sensitivity to
#' differences in acessibility.
#' Ensures population values do not inflate beyond initial values (max 100 *).
#' 
#' @param x Lookup dataframe created by "PopProj.R".
#' 
#' @details Obtains results from "PopProj.R" and combines to 
#' enable subsequent processing. 
#' Ensures population values do not inflate beyond initial values (max 100 *).
#'
#' @return List holding three dataframes, with results for each scenario.
#' 
#' @importFrom magrittr %>%
#' @importFrom utils read.csv
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' #1) River length coverage
#' lsf <- prepTabcover(pBasin = B, pBasinSp = Bsp, 
#' pBasinC = BC, riv = rin, rastAc = ras1, make_shape = FALSE)
#' rp <- system.file("shape/shapes_rivers3395", package="cmartr")
#' lt <- resTab(listsf = lsf, input_rp = rp, make_html = FALSE)
#' atest <- lt$rlcb
#' # add key and make sure all levels are represented
#' atest$namekey <- paste(atest$BASIN_N, atest$name, atest$subbasn, sep = "_")
#' riverl <- expand.grid(namekey = unique(atest$namekey), 
#'                      accessible = unique(atest$accessible) )
#' riverl <- merge(riverl, atest, all.x=TRUE)
#' selNA <- which(is.na(riverl$tot_km))
#' riverl[selNA, c("tot_km", "tot_notPA", "tot_PA", "tot_Ind", "tot_SP", "tot_use")] <- 0
#' 
#' #2) Data frame with population parameters created from "PopParam.R" 
#' dfpop <- readRDS("inst/other/dfpop.RDS") 
#' 
#' #3) Create list with population parameters for different scenarios
#' # across river lengths per geographic/political coverage class: 
#' # basin, country, subbasin, accessible, protected area...
#' l.gpop <- plyr::dlply(dfpop, c("akey"), PopPrep, riverl=riverl)
#' 
#' #4) Project scenarios across species range
#' # takes 5 hours and writes 10 GB of results.
#' # Do not run unless you really want to.....
#' dflup <- plyr::ldply(l.gpop, PopProj, write_csv = TRUE, 
#' write_db = FALSE)
#' 
#' #5) Get results from 3 scenarios
#' lscen <- PopScen(dflup)
#' }
PopScenAcess <- function(x){
  # x = dflup from "projPop.R"
  # 2-3 mins, time takes to load large .csv files
  # selects relevant projection data from modelled scenarios
  # x <- dflup
  f.base <- which(x$hunt==0 & x$increase==0.2)
  base <- read.csv(as.character(x[f.base,'fileout']))
  f.hunt <- which(x$hunt==10 & x$increase==0.1)
  hunt <- read.csv(as.character(x[f.hunt,'fileout']))
  f.manage <- which(x$hunt==10 & x$increase==0.5)
  manage <- read.csv(as.character(x[f.manage,'fileout']))
  
  # Business as usual, final year totals
  # Accessible populations with hunt
  selBA <- which(hunt$ayear == max(hunt$ayear) & 
                   hunt$accessible == "Yes" &
                   hunt$variable %in% c("tot_notPA", "tot_PA")# "tot_km" #&
                   #hunt$prop_km == 1
  )
  # Inaccessible at base rates.
  selBNA <- which(base$ayear == max(base$ayear) & 
                    base$accessible == "No" &
                    base$variable %in% c("tot_notPA", "tot_PA") #&
                    #base$prop_km == 1
  )
  
  h1 <- hunt[selBA, ]
  h1$rank_prop <- h1$prop_km
  h2 <- base[selBNA, ]
  h2$rank_prop <- (1 - h2$prop_km)
  
  dft.BAU <-merge(h1, h2, 
                  by = c("namekey", "variable", "rank_prop"))
  dft.BAU$namekeys <- dft.BAU$namekey
  #library(magrittr)
  dft.BAU <- dft.BAU %>% tidyr::separate(namekeys, 
                                  into = c("basin", "country", "subbasin"), sep="_")
  # limit population to 100 times original (i.e. 10*100 = 1000 adult females per km)
  if(dft.BAU$lambda.x[1] > 0.99999){
    dft.BAU$adult_females.x <- ifelse(dft.BAU$adult_females.x > (dft.BAU$fem_t0.x*100), 
                                      dft.BAU$fem_t0.x*100, 
                                      dft.BAU$adult_females.x)
  }
  
  if(dft.BAU$lambda.y[1] > 0.99999){
    dft.BAU$adult_females.y <- ifelse(dft.BAU$adult_females.y > (dft.BAU$fem_t0.y*100), 
                                      dft.BAU$fem_t0.y*100, 
                                      dft.BAU$adult_females.y)
  }
  
  dft.BAU$tot_fem <- dft.BAU$adult_females.x + dft.BAU$adult_females.y
  dft.BAU$tot_fem_t0 <- dft.BAU$fem_t0.x + dft.BAU$fem_t0.y
  dft.BAU$adult_female_diff <- round(((dft.BAU$tot_fem  - dft.BAU$tot_fem_t0) / dft.BAU$tot_fem_t0), 3)
  dft.BAU$change50_flag <- as.integer(ifelse(abs(dft.BAU$adult_female_diff) > 0.499, 1, 0))
  dft.BAU$double_flag <- as.integer(ifelse(dft.BAU$adult_female_diff > 0.999, 1, 0))
  
  # Strict protection, 
  #Accessible with PAs set to base
  selSPA.pa <- which(base$ayear == max(base$ayear) & 
                       base$accessible == "Yes" &
                       base$variable == "tot_PA" #&
                       #base$prop_km == 1
  ) # 80 rows 
  #Accessible not PAs set to nest collection (first year survival 0.1) and
  # adult harvest (10%).
  selSPA.npa <- which(hunt$ayear == max(hunt$ayear) & 
                        hunt$accessible == "Yes" &
                        hunt$variable == "tot_notPA" #&
                        #hunt$prop_km == 1 
                        )
  a1.sp <- base[selSPA.pa, ]
  a2.sp <- hunt[selSPA.npa, ]
  spa <- rbind(a1.sp, a2.sp)
  spa$rank_prop <- spa$prop_km
  
  # Inaccessible at base rates.
  selSPNA <- which(base$ayear == max(base$ayear) & 
                     base$accessible == "No" &
                     base$variable %in% c("tot_notPA", "tot_PA")
                     #base$prop_km == 1 
                     )
  spia <- base[selSPNA, ]
  spia$rank_prop <- (1 - spia$prop_km)
  
  dft.SP <- merge(spa, spia,
                  by = c("namekey", "variable", "rank_prop"))
  dft.SP$namekeys <- dft.SP$namekey
  dft.SP <- dft.SP %>% tidyr::separate(namekeys, 
                                  into = c("basin", "country", "subbasin"), sep="_")
  # limit population to 100 times original (i.e. 10*100 = 1000 adult females per km)
  if(dft.SP$lambda.x[1] > 0.99999){
    dft.SP$adult_females.x <- ifelse(dft.SP$adult_females.x > (dft.SP$fem_t0.x*100), 
                                     dft.SP$fem_t0.x*100, 
                                     dft.SP$adult_females.x)
  }
  if(dft.SP$lambda.y[1] > 0.99999){
    dft.SP$adult_females.y <- ifelse(dft.SP$adult_females.y > (dft.SP$fem_t0.y*100), 
                                     dft.SP$fem_t0.y*100, 
                                     dft.SP$adult_females.y)
  }
 
  dft.SP$tot_fem <- dft.SP$adult_females.x + dft.SP$adult_females.y
  dft.SP$tot_fem_t0 <- dft.SP$fem_t0.x + dft.SP$fem_t0.y
  dft.SP$adult_female_diff <- round(((dft.SP$tot_fem  - dft.SP$tot_fem_t0) / dft.SP$tot_fem_t0), 3)
  dft.SP$change50_flag <- as.integer(ifelse(abs(dft.SP$adult_female_diff) > 0.499, 1, 0))
  dft.SP$double_flag <- as.integer(ifelse(dft.SP$adult_female_diff > 0.999, 1, 0))
  
  # Community management
  #Accessible with PAs set to hunt
  selCMA.pa <- which(hunt$ayear == max(hunt$ayear) & 
                       hunt$accessible == "Yes" &
                       hunt$variable == "tot_PA" #&
                       #hunt$prop_km == 1 
                       ) 
  #Accessible not PAs to manage
  selCMA.npa <- which(manage$ayear == max(manage$ayear) & 
                        manage$accessible == "Yes" &
                        manage$variable == "tot_notPA" #&
                        #manage$prop_km == 1 
                        ) 
  a1.cm <- hunt[selCMA.pa, ]
  a2.cm <- manage[selCMA.npa, ]
  cma <- rbind(a1.cm, a2.cm)
  cma$rank_prop <- cma$prop_km
  
  # Inaccessible at base rates.
  selCMNA <- which(base$ayear == max(base$ayear) & 
                     base$accessible == "No" &
                     base$variable %in% c("tot_notPA", "tot_PA")
                     #base$prop_km == 1 
                     ) 
  cmi <- base[selCMNA, ]
  cmi$rank_prop <- (1 - cmi$prop_km)
  
  dft.CM <- merge(cma, cmi,
                  by = c("namekey", "variable", "rank_prop"))

  dft.CM$namekeys <- dft.CM$namekey
  dft.CM <- dft.CM %>% tidyr::separate(namekeys, 
                                into = c("basin", "country", "subbasin"), sep="_")
  # limit population to 100 times original (i.e. 10*100 = 1000 adult females per km)
  if(dft.CM$lambda.x[1] > 0.99999){
    dft.CM$adult_females.x <- ifelse(dft.CM$adult_females.x > (dft.CM$fem_t0.x*100), 
                                     dft.CM$fem_t0.x*100, 
                                     dft.CM$adult_females.x)
  }
  if(dft.CM$lambda.y[1] > 0.99999){
    dft.CM$adult_females.y <- ifelse(dft.CM$adult_females.y > (dft.CM$fem_t0.y*100), 
                                     dft.CM$fem_t0.y*100, 
                                     dft.CM$adult_females.y)
  }

  dft.CM$tot_fem <- dft.CM$adult_females.x + dft.CM$adult_females.y
  dft.CM$tot_fem_t0 <- dft.CM$fem_t0.x + dft.CM$fem_t0.y
  dft.CM$adult_female_diff <- round(((dft.CM$tot_fem  - dft.CM$tot_fem_t0) / dft.CM$tot_fem_t0), 3)
  dft.CM$change50_flag <- as.integer(ifelse(abs(dft.CM$adult_female_diff) > 0.499, 1, 0))
  dft.CM$double_flag <- as.integer(ifelse(dft.CM$adult_female_diff > 0.999, 1, 0))
  
 l1 <-  list(dft.BAU = dft.BAU, dft.SP = dft.SP, dft.CM = dft.CM)
  return(l1)
}
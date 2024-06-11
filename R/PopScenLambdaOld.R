#' Title
#' @title Extract and combine data for different scenarios.
#' @description Obtains results from "PopProj.R" and combines to 
#' enable subsequent processing. Adapted to evaluate lambda for different values.
#' Ensures population values do not inflate beyond initial values (max 100 *).
#' 
#' @param x Lookup dataframe created by "PopProj.R".
#' 
#' @details Obtains results from "PopProj.R" and returns min and max lambda . 
#'
#' @return Data.frame, with Lambda for each projection run.
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
#' #5) Get results
#' library(plyr)
#' dfl <- ddply(dflup, .(akey), .fun = PopScenLambda) 
#' min(dfl$lambda_min); max(dfl$lambda_max) 
#' [1] 0.465 [1] 1.163518
#' }
PopScenLambda <- function(x){
  # x = dflup from "projPop.R"
  # 2-3 mins, time takes to load large .csv files
  # selects relevant projection data from modelled scenarios
  # x <- dflup
  base <- read.csv(as.character(x[,'fileout'])) 
  minl <- min(base$lambda) 
  maxl <- max(base$lambda) 
  dfout <- data.frame(lambda_min = minl, lambda_max = maxl)

  return(dfout)
}
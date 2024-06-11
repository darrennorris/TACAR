#' @title Prepare data for population projections.
#' 
#' @description Creates list holding parameter matrix and river lengths. 
#'
#' @param x Data.frame with population parameters. Created by "PopParam.R".
#' @param riverl Data.frame with river lengths per category of interest. 
#' Created by "resTab.R".
#'
#' @return List with population parameter matrix and river lengths.
#' @export
#' 
#' @examples
#' \dontrun{
#' #1) River length coverage
#' # Run, 4 - 5 hours depending on computer and memory allocation
#' lsf <- prepTabcover(pBasin = B, pBasinSp = Bsp, 
#'           pBasinC = BC, riv = rin, rastAc = ras1, make_shape = FALSE)
#' 
#' # Folder with shapefiles with 1km rvier sections
#' rp <- system.file("shape/shapes_rivers3395", package="cmartr")
#' lt <- resTab(listsf = lsf, input_rp = rp, make_html = FALSE)
#' atest <- lt$rlcb
#' # add key and make sure all levels present
#' atest$namekey <- paste(atest$BASIN_N, atest$name, atest$subbasn, sep = "_")
#' riverl <- expand.grid(namekey = unique(atest$namekey), 
#'                      accessible = unique(atest$accessible) )
#' riverl <- merge(riverl, atest, all.x=TRUE)
#' selNA <- which(is.na(riverl$tot_km))
#' coln <- c("tot_km", "tot_notPA", "tot_PA", "tot_Ind", "tot_SP", "tot_use")
#' riverl[selNA, coln] <- 0
#' 
#' #2) Data.frame with population parameters
#' dfpop <- PopParam(species = "Podocnemis unifilis", make_rds = FALSE)
#' 
#' #3) Create list with population parameters for different scenarios
#' # Run. Few seconds 33.5 Mb
#' # Create list with population parameters for different scenarios
#' # across river lengths per geographic/political coverage class: 
#' # basin, country, subbasin, accessible, protected area...
#' l.gpop <- plyr::dlply(dfpop, c("akey"), PopPrep, riverl=riverl)
#' }

PopPrep <- function(x = NA, riverl = NA){
  # function to prepare data for projection
  vpop <- unlist(x[ ,4:19])
  tracaja <- matrix(vpop, byrow = TRUE, ncol=4)
  dimnames(tracaja) <- list(c("a", "b", "c", "d"),
                            c(1,2,3,4))
  
  #numeric vector of individuals at different age stage
  dft <- reshape2::melt(riverl, id.vars = c('namekey', 'accessible'), 
                        measure.vars = c('tot_km', 'tot_notPA', 'tot_PA',
                                         'tot_Ind', 'tot_SP', 'tot_use'), 
                        value.name = c('distKMa'))
  
  myq <- function(x){
    myquant <- seq(0.05,1,by=0.05)
    dfout <- data.frame(prop_km = myquant, dist_km = x$distKMa * myquant)
    dfout
  }
  
  dft2 <- plyr::ddply(dft, c("namekey", "accessible", "variable"), myq)
  adultF.d <- 10 # adult female density per river km
  dft2$adultF.n <- trunc(adultF.d * dft2$dist_km)
  dft2 <- data.frame(species = x$species, hunt = x$type, 
                     increase = x$increase, dft2)
  l1 <- list(rdata = dft2, tracajam = tracaja)
  
}
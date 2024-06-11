#' @title Create data.frame with demographic parameters
#' 
#' @description Creates data.frame with parameters for use in subsequent
#' population projection function ("PopProj.R") 
#'
#' @param species Name of species. Default is "Podocnemis unifilis"
#' @param make_rds Logical (TRUE/FALSE). Should .RDS format be written.
#' Default is FALSE
#'
#' @details This function holds the demographic parameters used 
#' in population projections. When called by a user, it will create 
#' a data.frame with these demogrpahic parameters.
#' @return Creates data.frame with parameters used in "PopProj.R".
#' @export
#'
#' @examples
#' \dontrun{
#' dfpop <- PopParam(species = "Podocnemis unifilis", make_rds = FALSE)
#' }
#' 
pop01_param_poun <- function(species = "Podocnemis unifilis", make_rds = FALSE){
  
  sp <- species
  # parameters
  # female fecundity
  cfreq = 1.10 #clutch frequency
  csize = 21.00 # clutch size James
  csize.araguari = 14.00 # clutch size Araguari
  sratio = 0.5 # sex ratio
  ffec <- (cfreq * csize * sratio) # female fecundity
  
  #"Graduation" (First year survival and transition)
  esur <- 0.10
  emonths <- 2
  hsur <- 0.22
  hmonths <- (12 - emonths)
  ehsur <- ((esur*emonths)+(hsur*hmonths))/(emonths+hmonths) #0.2
  
  # define first year Graduation for projections
  # use 0.0001 so matrix remains irreducible
  ehsur1 <- c(0.0001, 0.1, ehsur, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  ehsur_text <- c(0, 0.1, ehsur, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
  
  # juvenile
  jearlysur <- 0.5
  jearlyears <- 2
  jlatesur <- 0.4
  jlateyears <- 2
  
  # early juvenile survive and stay (P) / survive and transition (G)
  jePn <- 1 - (jearlysur^(jearlyears - 1))
  jePd <- (1 - (jearlysur^jearlyears))
  jeP <- (jePn / jePd)*jearlysur
  jeG <- ((jearlysur^jearlyears)*(1-jearlysur)) / (1 - (jearlysur^jearlyears))
  
  # late juvenile survive and stay (P) / survive and transition (G)
  jlPn <- 1 - (jlatesur^(jlateyears - 1))
  jlPd <- (1 - (jlatesur^jlateyears))
  jlP <- (jlPn / jlPd)*jlatesur
  jlG <- ((jlatesur^jlateyears)*(1-jlatesur)) / (1 - (jlatesur^jlateyears))
  
  #adult survival
  asur <- 0.93
  
  
  # make data frame with row holding parameters for each population 
  dfpop <- 
    rbind(
      data.frame(
        species = sp,
        type = "early-protection",   #set PPM element values,
        increase = as.character(ehsur_text),
        a1 = 0,   #leave = 0
        a2 = 0,  #leave = 0,
        a3 = 0,
        a4 = (ffec * asur),  #fecundity per female,
        b1 = ehsur1, 
        b2 = 0.333333,  #juv_early survival adjusted for years as juv_early,
        b3 = 0,  # leave = zero,
        b4 = 0,
        c1 = 0,
        c2 = 0.166667,
        c3 = 0.285714, #juv_early ~graduation~ probability
        c4 = 0,
        d1 = 0,  #leave = zero,
        d2 = 0, #leave = zero
        d3 = 0.114286, #juv_late "graduation" probability,
        d4 = asur  #default = 0.93
      ),
      data.frame(
        species = sp,
        type = "early-protection, female-hunt 2.5%",   #set PPM element values
        increase = as.character(ehsur_text),
        a1 = 0,   #leave = 0
        a2 = 0,  #leave = 0,
        a3 = 0,
        a4 = (ffec * (asur * 0.975)),  #fecundity per female,
        b1 = ehsur1, 
        b2 = 0.333333,  #juv_early survival adjusted for years as juv_early,
        b3 = 0,  # leave = zero,
        b4 = 0,
        c1 = 0,
        c2 = 0.166667,
        c3 = 0.285714, #juv_early ~graduation~ probability
        c4 = 0,
        d1 = 0,  #leave = zero,
        d2 = 0, #leave = zero
        d3 = 0.114286, #juv_late "graduation" probability,
        d4 = (asur * 0.975) #default = 0.93
      ),
      data.frame(
        species = sp,
        type = "early-protection, female-hunt 10%",   #set PPM element values
        increase = as.character(ehsur_text),
        a1 = 0,   #leave = 0
        a2 = 0,  #leave = 0,
        a3 = 0,
        a4 = (ffec * (asur * 0.90)),  #fecundity per female,
        b1 = ehsur1, 
        b2 = 0.333333,  #juv_early survival adjusted for years as juv_early,
        b3 = 0,  # leave = zero,
        b4 = 0,
        c1 = 0,
        c2 = 0.166667,
        c3 = 0.285714, #juv_early ~graduation~ probability
        c4 = 0,
        d1 = 0,  #leave = zero,
        d2 = 0, #leave = zero
        d3 = 0.114286, #juv_late "graduation" probability,
        d4 = (asur * 0.90) #default = 0.93
      ),
      data.frame(
        species = sp,
        type = "early-protection, female-hunt 25%",   #set PPM element values
        increase = as.character(ehsur_text),
        a1 = 0,   #leave = 0
        a2 = 0,  #leave = 0,
        a3 = 0,
        a4 = (ffec * (asur * 0.75)),  #fecundity per female,
        b1 = ehsur1, 
        b2 = 0.333333,  #juv_early survival adjusted for years as juv_early,
        b3 = 0,  # leave = zero,
        b4 = 0,
        c1 = 0,
        c2 = 0.166667,
        c3 = 0.285714, #juv_early ~graduation~ probability
        c4 = 0,
        d1 = 0,  #leave = zero,
        d2 = 0, #leave = zero
        d3 = 0.114286, #juv_late "graduation" probability,
        d4 = (asur * 0.75)  #default = 0.93
      ),
      data.frame(
        species = sp,
        type = "early-protection, female-hunt 50%",   #set PPM element values
        increase = as.character(ehsur_text),
        a1 = 0,   #leave = 0
        a2 = 0,  #leave = 0,
        a3 = 0,
        a4 = (ffec * (asur * 0.5)),  #fecundity per female,
        b1 = ehsur1, 
        b2 = 0.333333,  #juv_early survival adjusted for years as juv_early,
        b3 = 0,  # leave = zero,
        b4 = 0,
        c1 = 0,
        c2 = 0.166667,
        c3 = 0.285714, #juv_early ~graduation~ probability
        c4 = 0,
        d1 = 0,  #leave = zero,
        d2 = 0, #leave = zero
        d3 = 0.114286, #juv_late "graduation" probability,
        d4 = (asur * 0.5)  #default = 0.93
      )
    )
  
  dfpop$akey <- paste(dfpop$species, dfpop$type, dfpop$increase, sep = "_")
  levels(dfpop$type) <- c("0", "2.5", "10", "25", "50") 
  dfpop$type <- as.character(dfpop$type)
  dfpop$increase <- as.numeric(as.character(dfpop$increase))
  
  if(make_rds!=FALSE){
    saveRDS(dfpop, "inst/other/dfpop_poun.RDS") 
  }
  return(dfpop)
  
}
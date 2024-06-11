#' Title
#'
#' @title Project populations.
#' 
#' @param x List with population parameters for different scenarios.
#' Creates using "PopPrep.R"
#' @param write_csv Logical (TRUE/FALSE). Should csv tables be written.
#' @param write_db Logical (TRUE/FALSE). Should postgre tables be written.
#' 
#' @description Projects populations through different scenarios across 
#' species range.
#'
#' @return Results as .csv and/or postgresql tables. A lookup table is returned for users that 
#' is then used for producing reaults.
#' @importFrom utils write.csv
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
#' coln <- c("tot_km", "tot_notPA", "tot_PA", "tot_Ind", "tot_SP", "tot_use")
#' riverl[selNA, coln] <- 0
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
#' }
PopProj <- function(x, write_csv = FALSE, write_db = FALSE){
  tracaja <- x$tracajam
  
  doproj <- function(x) {
    tracaja_n <-  x$adultF.n * c(11.1, 4, 2, 1) 
    
    # project PPM 
    pr_tracaja <- popdemo::project(tracaja, vector=tracaja_n, time=50)
    
    # data for plotting
    len <- length(pr_tracaja)
    Time.intervals <- 0:(len - 1)
    eggs <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[1])))
    eju <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[2])))
    lju <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[3])))
    ad.fe <- as.integer(trunc(pr_tracaja * (popbio::stable.stage(tracaja)[4])))
    plambda = popbio::lambda(tracaja)
    
    # make dataframe 
    dfout <- data.frame(lambda = plambda,
                        ayear = Time.intervals, 
                        individuals = as.integer(trunc(pr_tracaja)),
                        ss_egghatchling = round(as.numeric(popbio::stable.stage(tracaja)[1]),3),
                        ss_earlyjuven = round(as.numeric(popbio::stable.stage(tracaja)[2]),3),
                        ss_latejuven = round(as.numeric(popbio::stable.stage(tracaja)[3]),3),
                        ss_adultfemale = round(as.numeric(popbio::stable.stage(tracaja)[4]),3),
                        egghatch = eggs,
                        early_juven = eju,
                        late_juven = lju,
                        adult_females = ad.fe
    )
    fem0 <- dfout[(dfout$ayear == 0), 'adult_females']
    dft <- data.frame(dfout, fem_t0 = fem0)
    dft$adult_female_diff <- round(((dft$adult_females - dft$fem_t0) / dft$fem_t0), 3)
    dft$change50_flag <- as.integer(ifelse(abs(dft$adult_female_diff) > 0.499, 1, 0))
    dft$double_flag <- as.integer(ifelse(dft$adult_female_diff > 0.999, 1, 0))
    dft
  }
  
  dfin <- x$rdata
  #dfin <- l.gpop$`Podocnemis unifilis.headstart.0`$rdata
  dout <- plyr::ddply(dfin, c("species", "hunt", "increase","namekey", "accessible",
                              "variable", "prop_km", "dist_km"), doproj)
  # multiply by 10 to avoid decimal points in file names
  fn <- paste(dfin[1,'species'], dfin[1,'hunt'] * 10, dfin[1,'increase'] * 10, sep="_")
  fname <- paste(fn,"csv", sep = ".")
  csvout <- paste("inst/other/dataproj/",fname, sep = "")
  
  if(write_csv != FALSE){
    write.csv(dout, csvout, row.names = FALSE)
  }
  
  if(write_db != FALSE){
    drv <- RPostgreSQL::PostgreSQL() 
    db <- 'postgres'  
    host_db <- 'localhost'  
    db_port <- '5432'  
    db_user <- 'postgres'  
    db_password <- 'bob1975'
    
    conn <- RPostgreSQL::dbConnect(drv, dbname=db, host=host_db, 
                      port=db_port, user=db_user, 
                      password=db_password)

    if(RPostgreSQL::dbExistsTable(conn, c("turtles", "test"))){
      RPostgreSQL::dbWriteTable(conn, c("turtles", "test"), dout, 
                   row.names = FALSE, 
                   append = T)
    }else{
      RPostgreSQL::dbWriteTable(conn, c("turtles", "test"), dout, 
                   row.names = FALSE)
      sql_command <- "ALTER TABLE turtles.test SET UNLOGGED;"
      RPostgreSQL::dbGetQuery(conn, sql_command)
    }
    RPostgreSQL::dbDisconnect(conn) 
  }
  rm("dout")
  dflup <- data.frame(species = dfin[1,'species'], hunt = dfin[1,'hunt'], 
                      increase = dfin[1,'increase'],
             fileout = csvout, csvname = fname)
  dflup
}
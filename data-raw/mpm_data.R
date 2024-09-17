library(Rcompadre)
# Below adapted from cdb_fetch. Vignette still wont build with load.
env <- new.env()
dbfile <- system.file("other\\COMADRE_v.4.23.3.1.RData", package = "TACAR")
x <- load(dbfile, env)[1]
dbFetch <- env[[x]]
# convert to CompadreDB
if (inherits(dbFetch, "CompadreDB")) {
  comadre <- dbFetch
} else {
  comadre <- Rcompadre::as_cdb(dbFetch)
}


# data for turtles and tortoises
# A COM(P)ADRE database ('CompadreDB') object with 18 SPECIES and 89 MATRICES.
comadre_test <- subset(comadre, Order == "Testudines")
# A COM(P)ADRE database ('CompadreDB') object with 2 SPECIES and 2 MATRICES.
comadre_podoc <- subset(comadre_test, Family == "Podocnemididae")
usethis::use_data(comadre_podoc)

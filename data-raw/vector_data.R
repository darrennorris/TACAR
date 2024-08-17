## code to prepare vector data goes here
library(dplyr)
library(sf)
library(usethis)

# .rds with Free-flowing river points made by TACAR https://github.com/darrennorris/TACAR
# made by test_datajoin.Rmd
mypath <- "yourpathhere"
mypath <- "C:\\Users\\user\\Documents\\Articles\\2024_Norris_Greenstatus\\TACAR\\inst\\other\\scenario_res_ffr1a5.rds"
# points
inpoints <- readRDS(mypath)
# apply cieling threshold
nf <- 10
ceiling_threshold <- nf + (nf * 0.2)
# keep only columns needed for plotting
points_bau_ffr <- inpoints |>
  filter(model_name == "modelkey_BAU") |>
  mutate(flag_50_35y = factor(if_else(fem_diff_t35 <= -0.5, 1, 0))) |>
  select(BASIN_NAME, subbasin, SUBBASIN_FLAG,
         BAS_NAME, COUNTRY, RIV_ORD, BB_ID, BB_NAME, REACH_ID,
         Protected, Protected_cat, Accessible, Free_flowing,
         fem_t0, fem_t35, fem_t100, fem_diff_t35, flag_50_35y,
         geom) |>
  mutate(fem_t35 = ifelse(fem_t35 > ceiling_threshold, ceiling_threshold, fem_t35),
         fem_t100 = ifelse(fem_t100 > ceiling_threshold, ceiling_threshold, fem_t100)) |>
  sf::st_as_sf()
sf::st_crs(points_bau_ffr) <- NA
points_bau_ffr <- points_bau_ffr |> data.frame()
# export
usethis::use_data(points_bau_ffr, overwrite = TRUE)

# subset for mapping
points_bau_ffr_map <- points_bau_ffr |>
  arrange(BASIN_NAME, subbasin, BB_ID, REACH_ID) |>
  filter(row_number() %% 10 == 1)
# export
usethis::use_data(points_bau_ffr_map, overwrite = TRUE)

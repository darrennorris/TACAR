#' Points along Free-flowing rivers.
#'
#' A lot of points along rivers.
#'
#' @description
#'
#' The points show projected population changes along rivers.
#' Projections generated in TACAR (https://github.com/darrennorris/TACAR).
#' The population projections are an updated and extended version of those
#' published by Norris et. al. 2019 (https://doi.org/10.1016/j.biocon.2019.02.022).
#' \describe{
#'   \item{location}{South America.}
#'   \item{coord. ref}{WGS 84 (EPSG:3395).}
#' }
#'
#' @format `points_bau_ffr`
#' is a dataframe with 353437 points and 15 fields, including:
#' \describe{
#'   \item{BASIN_NAME}{Name of major basin.}
#'   \item{subbasn}{Name of subbasn.}
#'   \item{geom}{sf geometry column.}
#' }
#'
#' @references
#' Norris, D., et. al. (2019).
#' Prospects for freshwater turtle population recovery are catalyzed
#' by pan-Amazonian community-based management.
#' Biological Conservation, 233, https://doi.org/10.1016/j.biocon.2019.02.022.
#'
#'Grill, G., Lehner, B., Thieme, M. et al. (2019)
#'Mapping the world’s free-flowing rivers.
#'Nature 569, 215–221, https://doi.org/10.1038/s41586-019-1111-9.
#'
#' @source <https://github.com/darrennorris/TACAR>
#' @importFrom sf st_as_sf
#' @examples
#' \dontrun{
#' # convert to sf object
#' sf_points_bau_ffr <- sf::st_as_sf(points_bau_ffr, crs = 3395)
#' }
"points_bau_ffr"
NULL

#' Subset of points along Free-flowing rivers.
#'
#' A 10% subset, includes one in ten rows of points_bau_ffr.
#'
#' @description
#'
#' The points show projected population changes along rivers.
#' Projections generated in TACAR (https://github.com/darrennorris/TACAR).
#' The population projections are an updated and extended version of those
#' published by Norris et. al. 2019 (https://doi.org/10.1016/j.biocon.2019.02.022).
#' \describe{
#'   \item{location}{South America.}
#'   \item{coord. ref}{WGS 84 (EPSG:3395).}
#' }
#'
#' @format `points_bau_ffr_map`
#' is a dataframe with 35344 points and 15 fields, including:
#' \describe{
#'   \item{BASIN_NAME}{Name of major basin.}
#'   \item{subbasn}{Name of subbasn.}
#'   \item{geom}{sf geometry column.}
#' }
#'
#' @references
#' Norris, D., et. al. (2019).
#' Prospects for freshwater turtle population recovery are catalyzed
#' by pan-Amazonian community-based management.
#' Biological Conservation, 233, https://doi.org/10.1016/j.biocon.2019.02.022.
#'
#'Grill, G., Lehner, B., Thieme, M. et al. (2019)
#'Mapping the world’s free-flowing rivers.
#'Nature 569, 215–221, https://doi.org/10.1038/s41586-019-1111-9.
#'
#' @source <https://github.com/darrennorris/TACAR>
#' @importFrom sf st_as_sf
#' @examples
#' \dontrun{
#' # convert to sf object
#' sf_points_bau_ffr_map <- sf::st_as_sf(points_bau_ffr_map, crs = 3395)
#' }
"points_bau_ffr_map"
NULL

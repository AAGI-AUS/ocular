# =========================================================================
# Core: geometry helpers
# =========================================================================

#' Convert metres to decimal degrees
#'
#' Converts a ground distance to its east-west and north-south angular spans at
#' a given latitude using the WGS-84 ellipsoid. The ellipsoid radii are
#' Simpson-averaged across a latitude band of width \code{distance}; this band
#' controls the averaging interval, not the distance being converted.
#'
#' @param latitude Centre latitude in decimal degrees.
#' @param distance Full north-south width, in metres, of the latitude band used
#'   for averaging.
#' @param pixel_size Ground distance in metres to convert.
#' @returns A named numeric vector with \code{x_degrees} along the parallel and
#'   \code{y_degrees} along the meridian.
#' @examples
#' metres_to_degrees(latitude = -32, distance = 1000, pixel_size = 10)
#' @export
metres_to_degrees <- function(latitude, distance = 1, pixel_size = 1){

  if( !.isScalarNumber(latitude) || latitude <= -90 || latitude >= 90 )
    stop("latitude must be a finite numeric strictly between -90 and 90.",
         call. = FALSE)
  if( !.isScalarNumber(distance) || distance < 0 )
    stop("distance must be a finite non-negative numeric.", call. = FALSE)
  if( !.isScalarNumber(pixel_size) || pixel_size <= 0 )
    stop("pixel_size must be a finite positive numeric.", call. = FALSE)

  a  <- 6378137.0
  e2 <- 0.00669437999014
  lat_c_rad <- latitude * pi / 180
  m_c       <- a * (1 - e2) / (1 - e2 * sin(lat_c_rad)^2)^1.5
  half_lat  <- (distance / 2) * 180 / (pi * m_c)
  if( abs(latitude) + half_lat >= 90 )
    stop("distance makes the averaging band cross a pole.", call. = FALSE)
  lat_rad   <- (latitude + c(-half_lat, 0, half_lat)) * pi / 180
  denom     <- 1 - e2 * sin(lat_rad)^2
  n         <- a / sqrt(denom)
  m         <- a * (1 - e2) / denom^1.5
  w         <- c(1, 4, 1)
  return(c(x_degrees = pixel_size * 180 / (pi * sum(n * cos(lat_rad) * w) / 6),
           y_degrees = pixel_size * 180 / (pi * sum(m * w) / 6)))
}

#' Bounding box from a point
#'
#' Creates a WGS-84 bounding box centred on a point, using the requested full
#' width and height. The calculation requires \code{sf::sf_use_s2(TRUE)} so
#' metre distances are interpreted geodesically.
#'
#' @param longitude,latitude Centre point in WGS-84 decimal degrees.
#' @param x_metres,y_metres Full width and height in metres.
#' @details Requests crossing the antimeridian are rejected because a single
#'   WGS-84 bounding box cannot represent a wrapped interval.
#' @returns An \code{sf::st_bbox} object in WGS-84.
#' @examples
#' if (sf::sf_use_s2()) {
#'   point_to_bbox(117.8, -32, x_metres = 1000, y_metres = 500)
#' }
#' @export
point_to_bbox <- function(longitude, latitude, x_metres, y_metres){

  if( !.isScalarNumber(longitude) || longitude < -180 || longitude > 180 )
    stop("longitude must be a finite numeric in [-180, 180].", call. = FALSE)
  if( !.isScalarNumber(latitude) || latitude <= -90 || latitude >= 90 )
    stop("latitude must be a finite numeric strictly between -90 and 90.",
         call. = FALSE)
  if( !.isScalarNumber(x_metres) || !.isScalarNumber(y_metres) ||
      x_metres <= 0 || y_metres <= 0 )
    stop("x_metres and y_metres must be finite positive numeric scalars.",
         call. = FALSE)
  if( !sf::sf_use_s2() )
    stop("point_to_bbox: requires sf::sf_use_s2(TRUE) so metre distances are ",
         "interpreted geodesically; s2 is currently disabled (the buffer would ",
         "be taken in degrees). Re-enable with sf::sf_use_s2(TRUE).",
         call. = FALSE)
  half_width_deg <- metres_to_degrees(latitude,
                                      distance = y_metres,
                                      pixel_size = x_metres / 2)[["x_degrees"]]
  if( abs(longitude) + half_width_deg >= 180 )
    stop("point_to_bbox: the requested area crosses the antimeridian; ",
         "split it into two WGS-84 bounding boxes.", call. = FALSE)
  point  <- sf::st_sfc(sf::st_point(c(longitude, latitude)), crs = sf::st_crs(4326))
  bbox_x <- sf::st_bbox(sf::st_buffer(point, dist = x_metres / 2))
  bbox_y <- sf::st_bbox(sf::st_buffer(point, dist = y_metres / 2))
  return(sf::st_bbox(c(bbox_x[1], bbox_y[2], bbox_x[3], bbox_y[4]),
                     crs = sf::st_crs(4326)))
}

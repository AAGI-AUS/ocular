#' ocular: Spectral Indices and Sample-Free Field Boundary Delineation
#'
#' ocular retrieves Landsat and Sentinel-2 satellite imagery through the
#' Microsoft Planetary Computer SpatioTemporal Asset Catalog (STAC) and stores
#' the retrieved scenes in an \code{ocular} object. After retrieval, users can
#' configure and apply reusable processing and output pipelines, including
#' sample-free field boundary delineation. The delineation workflow uses a
#' supplied point and a multi-temporal spectral feature stack, minimising the
#' need for labelled training data.
#'
#' \code{\link{add_modis}} is an experimental option that uses daily MCD43A4
#' NBAR data to estimate additional values for vegetative indices (VI) not
#' covered by Landsat.
#' These estimates are not satellite observations. MCD43A4
#' product-specific quality layers are not yet applied. The leave-one-out
#' agreement metrics from \code{\link{validate_data_fusion}} apply to the
#' selected data, place, and period; they do not establish general accuracy or
#' replace independent validation.
#'
#' The main entry points are \code{\link{get_rs}} for retrieval,
#' \code{\link{segment_area}}, \code{\link{trace_perimeter}},
#' \code{\link{segment_interior}}, and \code{\link{split_area}} for custom
#' field boundary delineation, and \code{\link{boundary_delineation}} for the
#' default sequence. \code{\link{as_raster}} and
#' \code{\link{as_time_series}} produce outputs, while
#' \code{\link{batch_rs}} applies retrieval and processing across many sites.
#'
#' \code{\link{add_ftw_prior}} can attach compatible field boundary GeoParquet
#' data, including Fields of The World (FTW), as an optional prior. The default
#' column names follow the fiboa schema; fiboa describes the data structure and
#' is not another name for FTW. \code{\link{diagnose_against_ftw}} reports
#' agreement with a field boundary reference. Accuracy claims require a
#' separate, held-out reference when the attached polygon influenced the
#' delineation.
#'
#' Whether the prior improves field boundary delineation requires independent
#' validation. This is therefore an experimental feature.
#'
#' @keywords internal
#' @importFrom stats median
"_PACKAGE"

## The CQL2 filter expressions passed to rstac::ext_filter() reference DSL
## symbols that rstac evaluates by non-standard evaluation; they are not
## R-level bindings. Declare them so R CMD check does not flag them as
## undefined.
utils::globalVariables(c("collection", "datetime", "geometry",
                         "t_intersects", "s_intersects", "eo:cloud_cover"))

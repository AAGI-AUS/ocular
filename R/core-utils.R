# =========================================================================
# Core: shared constants and utilities
# =========================================================================

# Derive unsupported fusion indices when called, after both band lookups exist.
.fusionUnsupportedIndices <- function(){

  ln <- tryCatch(names(landsat_index_list), error = function(e){
    return(character(0))
  })
  mn <- tryCatch(names(mcd43a4_index_list), error = function(e){
    return(character(0))
  })
  return(setdiff(ln, mn))
}

.S2_START <- as.Date("2015-07-01")
.L8_START <- as.Date("2013-04-01")
.STAC_ENDPOINT <- "https://planetarycomputer.microsoft.com/api/stac/v1"

`%||%` <- function(a, b){

  if( is.null(a) ) { b } else { a }
}

# Ratio indices can produce NaN or Inf when their denominator is zero.
.finiteIndex <- function(r){

  return(terra::ifel(is.finite(r), r, NA))
}

# Serialisation preserves names, types and nested values in cache keys.
# Four rolling hashes avoid temporary files and an extra package dependency.
.stableHash <- function(payload){

  bytes <- as.integer(serialize(payload, NULL, version = 2))
  modulus <- c(2147483629, 2147483587, 2147483579, 2147483563)
  multiplier <- c(257, 263, 269, 271)
  seed <- c(17, 29, 43, 61)
  pieces <- vapply(seq_along(modulus), function(i){
    h <- seed[[i]]
    for( b in bytes ){
      h <- (h * multiplier[[i]] + b + 1) %% modulus[[i]]
    }
    return(sprintf("%08x", as.integer(h)))
  }, character(1L))
  return(paste0(pieces, collapse = ""))
}

.batch_test_sites <- function(ids = c("a", "b")){
  data.frame(id = ids, longitude = 117 + seq_along(ids) / 10,
              latitude = -32, start_date = "2021-01-01", end_date = "2021-02-01")
}

.batch_test_rs <- function(longitude = 117.1){
  rs <- ocular:::.newOcular()
  rs$geom$source <- "landsat-8"
  rs$geom$longitude <- longitude
  rs$scenes <- list(list())
  return(rs)
}

.batch_test_pipeline <- function(rs){
  data.frame(date = as.Date("2021-01-05"), value = rs$geom$longitude)
}

.batch_test_transient <- function(message = "Service unavailable"){
  errorCondition(message, class = "ocular_stac_error", retryable = TRUE)
}

.batch_test_interrupt <- function(){
  stop(structure(list(message = "test interrupt", call = NULL),
                 class = c("interrupt", "condition")))
}

test_that("transient STAC classification is narrow and uses structured conditions", {
  for( status in c(408, 429, 500, 502, 503, 504) ){
    expect_true(ocular:::.isTransientStacError(
      errorCondition("request failed", status_code = status)))
    expect_true(ocular:::.isTransientStacError(
      simpleError(paste0("HTTP status '", status, "'. unavailable"))))
  }
  expect_true(ocular:::.isTransientStacError(errorCondition(
    "request failed", class = c("http_503", "http_500", "http_error"))))
  expect_true(ocular:::.isTransientStacError(errorCondition(
    "request failed", response = list(status_code = 429))))
  expect_true(ocular:::.isTransientStacError(errorCondition(
    "request failed", class = "curl_error_operation_timedout")))
  for( status in c(400, 401, 403, 404, 422, 501, 505) ){
    expect_false(ocular:::.isTransientStacError(
      simpleError(paste0("HTTP status '", status, "'. final"))))
  }
  for( msg in c("transient", "HTTP status 503: not rstac format",
                "user pipeline: HTTP status '503'. error",
                "Timeout was reached", "No scenes returned for this query.") ){
    expect_false(ocular:::.isTransientStacError(simpleError(msg)))
  }
  prefix <- "Error while requesting 'https://example.invalid/stac'. \n"
  expect_true(ocular:::.isTransientStacError(simpleError(
    paste0(prefix, "Timeout was reached [example.invalid]: operation timed out"))))
  expect_true(ocular:::.isTransientStacError(simpleError(
    paste0(prefix, "Could not resolve hostname [example.invalid]: lookup failed"))))
  expect_false(ocular:::.isTransientStacError(simpleError(
    paste0(prefix, "SSL peer certificate or SSH remote key was not OK"))))
})

test_that("batch rounds remove the nested STAC retries", {
  visits <- numeric()
  clouds <- numeric()
  pauses <- numeric()
  requests <- 0L
  testthat::local_mocked_bindings(
    get_rs = function(..., longitude, max_cloud_cover){
      visits <<- c(visits, longitude)
      clouds <<- c(clouds, max_cloud_cover)
      ocular:::.fetchStacItems("https://example.invalid/stac", "collection",
                               c(0, 0, 1, 1), "2021-01-01", "2021-02-01",
                               max_cloud_cover, TRUE,
                               initial_pause = 0, max_pause = 0)
    },
    .batchPause = function(round){ pauses <<- c(pauses, round) },
    .package = "ocular")
  testthat::local_mocked_bindings(post_request = function(...){
    requests <<- requests + 1L
    stop("HTTP status '503'. Service unavailable")
  }, .package = "rstac")
  out <- batch_rs(.batch_test_sites(), cache = FALSE, quiet = TRUE, retries = 2L)
  expect_equal(visits, c(117.1, 117.2, 117.1, 117.2))
  expect_equal(clouds, c(50, 50, 80, 80))
  expect_equal(requests, 4L)
  expect_equal(pauses, 1)
  expect_equal(attr(out, "failures")$attempts, c(2, 2))
  expect_match(attr(out, "failures")$last_error, "after 1 attempt")
})

test_that("rounds retry service errors and keep site order, final failures and empty successes", {
  sites <- .batch_test_sites(letters[1:7])
  sites$longitude[[4L]] <- Inf
  visits <- integer()
  seen <- integer(7L)
  pipelines <- integer()
  testthat::local_mocked_bindings(
    get_rs = function(..., longitude){
      i <- round((longitude - 117) * 10)
      visits <<- c(visits, i)
      seen[[i]] <<- seen[[i]] + 1L
      if( i == 1L && seen[[i]] < 3L ) stop(.batch_test_transient())
      if( i == 3L ) stop("invalid request")
      if( i == 5L ) stop("No scenes returned for this query.")
      .batch_test_rs(longitude)
    }, .batchPause = function(round) NULL, .package = "ocular")
  pipeline <- function(rs){
    i <- round((rs$geom$longitude - 117) * 10)
    pipelines <<- c(pipelines, i)
    if( i == 6L ) stop(.batch_test_transient("pipeline failure"))
    if( i == 7L ) return(data.frame(date = character(), value = numeric()))
    .batch_test_pipeline(rs)
  }
  out <- batch_rs(sites, cache = FALSE, quiet = TRUE, pipeline = pipeline)
  expect_equal(visits, c(1, 2, 3, 5, 6, 7, 1, 1))
  expect_equal(pipelines, c(2, 6, 7, 1))
  expect_equal(out$id, c("a", "b"))
  expect_equal(attr(out, "failures")$id, c("c", "d", "e", "f"))
  expect_equal(attr(out, "failures")$attempts, c(1, 0, 1, 1))
})

test_that("infinite rounds retain the requested cloud ceiling and return on interrupt", {
  visits <- numeric()
  clouds <- numeric()
  pauses <- numeric()
  count_a <- 0L
  testthat::local_mocked_bindings(
    get_rs = function(..., longitude, max_cloud_cover){
      visits <<- c(visits, longitude)
      clouds <<- c(clouds, max_cloud_cover)
      if( longitude == 117.1 ){
        count_a <<- count_a + 1L
        if( count_a == 3L ) return(.batch_test_rs(longitude))
      }
      stop(.batch_test_transient("last transient error"))
    },
    .batchPause = function(round){
      pauses <<- c(pauses, round)
      if( round == 3 ) .batch_test_interrupt()
    }, .package = "ocular")
  out <- batch_rs(.batch_test_sites(), retries = Inf, cloud_relax = 90,
                   cache = FALSE, quiet = TRUE, pipeline = .batch_test_pipeline,
                   on_error = "stop")
  expect_equal(visits, rep(c(117.1, 117.2), 3))
  expect_equal(clouds, rep(50, 6))
  expect_equal(pauses, 1:3)
  expect_equal(out$id, "a")
  failure <- attr(out, "failures")
  expect_equal(failure$id, "b")
  expect_equal(failure$attempts, 3)
  expect_equal(failure$last_error, "last transient error")
  expect_match(failure$reason, "Interrupted")
})

test_that("on_error stop waits for transient exhaustion but stops at final errors", {
  attempts <- 0L
  pauses <- 0L
  testthat::local_mocked_bindings(
    get_rs = function(...){
      attempts <<- attempts + 1L
      stop(.batch_test_transient())
    }, .batchPause = function(round){ pauses <<- pauses + 1L },
    .package = "ocular")
  expect_error(batch_rs(.batch_test_sites("a"), retries = 2, cache = FALSE,
                         quiet = TRUE, on_error = "stop"), "Service unavailable")
  expect_equal(attempts, 2L)
  expect_equal(pauses, 1L)
  testthat::local_mocked_bindings(get_rs = function(...){
    attempts <<- attempts + 1L
    stop("No scenes returned for this query.")
  }, .package = "ocular")
  expect_error(batch_rs(.batch_test_sites(), retries = Inf, cache = FALSE,
                         quiet = TRUE, on_error = "stop"), "No scenes")
  expect_equal(attempts, 3L)
})

test_that("interrupts during retrieval or the pipeline preserve completed rows", {
  for( stage in c("retrieval", "pipeline") ){
    local({
      testthat::local_mocked_bindings(get_rs = function(..., longitude){
        if( stage == "retrieval" && longitude == 117.2 ) .batch_test_interrupt()
        .batch_test_rs(longitude)
      }, .package = "ocular")
      pipeline <- function(rs){
        if( stage == "pipeline" && rs$geom$longitude == 117.2 ) .batch_test_interrupt()
        .batch_test_pipeline(rs)
      }
      out <- batch_rs(.batch_test_sites(c("a", "b", "c")), cache = FALSE,
                       quiet = TRUE, pipeline = pipeline, on_error = "stop")
      expect_equal(out$id, "a")
      expect_equal(attr(out, "failures")$id, c("b", "c"))
      expect_equal(attr(out, "failures")$attempts, c(1, 0))
      expect_equal(attr(out, "failures")$last_error, c("", ""))
    })
  }
})

test_that("interrupts during cache writes keep the just-completed site", {
  cache_dir <- .test_artifact_dir("batch-interrupt-write")
  testthat::local_mocked_bindings(
    get_rs = function(..., longitude) .batch_test_rs(longitude),
    .package = "ocular")
  testthat::local_mocked_bindings(saveRDS = function(...) .batch_test_interrupt(),
                                  .package = "base")
  out <- batch_rs(.batch_test_sites(), cache_dir = cache_dir, quiet = TRUE,
                   pipeline = .batch_test_pipeline, pipeline_cache_key = "interrupt")
  expect_equal(out$id, "a")
  expect_equal(attr(out, "failures")$id, "b")
  expect_equal(attr(out, "failures")$attempts, 0)
})

test_that("cache hits survive interrupts while reading a later site's entry", {
  cache_dir <- .test_artifact_dir("batch-interrupt-read")
  fetched <- 0L
  testthat::local_mocked_bindings(get_rs = function(..., longitude){
    fetched <<- fetched + 1L
    .batch_test_rs(longitude)
  }, .package = "ocular")
  run <- function() batch_rs(.batch_test_sites(), cache_dir = cache_dir,
                              quiet = TRUE, pipeline = .batch_test_pipeline,
                              pipeline_cache_key = "read-interrupt")
  first <- run()
  real_read <- ocular:::.rsCacheRead
  reads <- 0L
  testthat::local_mocked_bindings(.rsCacheRead = function(path){
    reads <<- reads + 1L
    if( reads == 2L ) .batch_test_interrupt()
    real_read(path)
  }, .package = "ocular")
  out <- run()
  expect_equal(fetched, 2L)
  expect_equal(out$id, first$id[[1L]])
  expect_equal(attr(out, "failures")$attempts, 0)
})

test_that("completed results retain attempts in cache and are reused without retrieval", {
  cache_dir <- .test_artifact_dir("batch-rounds-cache")
  attempts <- 0L
  testthat::local_mocked_bindings(get_rs = function(..., longitude){
    attempts <<- attempts + 1L
    if( attempts < 3L ) stop(.batch_test_transient())
    .batch_test_rs(longitude)
  }, .batchPause = function(round) NULL, .package = "ocular")
  run <- function() batch_rs(.batch_test_sites("a"), cache_dir = cache_dir,
                              quiet = TRUE, pipeline = .batch_test_pipeline,
                              pipeline_cache_key = "rounds-cache")
  first <- run()
  entries <- ocular:::.cacheEntryFiles(cache_dir, "^[0-9a-f]{24}[.]rds$")
  expect_length(entries, 1L)
  meta <- readRDS(entries[[1L]])$meta
  expect_equal(meta$attempts, 3)
  expect_equal(meta$cloud_used, 80)
  expect_equal(meta$retry_policy, "rounds-transient-v1")
  second <- run()
  expect_equal(attr(first, "outcomes")$attempts, 3)
  expect_equal(attr(second, "outcomes")$attempts, 0)
  expect_true(attr(second, "outcomes")$cached)
  attr(first, "outcomes") <- NULL
  attr(second, "outcomes") <- NULL
  expect_equal(first, second)
  expect_equal(attempts, 3L)
})

test_that("batch query option is restored before pipeline and after every exit", {
  for( prior in list(NULL, FALSE, "existing value") ){
    local({
      old <- options(ocular.batch_single_attempt = prior)
      on.exit(options(old))
      for( outcome in c("success", "error", "interrupt") ){
        local({
          testthat::local_mocked_bindings(get_rs = function(..., longitude){
            expect_true(getOption("ocular.batch_single_attempt"))
            if( outcome == "error" ) stop("bad input")
            if( outcome == "interrupt" ) .batch_test_interrupt()
            .batch_test_rs(longitude)
          }, .package = "ocular")
          pipeline <- function(rs){
            expect_identical(getOption("ocular.batch_single_attempt"), prior)
            .batch_test_pipeline(rs)
          }
          batch_rs(.batch_test_sites("a"), cache = FALSE, quiet = TRUE,
                     pipeline = pipeline)
          expect_identical(getOption("ocular.batch_single_attempt"), prior)
        })
      }
    })
  }
})

test_that("failure attributes, cache warnings, and STAC conditions redact URL queries", {
  secret <- paste0("synthetic", "-credential")
  url <- paste0("https://example.invalid/asset", "?", "sig", "=", secret)
  testthat::local_mocked_bindings(get_rs = function(...) stop(url), .package = "ocular")
  out <- batch_rs(.batch_test_sites("a"), cache = FALSE, quiet = TRUE)
  expect_false(any(grepl(secret, unlist(attr(out, "failures")), fixed = TRUE)))
  expect_error(batch_rs(.batch_test_sites("a"), cache = FALSE, quiet = TRUE,
                         on_error = "stop"), "https://example.invalid/asset$")
  testthat::local_mocked_bindings(post_request = function(...){
    stop(paste0("HTTP status '503'. ", url))
  }, .package = "rstac")
  e <- tryCatch(ocular:::.fetchStacItems("https://example.invalid/stac",
    "collection", c(0, 0, 1, 1), "2021-01-01", "2021-02-01", 50, TRUE,
    retry_max = 1L), error = identity)
  expect_s3_class(e, "ocular_stac_error")
  expect_s3_class(e, "simpleError")
  expect_true(e$retryable)
  expect_false(any(grepl(secret, unlist(e), fixed = TRUE)))
  testthat::local_mocked_bindings(
    get_rs = function(..., longitude) .batch_test_rs(longitude),
    .package = "ocular")
  testthat::local_mocked_bindings(saveRDS = function(...) stop(url),
                                  .package = "base")
  warning_text <- character()
  out <- withCallingHandlers(
    batch_rs(.batch_test_sites("a"),
               cache_dir = .test_artifact_dir("batch-redacted-warning"),
               quiet = TRUE, pipeline = .batch_test_pipeline,
               pipeline_cache_key = "redaction"),
    warning = function(w){
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_length(warning_text, 1L)
  expect_false(grepl(secret, warning_text, fixed = TRUE))
  expect_equal(out$id, "a")
})

test_that("backoff grows between rounds and remains capped", {
  seconds <- numeric()
  testthat::local_mocked_bindings(Sys.sleep = function(time){
    seconds <<- c(seconds, time)
  }, .package = "base")
  invisible(lapply(c(1:8, 1e6), ocular:::.batchPause))
  expect_equal(seconds, c(2, 4, 8, 16, 32, 60, 60, 60, 60))
})

test_that("cache identity separates finite and unlimited cloud and retry policies", {
  key <- function(retries, cloud_relax = 80){
    ocular:::.rsCacheKey(117.1, -32, "2021-01-01", "2021-02-01", "EVI2",
                          "auto", 1000, 1000, 50, rs_params(), retries,
                          cloud_relax, "policy")
  }
  expect_false(identical(key(3L), key(Inf)))
  expect_false(identical(key(Inf, 80), key(Inf, 90)))
  real_hash <- ocular:::.stableHash
  query <- NULL
  testthat::local_mocked_bindings(.stableHash = function(x){
    query <<- x
    real_hash(x)
  }, .package = "ocular")
  key(Inf)
  expect_identical(query$relaxed_final_attempt, FALSE)
  expect_identical(query$retry_policy, "rounds-transient-v1")
  key(3L)
  expect_identical(query$relaxed_final_attempt, TRUE)
})

test_that("empty pipeline results are successful and cached", {
  cache_dir <- .test_artifact_dir("batch-rounds-empty")
  attempts <- 0L
  testthat::local_mocked_bindings(get_rs = function(..., longitude){
    attempts <<- attempts + 1L
    .batch_test_rs(longitude)
  }, .package = "ocular")
  run <- function() batch_rs(.batch_test_sites("a"), cache_dir = cache_dir,
                              quiet = TRUE, pipeline_cache_key = "empty",
                              pipeline = function(rs){
                                data.frame(date = character(), value = numeric())
                              })
  first <- run()
  second <- run()
  expect_equal(nrow(first), 0L)
  expect_null(attr(first, "failures"))
  expect_identical(attr(first, "outcomes")$status, "empty")
  expect_identical(attr(first, "outcomes")$dates, 0L)
  expect_equal(attr(first, "outcomes")$attempts, 1)
  expect_equal(attr(second, "outcomes")$attempts, 0)
  expect_true(attr(second, "outcomes")$cached)
  attr(first, "outcomes") <- NULL
  attr(second, "outcomes") <- NULL
  expect_identical(first, second)
  expect_equal(attempts, 1L)
})

test_that("interrupts during cache setup account for all unvisited sites", {
  testthat::local_mocked_bindings(.rsCacheDir = function(...) .batch_test_interrupt(),
                                  .package = "ocular")
  out <- batch_rs(.batch_test_sites(), quiet = TRUE)
  expect_equal(nrow(out), 0L)
  expect_equal(attr(out, "failures")$attempts, c(0, 0))
  expect_equal(attr(out, "failures")$id, c("a", "b"))
  expect_identical(attr(out, "outcomes")$status, c("interrupted", "interrupted"))
  expect_identical(attr(out, "outcomes")$dates, rep(NA_integer_, 2L))
  expect_identical(attr(out, "outcomes")$attempts, c(0, 0))
  expect_identical(attr(out, "outcomes")$cached, c(FALSE, FALSE))
})

test_that("unreadable-cache progress is retained once per site and respects quiet", {
  cache_dir <- .test_artifact_dir("batch-corrupt-rounds")
  attempts <- 0L
  testthat::local_mocked_bindings(get_rs = function(..., longitude){
    attempts <<- attempts + 1L
    if( attempts %% 2L == 0L ) stop(.batch_test_transient())
    .batch_test_rs(longitude)
  }, .batchPause = function(round) NULL, .package = "ocular")
  run <- function(quiet) batch_rs(.batch_test_sites("a"), cache_dir = cache_dir,
                                   quiet = quiet, pipeline = .batch_test_pipeline,
                                   pipeline_cache_key = "corrupt")
  run(TRUE)
  entry <- ocular:::.cacheEntryFiles(cache_dir, "^[0-9a-f]{24}[.]rds$")
  writeLines("corrupt cache entry", entry)
  messages <- character()
  withCallingHandlers(run(FALSE), message = function(m){
    messages <<- c(messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  expect_equal(sum(grepl("ignoring unreadable cache entry for", messages)), 1L)
  expect_equal(attempts, 3L)
  writeLines("corrupt cache entry", entry)
  expect_message(run(TRUE), NA)
  expect_equal(attempts, 5L)
})

test_that("unlimited retries admit only positive real infinity", {
  for( budget in list(Inf + 0i, 3 + 0i, NaN, -Inf, c(Inf, Inf), "Inf") ){
    expect_error(batch_rs(.batch_test_sites(), retries = budget,
                           cache = FALSE, quiet = TRUE),
                   "single positive integer or Inf")
  }
})

test_that("batch outcomes retain every input row and distinguish empty results", {
  sites <- .batch_test_sites(c("duplicate", "duplicate", "invalid", "failed", "active", "unvisited"))
  sites$longitude[[3L]] <- Inf
  testthat::local_mocked_bindings(get_rs = function(..., longitude){
    if( longitude == 117.4 ){
      stop("invalid request")
    }
    if( longitude == 117.5 ){
      .batch_test_interrupt()
    }
    .batch_test_rs(longitude)
  }, .package = "ocular")
  pipeline <- function(rs){

    if( rs$geom$longitude == 117.1 ){
      return(data.frame(date = character(), value = numeric()))
    }
    .batch_test_pipeline(rs)
  }
  out <- batch_rs(sites, pipeline = pipeline, cache = FALSE, quiet = TRUE)
  outcomes <- attr(out, "outcomes")
  expect_named(outcomes, c("id", "index", "status", "dates", "attempts", "cached"))
  expect_identical(outcomes$id, sites$id)
  expect_identical(outcomes$index, rep("EVI2", 6L))
  expect_identical(outcomes$status,
                   c("empty", "success", "failure", "failure", "interrupted", "interrupted"))
  expect_identical(outcomes$dates, c(0L, 1L, rep(NA_integer_, 4L)))
  expect_identical(outcomes$attempts, c(1, 1, 0, 1, 1, 0))
  expect_identical(outcomes$cached, rep(FALSE, 6L))
  expect_identical(out$id, "duplicate")
  expect_equal(out$value, 117.2)
  expect_identical(attr(out, "failures")$id, sites$id[3:6])
})

test_that("old cache data produces empty outcomes without extra metadata", {
  cache_dir <- .test_artifact_dir("batch-outcomes-old-cache")
  calls <- 0L
  testthat::local_mocked_bindings(get_rs = function(..., longitude){
    calls <<- calls + 1L
    .batch_test_rs(longitude)
  }, .package = "ocular")
  run <- function(quiet = TRUE){

    batch_rs(.batch_test_sites(), pipeline = function(rs){
      data.frame(date = character(), value = numeric())
    }, pipeline_cache_key = "old-empty", cache_dir = cache_dir, quiet = quiet)
  }
  first <- run()
  expect_identical(attr(first, "outcomes")$status, c("empty", "empty"))
  expect_equal(nrow(first), 0L)
  expect_null(attr(first, "failures"))
  entries <- ocular:::.cacheEntryFiles(cache_dir, "^[0-9a-f]{24}[.]rds$")
  for( path in entries ){
    entry <- readRDS(path)
    expect_null(attr(entry$data, "outcomes"))
    expect_false("outcomes" %in% names(entry))
    saveRDS(list(data = entry$data), path)
  }
  messages <- character()
  second <- withCallingHandlers(run(FALSE), message = function(m){
    messages <<- c(messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  expect_equal(calls, 2L)
  expect_identical(attr(second, "outcomes")$status, c("empty", "empty"))
  expect_identical(attr(second, "outcomes")$dates, c(0L, 0L))
  expect_identical(attr(second, "outcomes")$attempts, c(0, 0))
  expect_identical(attr(second, "outcomes")$cached, c(TRUE, TRUE))
  expect_equal(sum(grepl("(cached, 0 dates)", messages, fixed = TRUE)), 2L)
  attr(first, "outcomes") <- NULL
  attr(second, "outcomes") <- NULL
  expect_identical(first, second)
})

test_that("completed empty outcomes survive cache-write and progress interrupts", {
  for( stage in c("write", "progress", "cached-progress") ){
    local({
      cache_dir <- .test_artifact_dir(paste0("batch-empty-", stage))
      testthat::local_mocked_bindings(
        get_rs = function(..., longitude) .batch_test_rs(longitude), .package = "ocular")
      pipeline <- function(rs){

        return(data.frame(date = character(), value = numeric()))
      }
      if( stage == "cached-progress" ){
        batch_rs(.batch_test_sites("a"), cache_dir = cache_dir, quiet = TRUE,
                   pipeline = pipeline, pipeline_cache_key = "empty-interrupt")
      }
      if( stage == "write" ){
        testthat::local_mocked_bindings(saveRDS = function(...) .batch_test_interrupt(),
                                        .package = "base")
      }
      if( stage != "write" ){
        testthat::local_mocked_bindings(message = function(...) .batch_test_interrupt(),
                                        .package = "base")
      }
      out <- batch_rs(.batch_test_sites(), cache_dir = cache_dir,
                        cache = stage != "progress", quiet = stage == "write",
                        pipeline = pipeline, pipeline_cache_key = "empty-interrupt",
                        on_error = "stop")
      expect_equal(nrow(out), 0L)
      expect_identical(attr(out, "outcomes")$status, c("empty", "interrupted"))
      expect_identical(attr(out, "outcomes")$dates, c(0L, NA_integer_))
      expect_identical(attr(out, "outcomes")$attempts,
                       c(if( stage == "cached-progress" ) { 0 } else { 1 }, 0))
      expect_identical(attr(out, "outcomes")$cached, c(stage == "cached-progress", FALSE))
      expect_identical(attr(out, "failures")$id, "b")
    })
  }
})

test_that("an interrupted cache hit is not complete before its rows are stored", {
  cache_dir <- .test_artifact_dir("batch-hit-before-commit")
  sites <- .batch_test_sites()
  testthat::local_mocked_bindings(
    get_rs = function(..., longitude) .batch_test_rs(longitude), .package = "ocular")
  pipeline <- function(rs){

    return(data.frame(date = character(), value = numeric()))
  }
  batch_rs(sites[1L, , drop = FALSE], cache_dir = cache_dir, quiet = TRUE,
             pipeline = pipeline, pipeline_cache_key = "commit")
  out <- testthat::with_mocked_bindings(
    batch_rs(sites, cache_dir = cache_dir, quiet = TRUE,
               pipeline = pipeline, pipeline_cache_key = "commit"),
    cbind = function(...) .batch_test_interrupt(), .package = "base")
  expect_identical(attr(out, "outcomes")$status, c("interrupted", "interrupted"))
  expect_identical(attr(out, "outcomes")$dates, rep(NA_integer_, 2L))
  expect_identical(attr(out, "outcomes")$attempts, c(0, 0))
  expect_identical(attr(out, "outcomes")$cached, c(FALSE, FALSE))
  expect_identical(attr(out, "failures")$id, c("a", "b"))
})

test_that("failure progress interrupts preserve the completed failure outcome", {
  testthat::local_mocked_bindings(get_rs = function(...) stop("invalid request"),
                                  .package = "ocular")
  out <- testthat::with_mocked_bindings(
    batch_rs(.batch_test_sites(), cache = FALSE),
    message = function(...) .batch_test_interrupt(), .package = "base")
  expect_identical(attr(out, "outcomes")$status, c("failure", "interrupted"))
  expect_identical(attr(out, "outcomes")$dates, rep(NA_integer_, 2L))
  expect_identical(attr(out, "outcomes")$attempts, c(1, 0))
  expect_identical(attr(out, "outcomes")$cached, c(FALSE, FALSE))
  expect_identical(attr(out, "failures")$reason,
                   c("invalid request", "Interrupted before completion."))
})

test_that("outcomes count current attempts when another site supplies the cache", {
  sites <- .batch_test_sites()
  sites$longitude[[2L]] <- sites$longitude[[1L]]
  calls <- 0L
  testthat::local_mocked_bindings(get_rs = function(..., longitude){
    calls <<- calls + 1L
    if( calls == 1L ){
      stop(.batch_test_transient())
    }
    .batch_test_rs(longitude)
  }, .batchPause = function(round) NULL, .package = "ocular")
  out <- batch_rs(sites, cache_dir = .test_artifact_dir("batch-shared-outcome-cache"),
                   pipeline = .batch_test_pipeline, pipeline_cache_key = "shared",
                   quiet = TRUE)
  expect_equal(calls, 2L)
  expect_identical(out$id, c("a", "b"))
  expect_identical(attr(out, "outcomes")$status, c("success", "success"))
  expect_identical(attr(out, "outcomes")$dates, c(1L, 1L))
  expect_identical(attr(out, "outcomes")$attempts, c(1, 1))
  expect_identical(attr(out, "outcomes")$cached, c(TRUE, FALSE))
})

#!/usr/bin/env Rscript

# Record the checked tarball and the fresh installed ocular path immediately
# after R CMD INSTALL. This script installs nothing and makes no network call.

.enabled <- function(x) {
  identical(tolower(trimws(x)), "true")
}

.assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

.inside <- function(path, root, must_work = FALSE) {
  path <- normalizePath(path, winslash = "/", mustWork = must_work)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  identical(path, root) || startsWith(path, paste0(root, "/"))
}

.is_symlink <- function(path) {
  target <- Sys.readlink(path)
  length(target) == 1L && !is.na(target) && nzchar(target)
}

.required_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  .assert(nzchar(value), paste0("Missing required environment variable `",
                                name, "`."))
  value
}

.run <- function() {
  .assert(
    .enabled(Sys.getenv("OCULAR_WRITE_INSTALL_RECEIPT", unset = "")),
    paste0("Receipt creation is opt-in. Set ",
           "OCULAR_WRITE_INSTALL_RECEIPT=true as shown in ",
           "tests/manual/README.md.")
  )

  package_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  workspace_root <- dirname(package_root)
  description_path <- file.path(package_root, "DESCRIPTION")
  .assert(isTRUE(file_test("-f", description_path)),
          "Run this script from the ocular package root.")
  description <- read.dcf(description_path)
  .assert(identical(unname(description[1L, "Package"]), "ocular"),
          "The working directory is not the ocular package root.")
  source_version <- unname(description[1L, "Version"])
  .assert(.inside(tempdir(), package_root, must_work = TRUE),
          paste0("R's session tempdir is outside the ocular package. Set ",
                 "TMPDIR as shown in tests/manual/README.md before starting R."))

  raw_tarball <- .required_env("OCULAR_CHECKED_TARBALL")
  expected_md5 <- .required_env("OCULAR_CHECKED_TARBALL_MD5")
  raw_library <- .required_env("OCULAR_VALIDATION_LIBRARY")
  raw_receipt <- .required_env("OCULAR_INSTALL_RECEIPT")

  .assert(!.is_symlink(raw_tarball),
          "OCULAR_CHECKED_TARBALL must not be a symbolic link.")
  tarball <- normalizePath(raw_tarball, winslash = "/", mustWork = TRUE)
  .assert(isTRUE(file_test("-f", tarball)) &&
            .inside(tarball, workspace_root, must_work = TRUE),
          "The checked tarball must be a regular file inside the workspace.")
  .assert(identical(basename(tarball),
                    paste0("ocular_", source_version, ".tar.gz")),
          "The checked tarball name does not match the source version.")
  tarball_md5 <- unname(tools::md5sum(tarball))
  .assert(length(tarball_md5) == 1L && !is.na(tarball_md5) &&
            identical(tarball_md5, expected_md5),
          paste0("The checked tarball checksum changed after it was built; ",
                 "do not create an install receipt."))

  .assert(!.is_symlink(raw_library),
          "OCULAR_VALIDATION_LIBRARY must not be a symbolic link.")
  validation_library <- normalizePath(
    raw_library, winslash = "/", mustWork = TRUE)
  .assert(dir.exists(validation_library) &&
            .inside(validation_library, package_root, must_work = TRUE),
          "The validation library must be a directory inside the package.")
  raw_installed_path <- file.path(validation_library, "ocular")
  .assert(!.is_symlink(raw_installed_path),
          "The installed ocular directory must not be a symbolic link.")
  installed_path <- normalizePath(
    raw_installed_path, winslash = "/", mustWork = TRUE)
  .assert(identical(dirname(installed_path), validation_library),
          "The installed ocular directory is not in the validation library.")
  installed_description <- read.dcf(
    file.path(installed_path, "DESCRIPTION"))
  installed_version <- unname(installed_description[1L, "Version"])
  .assert(identical(unname(installed_description[1L, "Package"]), "ocular") &&
            identical(installed_version, source_version),
          "The installed package identity does not match the ocular source.")

  .assert(!file.exists(raw_receipt) && !dir.exists(raw_receipt) &&
            !.is_symlink(raw_receipt),
          "Refusing to replace an existing install receipt.")
  receipt <- normalizePath(raw_receipt, winslash = "/", mustWork = FALSE)
  .assert(.inside(receipt, package_root) &&
            identical(dirname(receipt), dirname(validation_library)),
          paste0("The receipt must be in the timestamped package directory ",
                 "that contains the validation library."))

  utils::write.csv(data.frame(
    item = c(
      "checked_tarball", "checked_tarball_md5",
      "installed_package_version", "installed_package_path",
      "validation_library"
    ),
    value = c(
      tarball, tarball_md5, installed_version, installed_path,
      validation_library
    ),
    stringsAsFactors = FALSE
  ), receipt, row.names = FALSE)
  .assert(isTRUE(file_test("-f", receipt)) && !.is_symlink(receipt),
          "The install receipt was not written as a regular file.")
  cat("Install receipt written to:\n", receipt, "\n", sep = "")
  invisible(receipt)
}

.run()

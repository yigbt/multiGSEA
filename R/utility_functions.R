#' Retrieve path to a cached file.
#'
#' The function retrieves the path to a file that is cached in the archive
#' directory.
#'
#' @param filename Name of the file.
#'
#' @return String containing the path to the file.
archivePath <- function(filename) {
  ad <- archiveDir()
  filename <- paste0(filename, ".rds")

  return(file.path(ad, filename, fsep = .Platform$file.sep))
}





#' Retrieve the path to the cache directory.
#'
#' Retrieve the path to the cache directory for the multiGSEA package.
#' Create the cache directory if need be.
#'
#' @return String containing the path to the cache directory.
#'
#' @importFrom rappdirs user_cache_dir
archiveDir <- function() {
  ad <- rappdirs::user_cache_dir("multiGSEA")

  if (!file.exists(ad)) {
    if (!dir.create(ad, FALSE, TRUE)) {
      stop("An error occurred during creating the archive directory: ", ad,
        call. = FALSE
      )
    }
  }

  return(ad)
}





#' Read a local RDS file.
#'
#' Use the readRDS function to load the given file which should be in RDS
#' format.
#'
#' @param filename Path to the file to be read.
#'
#' @return Content of file.
#'
#' @importFrom methods is
loadLocal <- function(filename) {
  res <- try(readRDS(filename), silent = TRUE)

  if ("try-error" %in% methods::is(res)) {
    return(NULL)
  } else {
    return(res)
  }
}


#' Make a list of strings unique
#'
#' It might happen that there are duplicated strings in a list. With this
#' function we will rename those duplicated entries in a way that we simply add
#' the number of occurrences to the string. I.e., when the string foo occurs
#' three times in a list, it will be renamed to foo_1, foo_2, and foo_3,
#' respectively.
#'
#' @param names List of strings where duplicates should be renamed
#'
#' @return List where duplicates are renamed.
#'
#' @examples
#' l <- c("foo", "bar", "foo", "bars")
#' rename_duplicates(l)
#'
#' @export
rename_duplicates <- function(names) {
  tn <- table(names)
  for (name in names(tn[tn > 1])) {
    names[names == name] <- paste0(name, "_", 1:tn[names(tn) == name])
  }

  return(names)
}

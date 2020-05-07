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
    filename <- paste0( filename, ".rds")
    
    return( file.path( ad, filename, fsep = .Platform$file.sep))
    
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
    
    ad <- user_cache_dir("multiGSEA")

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
#' Use the readRDS function to load read the given file which should be in RDS
#' format.
#'
#' @param filename Path to the file to be read.
#'
#' @return Content of file.
#'
#' @importFrom methods is
loadLocal <- function(filename) {
    
    res <- try(readRDS(filename), silent = TRUE)

    if ("try-error" %in% is(res)) {
        return(NULL)
    } else {
        return(res)
    }
    
}

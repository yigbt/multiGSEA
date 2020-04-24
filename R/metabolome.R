#' Metabolomic data set that is used in the toy example provided by the `multiGSEA` package.
#'
#' Processed metabolomics data set that will be used throughout
#' the vignette provided by the `multiGSEA` package. The raw data
#' was originally published by [Quiros _et al._](http://doi.org/10.1083/jcb.201702058)
#' and can be accessed within the online supplementary material.
#'
#' @docType data
#'
#' @usage data(metabolome)
#'
#' @format A tibble with 4 variables and 4881 measured proteome features:
#' \describe{
#'    \item{HMDB}{HMDB identifier of measured metabolites.}
#'    \item{logFC}{Log2-transformed fold change between treatment and control.}
#'    \item{pValue}{P-value associated with the fold change.}
#'    \item{adj.pValue}{Adjusted p-value associated with the fold change.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(metabolome)
"metabolome"

#' Proteomic data set that is used in the toy example provided by the `multiGSEA` package.
#'
#' Processed proteomics data set that will be used throughout the vignette
#' provided by the `multiGSEA` package. The raw data was originally published
#' by [Quiros _et al._](http://doi.org/10.1083/jcb.201702058) and deposited at
#' [ProteomeXchange](http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD006293).
#' 
#' @docType data
#'
#' @usage data(proteome)
#'
#' @format A data frame with 4 variables and 8275 measured proteome features:
#' \describe{
#'    \item{Symbol}{HGNC symbol of measured proteins.}
#'    \item{logFC}{Log2-transformed fold change between treatment and control.}
#'    \item{pValue}{P-value associated with the fold change.}
#'    \item{adj.pValue}{Adjusted p-value associated with the fold change.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(proteome)
"proteome"
#' Transcriptomic data set that is used in the toy example provided by the `multiGSEA` package.
#'
#' Processed transcriptomics data set that will be used throughout the
#' vignette provided by the `multiGSEA` package. The raw data was originally
#' published by [Quiros _et al._](http://doi.org/10.1083/jcb.201702058) and
#' deposited at [NCBI Geo](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84631).
#' 
#' @docType data
#'
#' @usage data(transcriptome)
#'
#' @format A dataframe with 4 variables and 15174 measured transcriptome features:
#' \describe{
#'    \item{Symbol}{HGNC symbol of measured transcripts.}
#'    \item{logFC}{Log2-transformed fold change between treatment and control.}
#'    \item{pValue}{P-value associated with the fold change.}
#'    \item{adj.pValue}{Adjusted p-value associated with the fold change.}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(transcriptome)
"transcriptome"
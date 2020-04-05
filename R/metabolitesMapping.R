#' Mapping table for metabolites, derived and compiled from four different data sources.
#'
#' Four different sourcex of annotated metabolites, i.e., \code{HMDB},
#' \code{ChEBI}, \code{CompTox}, and \code{KEGG}, have been retrieved to
#' compile a comprehensive mapping of available metabolite IDs. ID formats
#' that are present in the mapping table are: 
#' DTXCID (Comptox),
#' DTXSID (Comptox), 
#' CAS-number,
#' CID (Pubchem), 
#' HMDB,
#' ChEBI,
#' and KEGG
#' 
#' @docType data
#'
#' @usage data(metabolitesMapping)
#'
#' @format A tibble with 7 variables and more then 1 Million metabolites:
#' \describe{
#'    \item{DTXCID}{DSSTox structure identifier, character}
#'    \item{DTXSID}{DSSTox substance identifier, character}
#'    \item{CAS}{CAS registry number, character}
#'    \item{CID}{Pubchem compound identifier, character}
#'    \item{HMDB}{Human Metabolome Database identifier (new format), character}
#'    \item{ChEBI}{Chemical Entities of Biological Interest identifier, character}
#'    \item{KEGG}{KEGG Compound identifier, character}
#' }
#'
#' @keywords datasets
#'
#' @examples
#' data(metabolitesMapping)
"metabolitesMapping"
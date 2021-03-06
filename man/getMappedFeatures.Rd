% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_processing.R
\name{getMappedFeatures}
\alias{getMappedFeatures}
\title{Wrapper to get feature mappings.}
\usage{
getMappedFeatures(
  pathways,
  returnID = "SYMBOL",
  organism = "hsapiens",
  which = "proteins",
  useLocal = TRUE
)
}
\arguments{
\item{pathways}{List of pathway definitions.}

\item{returnID}{String specifying the returned ID format.}

\item{organism}{String defining the organism of analysis.}

\item{which}{Mode to extract the features, either \'proteins\' or
\'metabolites\'.}

\item{useLocal}{Boolean specifying whether or not to use the local
preprocessed mapping.}
}
\value{
List of mapped features for an omics layer.
}
\description{
Feature mappings will be used from hard disk in case they have been
mapped before and `useLocal` is not set to be FALSE.
In other cases, a feature extraction will be done and the results are
stored for a following occasion.
}

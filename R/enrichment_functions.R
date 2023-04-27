#' Calculate pathway enrichment for multiple omics layer.
#'
#' This function calculates GSEA-based enrichments scores for multiple omics
#' layer at once. Input pathways or gene sets have to be prepared in advance by
#' means of the function \code{\link[multiGSEA]{initOmicsDataStructure}}. The function uses pre-
#' ranked lists for each omics layer to calculate the enrichment score. The
#' ranking can be calculated by means of the function
#' \link[multiGSEA]{rankFeatures}.
#'
#' @param pathways Nested list containing all pathway features for the
#'   respective omics layer.
#' @param ranks Nested list containing the measured and pre-ranked features for
#'   each omics layer.
#' @param eps This parameter sets the boundary for calculating the p value.
#'
#' @return Nested list containing the enrichment scores for each given pathway
#'   and omics layer.
#'
#' @examples
#'
#' # Download pathway definition and extract features
#' pathways <- getMultiOmicsFeatures(dbs = c("kegg"), layer = c("transcriptome", "proteome"))
#'
#' # load omics data and calculate ranks
#' data(transcriptome)
#' data(proteome)
#' ranks <- initOmicsDataStructure(c("transcriptome", "proteome"))
#' ranks$transcriptome <- rankFeatures(transcriptome$logFC, transcriptome$pValue)
#' names(ranks$transcriptome) <- transcriptome$Symbol
#' ranks$proteome <- rankFeatures(proteome$logFC, proteome$pValue)
#' names(ranks$proteome) <- proteome$Symbol
#'
#' ## run the enrichment
#' multiGSEA(pathways, ranks)
#' @importFrom fgsea fgseaMultilevel
#'
#' @export
multiGSEA <- function(pathways, ranks, eps = 0) {
  # Go through all omics layer.
  es <- lapply(names(pathways), function(omics) {
    fgsea::fgseaMultilevel(pathways[[omics]], ranks[[omics]], eps = eps)
  })

  names(es) <- names(pathways)

  return(es)
}



#' Create a reshaped data frame from multiGSEA output.
#'
#' This function reshapes the output from multiGSEA to get a single data frame
#' with columns for p-values and adjusted p-values for each omics layer. Each
#' row of the data frame represents one pathway.
#'
#' @param enrichmentScores Nested List of enrichment scores, calculated by
#'   multiGSEA function.
#' @param pathwayNames List containing Pathway names.
#'
#' @return Data frame where rows are pathways and columns are (adjusted)
#'   p-values for each omics layer.
#'
#' @examples
#' # Download pathway definition and extract features
#' pathways <- getMultiOmicsFeatures(dbs = c("kegg"), layer = c("transcriptome", "proteome"))
#'
#' # load omics data and calculate ranks
#' data(transcriptome)
#' data(proteome)
#' ranks <- initOmicsDataStructure(c("transcriptome", "proteome"))
#' ranks$transcriptome <- rankFeatures(transcriptome$logFC, transcriptome$pValue)
#' names(ranks$transcriptome) <- transcriptome$Symbol
#' ranks$proteome <- rankFeatures(proteome$logFC, proteome$pValue)
#' names(ranks$proteome) <- proteome$Symbol
#'
#' # run the enrichment
#' es <- multiGSEA(pathways, ranks)
#'
#' extractPvalues(
#'   enrichmentScores = es,
#'   pathwayNames = names(pathways[[1]])
#' )
#' @export
extractPvalues <- function(enrichmentScores, pathwayNames) {
  # Go through all the pathways
  res <- lapply(pathwayNames, function(name) {
    # Go through all the possible omics layer
    unlist(lapply(names(enrichmentScores), function(y) {
      df <- enrichmentScores[[y]][which(enrichmentScores[[y]]$pathway == name), c(2, 3)]
      if (nrow(df) == 0) {
        df <- data.frame(pval = NA, padj = NA)
      }
      names(df) <- paste0(y, "_", names(df))
      df
    }))
  })

  # Combine list elements to data frame
  # and assign pathway names as rownames
  res <- data.frame(do.call(rbind, res))

  return(res)
}



#' Calculate a combined p-value for multiple omics layer.
#'
#' This function applies the Stouffer method, the Edgington method or the
#' Fisher\'s combined probability test to combine p-values of independent tests
#' that are based on the same null hypothesis. The Stouffer method can also be
#' applied in a weighted fashion.
#'
#' @param df Data frame where rows represent a certain pathway or gene set and
#'   columns represent p-values derived from independent tests, e.g., different
#'   omics layer.
#' @param method String that specifies the method to combine multiple p-values.
#'   Default: "stouffer" Options: "stouffer", "fisher", "edgington"
#' @param weights List of weights that will be used in a weighted Stouffer
#'   method.
#'
#' @return Vector of length \code{nrow(df)} with combined p-values.
#'
#' @examples
#' df <- cbind(runif(5), runif(5), runif(5))
#' colnames(df) <- c("trans.pval", "prot.pval", "meta.pval")
#'
#' # run the unweighted summation of z values
#' combinePvalues(df)
#'
#' # run the weighted variant
#' combinePvalues(df, weights = c(10, 5, 1))
#'
#' # run the Fisher's combined probability test
#' combinePvalues(df, method = "fisher")
#'
#' # run the Edgington's method
#' combinePvalues(df, method = "edgington")
#' @importFrom metap sumz sumlog sump
#'
#' @export
combinePvalues <- function(df, method = "stouffer", weights = NULL) {
  method <- tolower(method)
  if (!method %in% c("stouffer", "fisher", "edgington")) {
    stop("You can chose between the 'stouffer', 'edgington',
              and 'fisher' method to combine p-values.",
      call. = FALSE
    )
  }

  cols <- grep("padj", colnames(df))

  pvals <- apply(df, 1, function(row) {
    row <- row[cols]
    row <- row[!is.na(row)]

    if (length(row) >= 2) {
      if (method == "fisher") {
        p <- metap::sumlog(row)
        p$p
      } else if (method == "edgington") {
        p <- metap::sump(row)
        p$p
      } else {
        ## sumz allows only p-values smaller than 1
        row <- row[row > 0 & row < 1]

        if (length(row) >= 2) {
          if (length(weights) > 0) {
            p <- metap::sumz(row, weights = weights)
          } else {
            p <- metap::sumz(row)
          }
          p$p
        } else if (length(row == 1)) {
          row[1]
        } else {
          NA
        }
      }
    } else if (length(row) == 1) {
      row[1]
    } else {
      NA
    }
  })

  return(pvals)
}

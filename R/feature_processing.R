#' Wrapper to extract features (nodes) from given pathways.
#'
#' Function to extract the features (nodes) from a given list of pathways. These
#' pathways have to be compiled with the
#' \code{\link[graphite:pathways]{pathways}} function. Features can only be
#' extracted for \'proteins\' or \'metabolites\'. Features will by default be
#' mapped to gene symbols.
#'
#' @param pathway A pathway created with \code{\link[graphite:pathways]{pathways}} command.
#' @param which Mode to extract the features, either \'proteins\' or
#'   \'metabolites\'.
#' @param org String specifying the organism, which is necessary for featureID
#'   mapping. Default: human
#' @param returntype String that specifies the returning ID type. Default:
#'   SYMBOL Options (genes/proteins): SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ
#'   Options (metabolites): HMDB, CAS, DTXCID, DTXSID, SID, CID, ChEBI, KEGG, Drugbank
#'
#' @return Feature list with gene symbols (genes/proteins) or CHEBI IDs
#'   (metabolites)
#'
#' @examples
#' pathways <- graphite::pathways("hsapiens", "kegg")[[1]]
#' getFeatures(pathways)
#' \donttest{
#' pathways <- graphite::pathways("mmusculus", "kegg")[[1]]
#' getFeatures(pathways, which = "metabolites", org = "mmusculus", returntype = "HMDB")
#'
#' pathways <- graphite::pathways("mmusculus", "kegg")
#' getFeatures(pathways, which = "proteins", org = "mmusculus", returntype = "SYMBOL")
#' }
#'
#' @importFrom graphite nodes
#'
#' @export
getFeatures <- function(pathway, which = "proteins", org = "hsapiens", returntype = "SYMBOL") {

    if (which != "proteins" && which != "metabolites") {
        stop("Only 'proteins' and 'metabolites' are supported for 'which'.",
            call. = FALSE
        )
    }

    org <- tolower(org)

    ## check for the correct organism
    if (!(org %in% getOrganisms())) {
        stop("Please insert a correct organism name! Use getOrganisms()
              to see all supported organisms.",
            call. = FALSE
        )
    }

    ## extract the features (genes/proteins or metabolites) from the pathway nodes.
    features <- graphite::nodes(pathway, which = which)

    if (length(features) == 0) {
        return(list())
    }

    if (which == "proteins") {

        ## extract the keytype of the ID format
        kt <- gsub(":.*", "", features[1])
        mapped <- gsub("[A-Z]+:", "", features)

        ## get the mapping from the keytype to the user-specific type
        mapped <- getGeneMapping(
            features = mapped, keytype = kt,
            org = org, returntype = returntype
        )
    }

    ## its special for metabolites, because sometimes there are different
    ## identifiers used in the same pathway
    if (which == "metabolites") {
        chebi <- mapIDType(
            features = features, keytype = "CHEBI",
            maptype = "ChEBI", returntype = returntype
        )

        kegg <- mapIDType(
            features = features, keytype = "KEGGCOMP",
            maptype = "KEGG", returntype = returntype
        )

        pubchem <- mapIDType(
            features = features, keytype = "PUBCHEM",
            maptype = "CID", returntype = returntype
        )

        cas <- mapIDType(
          features = features, keytype = "CAS",
          maptype = "CAS", returntype = returntype
        )
        
        hmdb <- mapIDType(
          features = features, keytype = "HMDB",
          maptype = "HMDB", returntype = returntype
        )

        mapped <- c(chebi, kegg, pubchem, cas, hmdb)
    }

    return(mapped)
}



#' Mapping between pathway encoded genes/proteins and different ID formats.
#'
#' Function to retrieve the gene/protein identifier mapping. Ongoing from
#' genes/proteins retrieved from pathway definitions, which often include two or
#' more ID formats or a format that is not present in your omics measurement,
#' this function maps those IDs to a given format. Depending on the organism,
#' additional packages have to be installed.
#'
#' @param features List of identifiers to be mapped.
#' @param keytype String specifying the ID type, e.g., "ENTREZID" or "UNIPROT".
#' @param org String that defines the organism. Default: hsapiens
#'   Options: see \code{\link[multiGSEA]{getOrganisms}}
#' @param returntype String that specifies the returning ID type. Default:
#'   SYMBOL, Options: SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ
#'
#' @return List containing mapped gene/protein IDs.
#'
#' @examples
#' features <- graphite::nodes(graphite::pathways("hsapiens", "kegg")[[1]])
#' features <- gsub("ENTREZID:", "", features)
#' keytype <- "ENTREZID"
#' getGeneMapping(features, keytype)
#'
#' getGeneMapping(features, keytype, returntype = "UNIPROT")
#' \donttest{
#' features <- graphite::nodes(graphite::pathways("rnorvegicus", "reactome")[[1]])
#' features <- gsub("UNIPROT:", "", features)
#' getGeneMapping(features, keytype = "UNIPROT", org = "rnorvegicus")
#'
#' getGeneMapping(features,
#'     keytype = "UNIPROT",
#'     org = "rnorvegicus",
#'     returntype = "ENSEMBL"
#' )
#' }
#'
#' @importFrom AnnotationDbi select
#'
#' @export
getGeneMapping <- function(features, keytype, org = "hsapiens", returntype = "SYMBOL") {

    org <- tolower(org)

    ## check for the correct organism
    if (!(org %in% getOrganisms())) {
        stop("Please insert a correct organism name! Use getOrganisms()
              to see all supported organisms.",
            call. = FALSE
        )
    }

    db <- getIDMappingDatabase(org)

    supportedIDs <- c("SYMBOL", "ENTREZID", "UNIPROT", "ENSEMBL", "REFSEQ")
    if (org != "dmelanogaster" && !returntype %in% supportedIDs) {
        stop("Insert one of the following IDs to be returned (returntype):
              SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ.",
            call. = FALSE
        )
    }

    supportedIDs <- c(supportedIDs, "FLYBASE", "FLYBASECG")
    if (org == "dmelanogaster" && !returntype %in% supportedIDs) {
        stop("Insert one of the following IDs to be returned (returntype):
              SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ, FLYBASE, FLYBASECG.",
            call. = FALSE
        )
    }

    ## design the columns field such that we create a triple ID mapping between
    ## ENTREZIDs, UNIPROT, and gene symbols
    if (keytype == "UNIPROT") {
        col <- unique(c("SYMBOL", "ENTREZID", returntype))
    }
    else {
        col <- unique(c("SYMBOL", "UNIPROT", returntype))
    }

    ## run the actual mapping of IDS and return a list of the user-given type
    map <- tryCatch(
        {
            map <- AnnotationDbi::select(db, keys = features,
	                                 columns = col, keytype = keytype)
            m <- match(unique(map[[keytype]]), map[[keytype]])
            map <- map[m, ]
            map[[returntype]][!is.na(map[[returntype]])]
        },
        error = function(cond) {
            return(list())
        }
    )

    return(map)

}


#' Mapping between pathway encoded metabolites and different metabolite ID
#' formats.
#'
#' Function to retrieve the metabolite identifier mapping. Ongoing from
#' metabolites retrieved from pathway definitions, which often include two or
#' more ID formats, this function maps those IDs to a given format. The complete
#' mapping table based on \href{https://comptox.epa.gov/dashboard}{Comptox
#' Dashboard}, \href{https://pubchem.ncbi.nlm.nih.gov/}{PubChem},
#' \href{https://hmdb.ca/}{HMDB}, and \href{https://www.ebi.ac.uk/chebi}{ChEBI}
#' is provided in the AnnotationHub package metaboliteIDmapping.
#'
#' @param features List of identifiers to be mapped.
#' @param keytype String specifying the ID type, e.g., "ChEBI" or "KEGGCOMP".
#'
#' @param returntype String that specifies the returning ID type.
#'        Default: HMDB
#'        Options: HMDB, CAS, DTXCID, DTXSID, SID, CID, ChEBI, KEGG, Drugbank
#'
#' @return List containing mapped gene/protein IDs.
#'
#' @examples
#' features <- graphite::nodes(graphite::pathways("hsapiens", "kegg")[[1]], which = "metabolites")
#' features <- gsub("KEGGCOMP:", "", features)
#' keytype <- "KEGG"
#'
#' getMetaboliteMapping(features, keytype)
#'
#' getMetaboliteMapping(features, keytype = "KEGG", returntype = "CID")
#'
#' @importFrom dplyr pull filter distinct
#' @importFrom metaboliteIDmapping metabolitesMapping
#' @importFrom magrittr %>%
#'
#' @export
getMetaboliteMapping <- function( features, keytype, returntype = "HMDB") {

  ## check for the correct metabolite mapping format
  supportedIDs <- c("HMDB", "ChEBI", "KEGG", "CAS", "DTXCID",
                    "DTXSID", "SID", "CID", "Drugbank")
  if (!returntype %in% supportedIDs) {
    stop( "Insert one of the following IDs to be returned (returntype):
              HMDB, CAS, ChEBI, KEGG, SID, CID, DTXCID, DTXSID, Drugbank, Name",
          call. = FALSE
    )
  }

  ## load the mapping table which is deposited in the
  ## metaboliteIDmapping package.
  if (!requireNamespace( "metaboliteIDmapping", quietly = TRUE)) {
    stop( "The necessary package metaboliteIDmapping is not installed.",
          call. = FALSE
          )
  }


  ## run the actual mapping of IDS and return a list of the user-given type
  map <- tryCatch(
    {

      ## to speed up the mapping, we need to subest the whole
      ## metabolitesIDmapping table in the first place to contain
      ## only thoses entries that match the given feature list
      SUBmappingTable <- metaboliteIDmapping::metabolitesMapping %>%
        dplyr::select( !!as.name( keytype), !!as.name( returntype)) %>%
        dplyr::filter( !!as.name( keytype) %in% unique( features)) %>%
        dplyr::distinct()
      colnames( SUBmappingTable) <- c("Original", "Mapped")

      SUBmappingTable %>% dplyr::pull( Mapped)

    },
    error = function(cond) {
      return(
        rep( "NA", length( features))
      )
    }
  )

  map <- map[ !is.na(map)]
  return(map)

}



#' Collect feature mapping for user given databases and omics layer.
#'
#' The functions makes use of the graphite R package to collect pathways from
#' user specified databases. Depending on the omics layer specified, the
#' function extracts either annotated genes/proteins (for transcriptome,
#' proteome layer) or metabolites (for metabolite layer). The data structure
#' that is returned is mandatory to calculate the multi-omics pathway
#' enrichment.
#'
#' @param dbs List of databases that should be queried for pathways. Default:
#'   all available databases
#' @param layer List of omics layer that should be addressed. Default: all three
#'   layer (transcriptome, proteome, metabolome)
#' @param returnTranscriptome String specifying the returned gene ID format.
#'   Default: SYMBOL Options: SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ
#' @param returnProteome String specifying the returned protein ID format.
#'   Default: SYMBOL Options: SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ
#' @param returnMetabolome String specifying the returned metabolite ID format.
#'   Default: HMDB Options: HMDB, CAS, DTXCID, DTXSID, SID, CID, ChEBI, KEGG, Drugbank
#' @param organism String specifying the organism of interest. This has direct
#'   influence on the available pathway databases. Default: "hsapiens"
#'   Options: see \code{\link[multiGSEA]{getOrganisms}}
#' @param useLocal Boolean to use local pathway/feature descriptions. In case
#'   useLocal is set to FALSE, pathway definitions and feature extraction
#'   will be recalculated. This could take several minutes depending on the
#'   database used. Pathbank, for example, contains nearly 50000 pathway
#'   definition that have to be re-mapped. useLocal has no effect when
#'   pathway definitions are retrieved for the first time. Default: TRUE
#'
#' @return Nested list with extracted and mapped pathway features.
#'
#' @examples
#'
#' getMultiOmicsFeatures(
#'     dbs = c("kegg"),
#'     layer = c("transcriptome", "proteome"),
#'     organism = "hsapiens"
#' )
#' \donttest{
#' getMultiOmicsFeatures(
#'     dbs = c("kegg", "reactome"),
#'     layer = c("transcriptome", "metabolome"),
#'     organism = "mmusculus"
#' )
#'
#' getMultiOmicsFeatures(
#'     dbs = c("reactome"),
#'     layer = c("proteome"),
#'     organism = "rnorvegicus",
#'     returnProteome = "ENTREZID"
#' )
#' }
#' @importFrom graphite pathwayDatabases pathways
#' @importFrom magrittr %>%
#' @importFrom dplyr filter pull
#' @importFrom rlang .data
#'
#' @export
getMultiOmicsFeatures <- function(dbs = c("all"), layer = c("all"),
                                  returnTranscriptome = "SYMBOL",
                                  returnProteome = "SYMBOL",
                                  returnMetabolome = "HMDB",
                                  organism = "hsapiens",
                                  useLocal = TRUE) {

    layers <- c("all", "metabolome", "proteome", "transcriptome")
    organism <- tolower(organism)

    returnTranscriptome <- toupper(returnTranscriptome)
    returnProteome <- toupper(returnProteome)
    returnMetabolome <- toupper(returnMetabolome)

    ## check for the correct transcriptome mapping format
    supportedIDs <- c("SYMBOL", "ENTREZID", "UNIPROT", "ENSEMBL", "REFSEQ")
    if (!returnTranscriptome %in% supportedIDs) {
        stop("Insert one of the following IDs to be returned (returnTranscriptome):
              SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ",
             call. = FALSE
        )
    }

    ## check for the correct proteome mapping format
    if (!returnProteome %in% supportedIDs) {
        stop("Insert one of the following IDs to be returned (returnProteome):
              SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ",
             call. = FALSE
        )
    }

    ## check for the correct metabolite mapping format
    supportedIDs <- c("HMDB", "ChEBI", "KEGG", "CAS", "DTXCID",
                      "DTXSID", "SID", "CID", "Drugbank")
    if (!returnMetabolome %in% supportedIDs) {
        stop("Insert one of the following IDs to be returned (returnMetabolome):
              HMDB, CAS, ChEBI, KEGG, SID, CID, DTXCID, DTXSID, Drugbank",
            call. = FALSE
        )
    }

    ## check if the given organism is supported
    if (!(organism %in% getOrganisms())) {
        stop("You entered an organism that is not supported!
              Use getOrganisms() to get a list of all suported organisms.",
            call. = FALSE
        )
    }

    if (sum(tolower(layer) %in% layers) != length(layer)) {
        stop("You entered wrong input for the omics layer specification.
              Options are: all, transcriptome, proteome, metabolome, or a combination thereof.",
            call. = FALSE
        )
    }

    pDBs <- graphite::pathwayDatabases()
    dbs0 <- pDBs %>%
        dplyr::filter(.data$species == organism) %>%
        dplyr::pull(.data$database)
    databases <- c("all", as.vector(dbs0))

    if (sum(tolower(dbs) %in% databases) != length(dbs)) {
        stop( paste0( "You entered wrong input for the omics layer specification.
              Options are: ", paste( databases, collapse = ' '), " or a combination thereof."),
              call. = FALSE
        )
    }

    if ("all" %in% dbs) dbs <- as.vector(dbs0)

    if ("all" %in% layer) {
          layer <- c("transcriptome", "proteome", "metabolome")
      }

    pathways <- lapply(dbs, function(x) {
        graphite::pathways(organism, x)
    })
    names(pathways) <- dbs

    features <- list()
    if ("transcriptome" %in% layer) {
        features$transcriptome <- getMappedFeatures(
            pathways = pathways,
            organism = organism,
            returnID = returnTranscriptome,
            useLocal = useLocal
        )

        ## adapt for duplicated pathways
        names( features$transcriptome) <- rename_duplicates( names( features$transcriptome))

    }

    if ("proteome" %in% layer) {
        if ("transcriptome" %in% layer && returnProteome == returnTranscriptome) {
            features$proteome <- features$transcriptome
        } else {
            features$proteome <- getMappedFeatures(
                pathways = pathways,
                organism = organism,
                returnID = returnProteome,
                useLocal = useLocal
            )

            ## adapt for duplicated pathways
            names( features$proteome) <- rename_duplicates( names( features$proteome))

        }
    }


    if ("metabolome" %in% layer) {
        features$metabolome <- getMappedFeatures(
            pathways = pathways,
            organism = organism,
            returnID = returnMetabolome,
            which = "metabolites",
            useLocal = useLocal
        )

        ## adapt for duplicated pathways
        names( features$metabolome) <- rename_duplicates( names( features$metabolome))

    }

    return(features)

}



#' Wrapper to get feature mappings.
#'
#' Feature mappings will be used from hard disk in case they have been
#' mapped before and `useLocal` is not set to be FALSE.
#' In other cases, a feature extraction will be done and the results are
#' stored for a following occasion.
#'
#' @param pathways List of pathway definitions.
#' @param returnID String specifying the returned ID format.
#' @param organism String defining the organism of analysis.
#' @param which Mode to extract the features, either \'proteins\' or
#'   \'metabolites\'.
#' @param useLocal Boolean specifying whether or not to use the local
#'              preprocessed mapping.
#'
#' @return List of mapped features for an omics layer.
getMappedFeatures <- function(pathways, returnID = "SYMBOL", organism = "hsapiens", which = "proteins", useLocal = TRUE) {

    feat <- unlist(lapply(names(pathways), function(db) {
        ap <- archivePath(paste0(organism, "_", db, "_", returnID))

        if (file.exists(ap) && useLocal) {
            loadLocal(ap)
        } else {
            tmp <- lapply(pathways[[db]], function(p) {
                getFeatures(
                    pathway = p, org = organism,
                    which = which,
                    returntype = returnID
                )
            })

            header <- rep(paste0("(", toupper(db), ") "), length(pathways[[db]]))
            names(tmp) <- paste0(header, names(tmp))
            saveRDS(tmp, file = ap)
            tmp
        }
    }), recursive = FALSE)

    return(feat)

}



#' Helper function to map only a subset of metabolite IDs
#'
#' This helper function becomes necessary since there are sometimes multiple ID
#' formats used in a single pathway definition.
#'
#' @param features List of metabolite feature IDs of the pathway.
#' @param keytype String specifying the ID format in pathway definition.
#' @param maptype String specifying the corresponding ID format in multiGSEA.
#' @param returntype String identifying the ID type that should be mapped.
#'
#' @return List of mapped metabolite IDs.
mapIDType <- function(features, keytype = "CHEBI", maptype = "ChEBI", returntype = "HMDB") {

    mapped <- c()
    ids <- gsub(paste0(keytype, ":"), "", features[grep(keytype, features)])

    if (returntype != maptype) {
        mapped <- getMetaboliteMapping(
            features = ids,
            keytype = maptype,
            returntype = returntype
        )
    } else {
        mapped <- ids
    }

    return(mapped)

}



#' Get list of supported organisms
#'
#' Get a list of organisms that are covered in our workflow through a supporting
#' `AnnotationDBi` package. Without such a package we would not be able to map
#' transcript and protein identifier between different formats. All the
#' organisms that are listed here have at lest kegg and or reactome pathway
#' annotation that can be queried by `graphite`.
#'
#' @return List of supported organisms
#'
#' @examples
#' getOrganisms()
#' @export
getOrganisms <- function() {

    orglist <- c(
        "hsapiens", "rnorvegicus", "mmusculus", "sscrofa",
        "btaurus", "celegans", "dmelanogaster", "drerio",
        "ggallus", "xlaevis", "cfamiliaris"
    )

    return(orglist)

}



#' Get the correct ID mapping database
#'
#' Check by means of the given organism name if the required `AnnotationDbi`
#' package is installed. Select the ID mapping table based on the organism name
#' and return it.
#'
#' @param organism String that defines the organism.
#'
#' @return AnnotationDbi database for ID mapping.
getIDMappingDatabase <- function(organism) {

    map <- c(
        hsapiens = "org.Hs.eg.db", rnorvegicus = "org.Rn.eg.db",
        mmusculus = "org.Mm.eg.db", sscrofa = "org.Ss.eg.db",
        btaurus = "org.Bt.eg.db", celegans = "org.Ce.eg.db",
        dmelanogaster = "org.Dm.eg.db", drerio = "org.Dr.eg.db",
        ggallus = "org.Gg.eg.db", xlaevis = "org.Xl.eg.db",
        cfamiliaris = "org.Cf.eg.db"
    )

    stopifnot(organism %in% names(map))
    pkg <- map[[organism]]

    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(paste0("The necessary package ", pkg, " is not installed."),
            call. = FALSE
        )
    }

    return(get(pkg, envir = getNamespace(pkg)))

}



#' Pre-rank features prior to calculating enrichment scores.
#'
#' Rank features based on the direction of their fold change and their magnitude
#' implicated through their assigned p-value.
#'
#' @param logFC Vector containing the log-transformed fold changes of features.
#' @param pvalues Vector containing the p-values associated with those logFCs.
#' @param base Integer specifying the base of the logarithm. Default: 10
#'
#' @return Vector of pre-ranked features, still unsorted
#'
#' @examples
#' logFC <- rnorm(10)
#' pvalues <- runif(10)
#' rankFeatures(logFC, pvalues)
#' @export
rankFeatures <- function(logFC, pvalues, base = 10) {

    return(sign(logFC) * -log(pvalues, base = base))

}



#' Create an empty data structure for measured omics features
#'
#' This function creates a data structure of nested but empty lists. One list
#' for each omics layer. By default all three supported omics layer are used to
#' create a data structures with three empty sublists: transcriptome, proteome,
#' and metabolome.
#'
#' @param layer List specifying the omics layer which should be created
#'
#' @return List with length(layer) empty sublists
#'
#' @examples
#' initOmicsDataStructure()
#' initOmicsDataStructure(c("transcriptome", "proteome"))
#' initOmicsDataStructure(c("Transcriptome", "Metabolome"))
#' @export
initOmicsDataStructure <- function(layer = c("transcriptome", "proteome", "metabolome")) {

    l <- c()
    layer <- tolower(layer)
    if (length(grep("trans*", layer)) > 0) l <- c(l, "transcriptome")
    if (length(grep("prote*", layer)) > 0) l <- c(l, "proteome")
    if (length(grep("metabo*", layer)) > 0) l <- c(l, "metabolome")

    if (length(l) != length(layer)) {
        stop("Not all your omics layer could be mapped to
              'transcriptome', 'proteome', or 'metabolome'. ",
            call. = TRUE
        )
    }

    empty_structure <- rep(list(list()), length(l))
    names(empty_structure) <- l

    return(empty_structure)

}



#' Helper function to get all different metabolite ID formats
#'
#' This helper function extracts all used ID formats in all pathways 
#' and returns a nested list for each pathway database.
#'
#' @param pathways List of pathway databases and their pathway definition.
#'
#' @return List of metabolite ID formats.
getMetaboliteIDformats <- function(pathways){
  
  n1 <- lapply(names(pathways), function( dbs){
    
    n2 <- lapply(pathways[[dbs]], function(p){
      
      unique(gsub( ":.*", "", graphite::nodes(p, which = "metabolites")))
      
      })
    
    unique(unlist(n2))

  })
  
  names(n1) <- names(pathways)
  return(n1)
}

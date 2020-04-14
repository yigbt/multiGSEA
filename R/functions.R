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
#'   Options (metabolites): HMDB, CAS, DTXCID, DTXSID, CID, SID, ChEBI, KEGG
#'
#' @return Feature list with gene symbols (genes/proteins) or CHEBI IDs
#'   (metabolites)
#'
#' @examples
#' pathways <- graphite::pathways( "hsapiens", "reactome")[[1]]
#' getFeatures( pathways)
#'
#' pathways <- graphite::pathways( "mmusculus", "kegg")[[1]]
#' getFeatures( pathways, which = "metabolites", org = "mmusculus", returntype = "HMDB")
#'
#' \donttest{
#' pathways <- graphite::pathways( "mmusculus", "kegg")
#' getFeatures( pathways, which = "proteins", org = "mmusculus", returntype = "SYMBOL")
#' }
#'
#' @importFrom graphite nodes
#'
#' @export
getFeatures <- function( pathway, which = "proteins", org = "hsapiens", returntype = "SYMBOL"){
    
    if( which != "proteins" & which != "metabolites"){
        stop( "Only 'proteins' and 'metabolites' are supported for 'which'.",
              call. = FALSE)
    }
    
    org <- tolower( org)
    
    ## check for the correct organism 
    if( !( org %in% getOrganisms() )){
        stop( "Please insert a correct organism name! Use getOrganisms() to see all supported organisms.", 
              call. = FALSE)
    }
    
    ## extract the features (genes/proteins or metabolites) from the pathway nodes.
    features <- nodes( pathway, which = which)
    
    if( length( features) == 0) return( list())
    
    if( which == "proteins"){

        # Remove indentifier strings and keep only the IDs
        if( startsWith( features[1], "ENTREZ")){ 
            mapped <- gsub( "ENTREZID:", "", features)
            mapped <- getGeneMapping( features = mapped, keytype = "ENTREZID",
                                      org = org, returntype = returntype)
        }
        if( startsWith( features[1], "UNIPROT")){
            mapped <- gsub( "UNIPROT:", "", features)
            mapped <- getGeneMapping( features = mapped, keytype = "UNIPROT",
                                      org = org, returntype = returntype) 
        }
    }
    
    ## its special for metabolites, because sometimes there are different
    ## identifiers used in the same pathway
    if( which == "metabolites"){

        chebi <- mapIDType( features = features, keytype = "CHEBI",
                            maptype = "ChEBI", returntype = returntype)

        kegg <- mapIDType( features = features, keytype = "KEGGCOMP",
                           maptype = "KEGG", returntype = returntype)

        pubchem <- mapIDType( features = features, keytype = "PUBCHEM",
                              maptype = "CID", returntype = returntype)
    }
    
    return( c( chebi, kegg, pubchem))
    
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
mapIDType <- function( features, keytype = "CHEBI", maptype = "ChEBI", returntype = "HMDB"){
    
    mapped <- c()
    ids <- gsub( paste0( keytype, ":"), "", features[ grep( keytype, features)])
    if( returntype != maptype){
        mapped <- c( mapped, getMetaboliteMapping( features = ids,
                                                   keytype = maptype,
                                                   returntype = returntype))
    }else{
        mapped <- c( mapped, ids)
    }
    
    return( mapped)
    
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
#' @param org String that defines the organism. Default: hsapiens Options: hsapiens,
#'   rnorvegicus, mmusculus
#' @param returntype String that specifies the returning ID type. Default:
#'   SYMBOL, Options: SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ
#'
#' @return List containing mapped gene/protein IDs.
#'
#' @examples
#' features <- graphite::nodes( graphite::pathways( "hsapiens", "kegg")[[1]])
#' features <- gsub( "ENTREZID:", "", features)
#' keytype <- "ENTREZID"
#' getGeneMapping( features, keytype)
#'
#' getGeneMapping( features, keytype, returntype = "UNIPROT")
#'
#' features <- graphite::nodes( graphite::pathways( "rnorvegicus", "reactome")[[1]])
#' features <- gsub( "UNIPROT:", "", features)
#' getGeneMapping( features, keytype = "UNIPROT", org = "rnorvegicus")
#'
#' getGeneMapping( features, keytype = "UNIPROT",
#'                 org = "rnorvegicus",
#'                 returntype = "ENSEMBL")
#'
#' @importFrom stringr str_detect regex
#' @importFrom AnnotationDbi select
#'
#' @export
getGeneMapping <- function( features, keytype, org = "hsapiens", returntype = "SYMBOL"){
    
    org <- tolower( org)
    
    ## check for the correct organism 
    if( !( org %in% getOrganisms() )){
        stop( "Please insert a correct organism name! Use getOrganisms() to see all supported organisms.", 
              call. = FALSE)
    }
    
    db <- getIDMappingDatabase( org)
    
    if( !str_detect( returntype, regex( "^SYMBOL$|^ENTREZID$|^UNIPROT$|^ENSEMBL$|^REFSEQ$"))){
        stop( "Insert one of the following ID types to be returned (returntype): SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ", call. = FALSE)
    }
    
    ## design the columns field such that we create a triple ID mapping between
    ## ENTREZIDs, UNIPROT, and gene symbols
    if( keytype == "UNIPROT"){ col <- unique( c( "SYMBOL", "ENTREZID", returntype))}
    else{ col <- unique( c( "SYMBOL", "UNIPROT", returntype))} 
    
    ## run the actual mapping of IDS and return a list of the user-given type
    map <- tryCatch({
        map <- select( db, keys = features, columns = col, keytype = keytype)
        m <- match( unique( map[[ keytype]]), map[[ keytype]])
        map <- map[ m, ]
        map[[returntype]][ !is.na( map[[returntype]])]
    }, 
    error=function( cond){
        return(list())
    })
    
    return( map)
    
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
#' 
#' @export
getOrganisms <- function(){
    
    orglist <- c( "hsapiens", "rnorvegicus", "mmusculus", "sscrofa",
                  "btaurus", "celegans", "dmelanogaster", "drerio",
                  "ggallus", "xlaevis", "cfamiliaris")
    
    return( orglist)
    
}


#' Get the correct ID mapping database
#' 
#' Check by means of the given organism name if the required `AnnotationDbi`
#' package is installed. Select the ID mapping table based on the organism name
#' and return it.
#' 
#' @param organism String that defines the organism.
#' 
#' @return AnnotationDbi databse for ID mapping.
getIDMappingDatabase <- function( organism){
    
    # check for the installed package to do the human ID mapping
    if( organism == "hsapiens") pkg <- "org.Hs.eg.db"
    
    ## check for the installed package to do the rat ID mapping
    if( organism == "rnorvegicus") pkg <- "org.Rn.eg.db"
    
    ## check for the installed package to do the mouse ID mapping
    if( organism == "mmusculus") pkg <- "org.Mm.eg.db"
    
    ## check for the installed package to do the pig ID mapping
    if( organism == "sscrofa") pkg <- "org.Ss.eg.db"
    
    ## check for the installed package to do the bovine ID mapping
    if( organism == "btaurus") pkg <- "org.Bt.eg.db"
    
    ## check for the installed package to do the _C.elegans_ ID mapping
    if( organism == "celegans") pkg <- "org.Ce.eg.db"
    
    ## check for the installed package to do the fruit fly ID mapping
    if( organism == "dmelanogaster") pkg <- "org.Dm.eg.db"
    
    ## check for the installed package to do the zebrafish ID mapping
    if( organism == "drerio") pkg <- "org.Dr.eg.db"
    
    ## check for the installed package to do the chicken ID mapping
    if( organism == "ggallus") pkg <- "org.Gg.eg.db"
    
    ## check for the installed package to do the frog ID mapping
    if( organism == "xlaevis") pkg <- "org.Xl.eg.db"
    
    ## check for the installed package to do the dog ID mapping
    if( organism == "cfamiliaris") pkg <- "org.Cf.eg.db"
    
    if( !requireNamespace( pkg, quietly = TRUE)){
        stop( paste0("The necessary package ", pkg, " is not installed."),
              call. = FALSE)
    }
    return( eval( parse(text=paste0(pkg, "::", pkg))))
    
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
#' is provided with the package.
#'
#' @param features List of identifiers to be mapped.
#' @param keytype String specifying the ID type, e.g., "ChEBI" or "KEGGCOMP".
#'
#' @param returntype String that specifies the returning ID type. 
#'        Default: HMDB
#'        Options: HMDB, CAS, DTXCID, DTXSID, CID, SID, ChEBI, KEGG
#'
#' @return List containing mapped gene/protein IDs.
#'
#' @examples
#' features <- graphite::nodes( graphite::pathways( "hsapiens", "kegg")[[1]], which = "metabolites")
#' features <- gsub( "KEGGCOMP:", "", features)
#' keytype <- "KEGG"
#'
#' getMetaboliteMapping( features, keytype)
#'
#' getMetaboliteMapping( features, keytype = "KEGG", returntype = "CID")
#'
#' @importFrom stringr str_detect regex
#'
#' @export
getMetaboliteMapping <- function( features, keytype, returntype = "HMDB"){
    
    ## check for the correct metabolite mapping format
    match <- "^HMDB$|^ChEBI|^KEGG$|^CAS$|^DTXCID$|^DTXSID$|^CID$|^SID$"
    if( !str_detect( returntype, regex( match) )){
        stop( "Insert one of the following IDs to be returned (returntype):
              HMDB, CAS, ChEBI, KEGG, CID, SID, DTXCID, DTXSID", 
              call. = FALSE)
    }
    
    ## run the actual mapping of IDS and return a list of the user-given type
    map <- tryCatch({
        m <- match( unique( features), metabolitesMapping[[ keytype]])
        m <- m[ !is.na( m)]
        metabolitesMapping[[returntype]][ m]
    }, 
    error=function( cond){
        return(list())
    })
    
    map <- map[ !is.na( map)]
    
    if( returntype == "DTXCID") map <- paste0( "DTXCID", map)
    if( returntype == "DTXSID") map <- paste0( "DTXSID", map)
    
    return( map)
    
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
#'   Default: HMDB Options: HMDB, CAS, DTXCID, DTXSID, CID, SID, ChEBI, KEGG
#' @param organism String specifying the organism of interest. This has direct
#'   influence on the available pathway databases. Default: "hsapiens" Options:
#'   "hsapiens", "mmusculus", "rnorvegicus"
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
#' getMultiOmicsFeatures( dbs = c( "kegg", "reactome"),
#'                        layer = c( "transcriptome", "metabolome"),
#'                        organism = "mmusculus")
#'
#' getMultiOmicsFeatures( dbs = c("reactome"),
#'                        layer = c( "proteome"),
#'                        organism = "rnorvegicus",
#'                        returnProteome = "ENTREZID")
#'
#' @importFrom graphite pathwayDatabases pathways
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect regex
#' @importFrom dplyr filter pull
#' @importFrom rlang .data
#'
#' @export
getMultiOmicsFeatures <- function( dbs = c("all"), layer = c("all"), 
                                   returnTranscriptome = "SYMBOL", 
                                   returnProteome = "SYMBOL", 
                                   returnMetabolome = "HMDB", 
                                   organism = "hsapiens", 
                                   useLocal = TRUE){
    
    layers = c( "all", "metabolome", "proteome", "transcriptome")
    organism <- tolower( organism)
    
    returnTranscriptome <- toupper( returnTranscriptome)
    returnProteome <- toupper( returnProteome)
    returnMetabolome <- toupper( returnMetabolome)
    
    match <- "^SYMBOL$|^ENTREZID$|^UNIPROT$|^ENSEMBL$|^REFSEQ$"
    ## check for the correct transcriptome mapping format
    if( !str_detect( returnTranscriptome,  regex( match) )){
        stop( "Insert one of the following IDs to be returned (returnTranscriptome):
              SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ",
              call. = FALSE)
    }
    
    ## check for the correct proteome mapping format
    if( !str_detect( returnTranscriptome,  regex( match) )){
        stop( "Insert one of the following IDs to be returned (returnProteome):
              SYMBOL, ENTREZID, UNIPROT, ENSEMBL, REFSEQ",
              call. = FALSE)
    }
    
    match <- "^HMDB$|^ChEBI|^KEGG$|^CAS$|^DTXCID$|^DTXSID$|^CID$|^SID$"
    ## check for the correct metabolite mapping format
    if( !str_detect( returnMetabolome,  regex( match) )){
        stop( "Insert one of the following IDs to be returned (returnMetabolome):
              HMDB, CAS, ChEBI, KEGG, CID, SID, DTXCID, DTXSID",
              call. = FALSE)
    }
    
    ## check if the given organism is supported
    if( !( organism %in% getOrganisms())){
        stop( "You entered an organism that is not supported! 
              Use getOrganisms() to get a list of all suported organisms.", 
              call. = FALSE)
    }
    
    if( sum( tolower(layer) %in% layers) != length( layer)){
        stop( "You entered wrong input for the omics layer specification.
              Options are: all, transcriptome, proteome, metabolome, or a combination thereof.",
              call. = FALSE)
    }
    
    pDBs <- pathwayDatabases()
    databases <- c( "all", as.vector( pDBs %>% filter( .data$species == organism) %>% pull( .data$database)))
    
    if( sum( tolower(dbs) %in% databases) != length( dbs)){
        stop( "You entered wrong input for the omics layer specification.
              Options are: all, kegg, pathbank, reactome, biocarta, nci, panther, 
              smpdb, pharmgkb, or a combination thereof.",
              call. = FALSE)
    }
    
    if( "all" %in% dbs){
        dbs <- as.vector( pDBs %>% filter( .data$species == organism) %>% pull( .data$database))
    }
    
    if( "all" %in% layer){
        layer <- c( "transcriptome", "proteome", "metabolome")
    }
    
    pathways <- lapply( dbs, function(x) { pathways( organism, x)})
    names( pathways) <- dbs
    
    features <- list()
    if( "transcriptome" %in% layer){
        
        trans <- unlist( lapply( names(pathways), function(dbname){
            
            ap <- archivePath( paste0( organism, "_", dbname, "_", returnTranscriptome))
            
            if( file.exists( ap) && useLocal){
                
                loadLocal( ap)
                
            } else{
                
                tmp <- lapply( pathways[[dbname]], function(p){ 
                    getFeatures( pathway = p, org = organism, 
                                 which = "proteins", 
                                 returntype = returnTranscriptome)
                })
                
                names( tmp) <- paste0( rep( paste0( "(", toupper(dbname), ") "), length( pathways[[dbname]])), names( tmp))
                saveRDS( tmp, file = ap)
                tmp
                
            }
            
        }), recursive = FALSE)
        
        features$transcriptome <- trans
        
    }  
    
    if( "proteome" %in% layer){
        
        if( returnProteome == returnTranscriptome){
            features$proteome <- trans
        } else{
            
            prot <- unlist( lapply( names(pathways), function(dbname){
                
                ap <- archivePath( paste0( organism, "_", dbname, "_", returnProteome))
                
                if( file.exists( ap) && useLocal){
                    
                    loadLocal( ap)
                    
                } else{
                    
                    tmp <- lapply( pathways[[dbname]], function(p){
                        getFeatures( pathway = p, org = organism, 
                                     which = "proteins",
                                     returntype = returnProteome)
                    })
                    
                    names( tmp) <- paste0( rep( paste0( "(", toupper(dbname), ") "), length( pathways[[dbname]])), names( tmp))
                    saveRDS( tmp, file = ap)
                    tmp
                    
                }
                
            }), recursive = FALSE)
            
            features$proteome <- prot
            
        }
        
    }
    
    
    if( "metabolome" %in% layer){
        
        meta <- unlist( lapply( names(pathways), function(dbname){
            
            ap <- archivePath( paste0( organism, "_", dbname, "_", returnMetabolome))
            
            if( file.exists( ap) && useLocal){
                
                loadLocal( ap)
                
            } else{
                
                tmp <- lapply( pathways[[dbname]], function(p){
                    getFeatures( pathway = p, org = organism, 
                                 which = "metabolites",
                                 returntype = returnMetabolome)
                })
                
                names( tmp) <- paste0( rep( paste0( "(", toupper(dbname), ") "), length( pathways[[dbname]])), names( tmp))
                saveRDS( tmp, file = ap)
                tmp
                
            }
            
        }), recursive = FALSE)
        
        features$metabolome <- meta
        
    }
    
    return( features)
    
}




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
#'
#' @return Nested list containing the enrichment scores for each given pathway
#'   and omics layer.
#'
#' @examples
#'
#' # Download pathway definition and extract features
#' pathways <- getMultiOmicsFeatures( dbs = c( "kegg"), layer = c("transcriptome", "proteome"))
#'
#' # load omics data and calculate ranks
#' data(transcriptome)
#' data(proteome)
#' ranks <- initOmicsDataStructure( c("transcriptome", "proteome"))
#' ranks$transcriptome <- rankFeatures( transcriptome$logFC, transcriptome$pValue)
#' names( ranks$transcriptome) <- transcriptome$Symbol
#' ranks$proteome <- rankFeatures( proteome$logFC, proteome$pValue)
#' names( ranks$proteome) <- proteome$Symbol
#' 
#' ## run the enrichment
#' multiGSEA( pathways, ranks)
#'
#'
#' @importFrom fgsea fgsea
#'
#' @export
multiGSEA <- function(pathways, ranks) {
    
    # Go through all omics layer.
    es <- lapply(names(pathways), function(omics) {
        fgsea(pathways[[omics]], ranks[[omics]], nperm = 1000, minSize = 5)
        
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
#' pathways <- getMultiOmicsFeatures( dbs = c( "kegg"), layer = c("transcriptome", "proteome"))
#'
#' # loamultiple of 4d omics data and calculate ranks
#' data(transcriptome)
#' data(proteome)
#' ranks <- initOmicsDataStructure( c("transcriptome", "proteome"))
#' ranks$transcriptome <- rankFeatures( transcriptome$logFC, transcriptome$pValue)
#' names( ranks$transcriptome) <- transcriptome$Symbol
#' ranks$proteome <- rankFeatures( proteome$logFC, proteome$pValue)
#' names( ranks$proteome) <- proteome$Symbol
#' 
#' # run the enrichment
#' es <- multiGSEA( pathways, ranks)
#' 
#' extractPvalues( enrichmentScores = es,
#'                 pathwayNames = names( pathways[[1]]))
#'
#' @export
extractPvalues <- function( enrichmentScores, pathwayNames){
    
    # Go through all the pathways
    res <- lapply( pathwayNames, function( name){ 
        
        # Go through all the possible omics layer
        unlist(lapply( names( enrichmentScores), function( y){
            
            df <- enrichmentScores[[y]][which( enrichmentScores[[y]]$pathway== name),c(2,3)]
            if( nrow( df) == 0){ df <- data.frame( pval = NA, padj = NA)}
            names(df) <- paste0( y, "_", names(df))
            df
            
        }))
        
    })
    
    # Combine the list elements to a data frame and assign the pathway names as rownames
    res <- data.frame( do.call( rbind, res))
    
    return( res)
    
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
#' initOmicsDataStructure( c("transcriptome", "proteome"))
#' initOmicsDataStructure( c("Transcriptome", "Metabolome"))
#'
#' @export
initOmicsDataStructure <- function( layer = c("transcriptome", "proteome", "metabolome")){
    
    l <- c()
    if( length( grep( "TRANS*", layer, ignore.case = TRUE)) > 0) l <- c( l, "transcriptome")
    if( length( grep( "PROTE*", layer, ignore.case = TRUE)) > 0) l <- c( l, "proteome")
    if( length( grep( "METABO*", layer, ignore.case = TRUE)) > 0) l <- c( l, "metabolome")
    
    if( length( l) != length( layer)){
        stop( "Not all your omics layer could be mapped to 'transcriptome', 'proteome', or 'metabolome'. ",
              call. = TRUE)
    }
    
    empty_structure <- rep( list( list()), length( l))
    names( empty_structure) <- l
    
    return( empty_structure)
    
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
#'   Default: "stouffer"no visible binding for global variable Options:
#'   "stouffer", "fisher", "edgington"
#' @param weights List of weights that will be used in a weighted Stouffer
#'   method.
#'
#' @return Vector of length \code{nrow(df)} with combined p-values.
#'
#' @examples
#' df <- cbind( runif(5), runif(5), runif(5))
#' colnames( df) <- c("trans.pval", "prot.pval", "meta.pval")
#'
#' # run the unweighted summation of z values
#' combinePvalues( df)
#'
#' # run the weighted variant
#' combinePvalues( df, weights = c( 10,5,1))
#'
#' # run the Fisher's combined probability test
#' combinePvalues( df, method = "fisher")
#'
#' # run the Edgington's method
#' combinePvalues( df, method = "edgington")
#'
#' @importFrom metap sumz sumlog sump
#'
#' @export
combinePvalues <- function( df, method = "stouffer", weights = NULL){
    
    method <- tolower( method)
    if( method != "stouffer" & method != "fisher" & method != "edgington"){
        stop( "You can chose between the 'stouffer', 'edgington', 
              and 'fisher' method to combine p-values.",
              call. = FALSE)
    }
    
    cols <- grep( "pval", colnames( df))
    
    pvals <- apply( df, 1, function( row){
        
        row <- row[ cols]
        row <- row[ !is.na( row)]
        
        if( length( row) >= 2) {
            if( method == "fisher"){
                p <- sumlog( row)
                p$p
            }else if( method == "edgington"){
                p <- sump( row)
                p$p
            }else {
                
                ## sumz allows only p-values smaller than 1
                row <- row[ row > 0 & row < 1]
                
                if( length( row) >= 2){
                    if( length( weights) > 0){
                        p <- sumz( row, weights = weights)
                    }else{
                        p <- sumz( row)
                    }
                    p$p
                }else if( length( row == 1)){
                    row[1]
                }else{
                    NA
                }
                
            }
        } else if( length( row) == 1){
            row[1]
        } else{
            NA
        }
    })
    
    return( pvals)
    
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
#' rankFeatures( logFC, pvalues)
#'
#' @export
rankFeatures <- function( logFC, pvalues, base = 10){
    
    return( sign(logFC) * -log( pvalues, base = base))
    
}



#' Retrieve path to a cached file.
#'
#' The function retrieves the path to a file that is cahed in the archive
#' directory.
#'
#' @param filename Name of the file.
#'
#' @return String containing the path to the file.
archivePath <- function( filename){
    
    ad <- archiveDir()
    return( paste0( ad, "/", filename, ".rds"))
    
}





#' Retrieve the path to the cache directory.
#'
#' Retrieve the path to the cache directory for the multiGSEA package. Create
#' the cache directory if need be.
#'
#' @return String containing the path to the cache directory.
#'
#' @importFrom rappdirs user_cache_dir
archiveDir <- function(){
    
    ad <- user_cache_dir( "multiGSEA")
    
    if ( !file.exists( ad)) {
        if ( !dir.create( ad, FALSE, TRUE)) 
            stop( "An error occured during creating the archive directory: ", ad, 
                  call. = FALSE)
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
#' @importFrom magrittr %>%
loadLocal <- function( filename){
    
    res <- try( readRDS( filename), silent = TRUE)
    
    if( "try-error" %in% is(res)){
        return( NULL)
    } else{
        return( res)
    }
    
}


# create binding for global variable of our mapping data frame
if(getRversion() >= "2.15.1")  globalVariables( c("metabolitesMapping"))

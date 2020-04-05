# The `multiGSEA` `R` package

- The [multiGSEA](https://github.com/yigbt/multiGSEA) package was
  designed to run a robust GSEA-based pathway enrichmemt for multiple
  omics layers.
- Pathway definition can be downloaded from up to eight different
  pathway databases by means of the `graphite` Bioconductor package.
- Feature mapping for transcripts and proteins towards Entrez Gene
  IDs, Uniprot, Gene Symbol, RefSeq, Ensembl
- Metabolite ID mapping to Comptox Dashboard IDs (DTXCID, DTXSID),
  CAS-numbers, Pubchem IDs (CID, SID), HMDB, KEGG, and ChEBI
- Enrichment scores and p-values are shown for each layer separately
  and aggregated
  
# Install

You can install the most up to date version easily with
[devtools](https://github.com/hadley/devtools):

    install.packages("devtools")
    devtools::install_github("https://github.com/yigbt/multiGSEA")


Once installed, just load the *multiGSEA* package with:

    library(multiGSEA)




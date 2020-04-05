# The `multiGSEA` `R` package

The `multiGSEA` package was designed to run a robust GSEA-based pathway
enrichment for multiple omics layers. The enrichment is calculated for each
omics layer separately and aggregated p-values are calculated afterwards to
derive a composite multi-omics pathway enrichment.

Pathway definitions can be downloaded from up to eight different
pathway databases by means of the
[`graphite`](http://bioconductor.org/packages/release/bioc/html/graphite.html)
Bioconductor package.

Features of the transcriptome and proteome level can be mapped to the following ID formats:

	* Entrez Gene ID
	* Uniprot IDs
	* Gene Symbols
	* RefSeq
	* Ensembl
	
Features of the metabolome layer can be mapped to:

	* Comptox Dashboard IDs (DTXCID, DTXSID)
	* CAS-RN numbers
	* Pubchem IDs (CID)
	* HMDB IDs
	* KEGG IDs
	* ChEBI IDs

  
# Install

You can install the most up to date version easily with
[devtools](https://github.com/hadley/devtools):

    install.packages("devtools")
    devtools::install_github("https://github.com/yigbt/multiGSEA")


Once installed, just load the *multiGSEA* package with:

    library(multiGSEA)


# Workflow

A common workflow is exemplified in the package vignette and is typically separated
in the following steps:

1. Load required libraries, including the `multiGSEA` package, and omics data sets.
2. Create data structure for enrichment analysis.
3. Download and customize the pathway definitions.
4. Run the pathway enrichment for each omics layer.
5. Calculate the aggregated pathway enrichment.


For more information please have a look in the vignette.

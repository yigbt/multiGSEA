<!-- badges: start -->
  [![R-CMD-check](https://github.com/yigbt/multiGSEA/actions/workflows/test.yaml/badge.svg)](https://github.com/yigbt/multiGSEA/actions/workflows/test.yaml)
  <!-- badges: end -->


# The `multiGSEA` `R` package

## Authors

Sebastian Canzler and Jörg Hackermüller

[multiGSEA: a GSEA-based pathway enrichment analysis for multi-omics data](https://doi.org/10.1186/s12859-020-03910-x), _BMC Bioinformatics_ 21, 561 (2020)

## Introduction

The `multiGSEA` package was designed to run a robust GSEA-based
pathway enrichment for multiple omics layers. The enrichment is
calculated for each omics layer separately and aggregated p-values are
calculated afterwards to derive a composite multi-omics pathway
enrichment.

Pathway definitions can be downloaded from up to eight different
pathway databases by means of the
[`graphite`](http://bioconductor.org/packages/release/bioc/html/graphite.html)
Bioconductor package.

Features of the transcriptome and proteome level can be mapped to the
following ID formats:

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
	* Drugbank IDs
	* Common names

  
Please note, that the mapping of metabolite IDs is accomplished
through the `metaboliteIDmapping` package.  This `AnnotationHub`
package provides a comprehensive mapping table with more than one
million compounds (`metaboliteIDmapping` on our [github
page](https://github.com/yigbt/metaboliteIDmapping) or at
[Bioconductor](http://bioconductor.org/packages/metaboliteIDmapping/)).
  
  
## Installation

There are two ways to install the `multiGSEA` package. For both you
have to install and start R in at least version 4.0:

(i) Use the Bioconductor framework:

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("multiGSEA")
```

(ii) Alternatively, you can install the most up to date version
(development) easily with
[devtools](https://github.com/hadley/devtools):

```R
install.packages("devtools")
devtools::install_github("https://github.com/yigbt/multiGSEA")
```

Once installed, just load the `multiGSEA` package with:

```R
library(multiGSEA)
```


# Workflow

A common workflow is exemplified in the package vignette and is
typically separated in the following steps:

1. Load required libraries, including the `multiGSEA` package, and
   omics data sets.
2. Create data structure for enrichment analysis.
3. Download and customize the pathway definitions.
4. Run the pathway enrichment for each omics layer.
5. Calculate the aggregated pathway enrichment.


For more information please have a look in the vignette at our
[Bioconductor
page](https://bioconductor.org/packages/release/bioc/vignettes/multiGSEA/inst/doc/multiGSEA.html).


# LICENSE

Copyright (C) 2011 - 2020 Helmholtz Centre for Environmental Research
UFZ.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the UFZ
License document for more details:
<https://github.com/yigbt/multiGSEA/blob/master/LICENSE.md>

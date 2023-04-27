context("Pathway enrichment")

# load transcriptomic features
data(transcriptome)
rF_trans <- rankFeatures(transcriptome$logFC, transcriptome$pValue)
names(rF_trans) <- transcriptome$Symbol

# load proteomic features
data(proteome)
rF_prot <- rankFeatures(proteome$logFC, proteome$pValue)
names(rF_prot) <- proteome$Symbol

# load metabolomic features
data(metabolome)
rF_meta <- rankFeatures(metabolome$logFC, metabolome$pValue)
names(rF_meta) <- metabolome$HMDB
names(rF_meta) <- gsub("HMDB", "HMDB00", names(rF_meta))

library(metaboliteIDmapping)

test_that("provided data can be ranked", {
  expect_is(rF_trans, "numeric")
  expect_equal(length(rF_trans), nrow(transcriptome))

  expect_is(rF_prot, "numeric")
  expect_equal(length(rF_prot), nrow(proteome))

  expect_is(rF_meta, "numeric")
  expect_equal(length(rF_meta), nrow(metabolome))
})

test_that("pathway enrichment works.", {
  dbs <- "kegg"
  df <- getMultiOmicsFeatures(dbs = dbs)

  omics_data <- initOmicsDataStructure()
  expect_is(omics_data, "list")
  expect_equal(length(omics_data), 3)
  expect_equivalent(names(omics_data), c("transcriptome", "proteome", "metabolome"))

  omics_data$transcriptome <- rF_trans
  omics_data$proteome <- rF_prot
  omics_data$metabolome <- rF_meta

  expect_warning((es <- multiGSEA(df, omics_data)))
  expect_equal(length(es), length(omics_data))
  expect_true("pathway" %in% colnames(es$transcriptome))
  expect_true("pval" %in% colnames(es$transcriptome))
  expect_true("padj" %in% colnames(es$transcriptome))

  expect_true("pathway" %in% colnames(es$proteome))
  expect_true("pval" %in% colnames(es$proteome))
  expect_true("padj" %in% colnames(es$proteome))

  expect_true("pathway" %in% colnames(es$metabolome))
  expect_true("pval" %in% colnames(es$metabolome))
  expect_true("padj" %in% colnames(es$metabolome))

  eP <- extractPvalues(
    enrichmentScores = es,
    pathwayNames = names(df[[1]])
  )

  expect_equal(ncol(eP), 2 * length(omics_data))
  expect_equal(nrow(eP), length(df$transcriptome))

  expect_lte(max(eP, na.rm = TRUE), 1)
  expect_gte(min(eP, na.rm = TRUE), 0)
})

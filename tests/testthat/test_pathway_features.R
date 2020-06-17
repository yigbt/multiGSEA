context("Pathways and feature extraction")

test_that("transcriptomic features get mapped", {
    layer <- c("transcriptome")
    dbs <- c("kegg")

    df <- getMultiOmicsFeatures(dbs = dbs, layer = layer)

    expect_is(df, "list")
    expect_equal(length(df), 1)
    expect_identical(names(df), c("transcriptome"))

    expect_match(names(df$transcriptome), "^\\(KEGG\\)")
})

test_that("proteomic features get mapped", {
    layer <- c("proteome")
    dbs <- c("kegg")

    df <- getMultiOmicsFeatures(dbs = dbs, layer = layer)

    expect_is(df, "list")
    expect_equal(length(df), 1)
    expect_identical(names(df), c("proteome"))

    expect_match(names(df$proteome), "^\\(KEGG\\)")
})

test_that("metabolomic features get mapped", {
    layer <- c("metabolome")
    dbs <- c("kegg")

    df <- getMultiOmicsFeatures(dbs = dbs, layer = layer)

    expect_is(df, "list")
    expect_equal(length(df), 1)
    expect_identical(names(df), c("metabolome"))

    expect_match(names(df$metabolome), "^\\(KEGG\\)")
})


test_that("each layer has equal length", {
    dbs <- c("kegg")
    df <- getMultiOmicsFeatures(dbs = dbs, organism = "hsapiens")

    expect_equal(length(df), 3)
    expect_equal(length(df$transcriptome), length(df$proteome))
    expect_equal(length(df$transcriptome), length(df$metabolome))
})

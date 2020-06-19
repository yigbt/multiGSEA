context("P-value aggregation")

test_that("method selection works", {
    df <- cbind(runif(5), runif(5), runif(5))
    colnames(df) <- c("trans.padj", "prot.padj", "meta.padj")

    method <- "random"
    expect_error(combinePvalues(df, method = method))
})


test_that("Fisher method works", {
    df <- cbind(runif(5), runif(5), runif(5))
    colnames(df) <- c("trans.padj", "prot.padj", "meta.padj")

    method <- "Fisher"
    expect_is(combinePvalues(df, method = method), "numeric")
    expect_equal(length(combinePvalues(df, method = method)), nrow(df))
})


test_that("Edgington method works", {
    df <- cbind(runif(5), runif(5), runif(5))
    colnames(df) <- c("trans.padj", "prot.padj", "meta.padj")

    method <- "Edgington"
    expect_is(combinePvalues(df, method = method), "numeric")
    expect_equal(length(combinePvalues(df, method = method)), nrow(df))
})


test_that("Stouffer method works", {
    df <- cbind(runif(5), runif(5), runif(5))
    colnames(df) <- c("trans.padj", "prot.padj", "meta.padj")

    method <- "Stouffer"
    expect_is(combinePvalues(df, method = method), "numeric")
    expect_equal(length(combinePvalues(df, method = method)), nrow(df))
})


test_that("weighted Stouffer method works", {
    df <- cbind(runif(5), runif(5), runif(5))
    colnames(df) <- c("trans.padj", "prot.padj", "meta.padj")

    weights <- sample(1:100, 3, replace = TRUE)
    expect_is(combinePvalues(df, weights = weights), "numeric")
    expect_equal(length(combinePvalues(df)), nrow(df))

    weights <- sample(1:100, 2, replace = TRUE)
    expect_error(combinePvalues(df, weights = weights))

    weights <- c(1, 1, 1)
    expect_equal(combinePvalues(df), combinePvalues(df, weights = weights))
})

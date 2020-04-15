context( "Pathways and feature extraction")

test_that( "transcriptomic features get mapped", {
    
    layer <- c( "transcriptome")
    dbs <- c( "biocarta")
    
    df <- getMultiOmicsFeatures( dbs = dbs, layer = layer)
    
    expect_is( df, "list")
    expect_equal( length( df), 1)
    expect_identical( names( df), c("transcriptome"))
    
    expect_match( names(df$transcriptome), "^\\(BIOCARTA\\)")
    
})

test_that( "proteomic features get mapped", {
    
    layer <- c( "proteome")
    dbs <- c( "reactome")
    
    df <- getMultiOmicsFeatures( dbs = dbs, layer = layer)
    
    expect_is( df, "list")
    expect_equal( length( df), 1)
    expect_identical( names( df), c("proteome"))
    
    expect_match( names(df$proteome), "^\\(REACTOME\\)")
    
})

test_that( "metabolomic features get mapped", {
    
    layer <- c( "metabolome")
    dbs <- c( "kegg")
    
    df <- getMultiOmicsFeatures( dbs = dbs, layer = layer)
    
    expect_is( df, "list")
    expect_equal( length( df), 1)
    expect_identical( names( df), c("metabolome"))
    
    expect_match( names( df$metabolome), "^\\(KEGG\\)")
    expect_match( df$metabolome[[1]], "^HMDB")
    
})


test_that( "each layer has equal length", {
    
    dbs <- c( "kegg", "reactome")
    df <- getMultiOmicsFeatures( dbs = dbs, organism = "mmusculus")
    
    expect_equal( length( df), 3)
    expect_equal( length( df$transcriptome), length( df$proteome))
    expect_equal( length( df$transcriptome), length( df$metabolome))
    
})


test_that( "all organisms can be used", {
    
    dbs <- c( "kegg", "reactome")
    organisms <- getOrganisms()
    organisms <- organisms[ organisms != "xlaevis"]
    lapply( organisms, function( org) {
        df <- getMultiOmicsFeatures( dbs = dbs, organism = org)
        expect_equal( length( df), 3)
        expect_equal( length( df$transcriptome), length( df$proteome))
        expect_equal( length( df$transcriptome), length( df$metabolome))
    })
    
    df <- getMultiOmicsFeatures( dbs = c( "kegg"), organism = "xlaevis")
    expect_equal( length( df), 3)
    expect_equal( length( df$transcriptome), length( df$proteome))
    expect_equal( length( df$transcriptome), length( df$metabolome))
    
})


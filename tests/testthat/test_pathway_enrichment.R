context("Pathway enrichment")

# load transcriptomic features
data( transcriptome)

# load proteomic features
data( proteome)

# load metabolomic features
data( metabolome)


test_that( "provided data can be ranked", {
    
    rF <- rankFeatures( transcriptome$logFC, transcriptome$pValue)
    expect_is( rF, "numeric")
    expect_equal( length( rF), nrow( transcriptome))
    
    rF <- rankFeatures( proteome$logFC, proteome$pValue)
    expect_is( rF, "numeric")
    expect_equal( length( rF), nrow( proteome))
    
    rF <- rankFeatures( metabolome$logFC, metabolome$pValue)
    expect_is( rF, "numeric")
    expect_equal( length( rF), nrow( metabolome))
    

})

test_that( "pathway enrichment works.", {
    
    dbs <- "kegg"
    df <- getMultiOmicsFeatures( dbs = dbs)

    rF <- rankFeatures( proteome$logFC, proteome$pValue)
    expect_is( rF, "numeric")
    expect_equal( length( rF), nrow( proteome))
        
})
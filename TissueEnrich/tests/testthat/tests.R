
data<-system.file("extdata", "combined-proteincodingGenedataCombine.txt", package = "TissueEnrich")
expressionData<-read.table(data,header=TRUE,row.names=1,sep='\t')
TSGenes<-teGenesRetrieval(expressionData)

test_that("checking TS Genes Matrix from Gene reterival", {
  expect_equal(ncol(TSGenes), 3)
})

test_that("checking number of genes in TS Genes Matrix from Gene reterival", {
  expect_equal(nrow(expressionData), length(unique(TSGenes$Gene)))
})

genes<-system.file("extdata", "inputGenes.txt", package = "TissueEnrich")
inputGenes<-scan(genes,character())
output<-teGeneEnrichment(inputGenes = inputGenes,geneFormat=2)

test_that("check the length of the output from teGeneEnrichment", {
  expect_equal(length(output), 3)
})

test_that("check number of columns of the enrichment matrix", {
  expect_equal(ncol(output[[1]]), 3)
})

test_that("check number of tissues in the enrichment matrix using Human protein atlas output", {
  expect_equal(nrow(output[[1]]), 35)
})



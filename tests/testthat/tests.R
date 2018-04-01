
data<-system.file("extdata", "test.expressiondata.txt", package = "TissueEnrich")
expressionData<-read.table(data,header=TRUE,row.names=1,sep='\t')
se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData = row.names(expressionData),colData = colnames(expressionData))
output<-teGeneRetrieval(se)
TSGenes<-data.frame(assay(output))
colnames(TSGenes)<-colData(output)[,1]

test_that("checking TS Genes Matrix from Gene retrieval", {
  expect_equal(ncol(TSGenes), 3)
})

test_that("checking number of genes in TS Genes Matrix from Gene retrieval", {
  expect_equal(nrow(expressionData), length(unique(TSGenes$Gene)))
})

genes<-system.file("extdata", "inputGenes.txt", package = "TissueEnrich")
inputGenes<-scan(genes,character())
gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
output<-teEnrichment(inputGenes = gs)
TSGEOutput<-data.frame(assay(output[[1]]))

test_that("check the length of the output from teGeneEnrichment", {
  expect_equal(length(output), 3)
})

test_that("check number of columns of the enrichment matrix", {
  expect_equal(ncol(TSGEOutput), 2)
})

test_that("check number of tissues in the enrichment matrix using Human protein atlas output", {
  expect_equal(nrow(TSGEOutput), 35)
})



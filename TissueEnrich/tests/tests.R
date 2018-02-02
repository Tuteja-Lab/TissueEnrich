test.examples <- function()
{
  data<-system.file("extdata", "combined-proteincodingGenedataCombine.txt", package = "TissueEnrich")
  expressionData<-read.table(data,header=TRUE,row.names=1,sep='\t')
  TSGenes<-tissueSpecificGenesRetrieval(expressionData)
  checkEqualsNumeric(3, ncol(TSGenes))
  checkEqualsNumeric(ncol(expressionData), length(unique(TSGenes$Gene)))
}

test.deactivation <- function()
{
  DEACTIVATED('Deactivating this test function')
}

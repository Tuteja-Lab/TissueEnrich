##To Supress Note
utils::globalVariables(c("dataset", "%>%","Gene","Gene.name","Gene.stable.ID"
                         ,".","Group", "geneIds","geneIdType","SimpleList"))

#' Define tissue-specific genes by using the algorithm from the Human Protein
#' Atlas
#' @description The teGeneRetrieval function is used to define tissue-specific
#' genes, using the algorithm
#' from the HPA (Uhlén et al. 2015). It takes a gene expression
#' SummarizedExperiment object as input
#' (rows as genes and columns as tissue) and classifies the genes into
#' different gene groups. The users also have the option of changing the
#' default thresholds to vary the degree of tissue specificity of genes. More
#' details about the gene groups and HPA thresholds are provided below. More
#' details about the gene groups are provided in the vignette.
#' @author Ashish Jain, Geetu Tuteja
#' @param expressionData A SummarizedExperiment object containing gene
#' expression values.
#' @param foldChangeThreshold A numeric Threshold of fold change, default 5.
#' @param maxNumberOfTissues A numeric Maximum number of tissues in a group for
#' group enriched genes, default 7.
#' @param expressedGeneThreshold A numeric Minimum gene expression cutoff for
#' the gene to be called as expressed, default 1.
#' @export
#' @return The output is a SummarizedExperiment object containing the
#' information about the tissue-specific genes with three columns:
#' Gene, Tissue, and Enrichment group of the gene in the given tissue.
#' @examples
#' library(SummarizedExperiment)
#' data<-system.file("extdata", "test.expressiondata.txt", package =
#' "TissueEnrich")
#' expressionData<-read.table(data,header=TRUE,row.names=1,sep='\t')
#' se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),
#' rowData = row.names(expressionData),colData = colnames(expressionData))
#' output<-teGeneRetrieval(se)
#' head(assay(output))


teGeneRetrieval <- function(expressionData,foldChangeThreshold=5,
                          maxNumberOfTissues=7,expressedGeneThreshold=1)
{
  ###Add checks for the conditions
  expressionData <- ensurer::ensure_that(expressionData, !is.null(.) &&
                                         class(.) == "SummarizedExperiment"
                                       && !is.null(assay(.))
                                       && (nrow(assay(.)) > 0)
                                       && (ncol(assay(.)) > 1)
                                       && (ncol(rowData(.)) == 1)
                                       && (ncol(colData(.)) == 1),
                                       err_desc = "expressionData should be a
                                       non-empty SummarizedExperiment object
                                       with atleast 1 gene and 2 tissues.")
  foldChangeThreshold <- ensurer::ensure_that(foldChangeThreshold, !is.null(.)
                                            && is.numeric(.) && (. >=1),
                                            err_desc =
                                            "foldChangeThreshold should be a
                                            numeric value greater than
                                            or equal to 1.")
  maxNumberOfTissues <- ensurer::ensure_that(maxNumberOfTissues, !is.null(.) &&
                                             is.numeric(.) && (. >=2),
                                             err_desc =
                                             "maxNumberOfTissues should be an
                                              integer value greater than or
                                              equal to 2.")

  expressedGeneThreshold <- ensurer::ensure_that(expressedGeneThreshold,
                                      !is.null(.) && is.numeric(.)
                                      && (. >=0),err_desc =
                                      "expressedGeneThreshold should be a
                                        numeric value greater than or
                                        equal to 0.")
  SEExpressionData <- expressionData
  minNumberOfTissues <- 2

  ###Creating the expression dataframe object from summarized object.
  expressionData <- setNames(data.frame(assay(expressionData),row.names =
                                        rowData(expressionData)[,1]),
                           colData(expressionData)[,1])
  gene<-1
  start.time <- Sys.time()
  geneList<-as.list(row.names(expressionData))
  x <- apply(expressionData,1, FUN =
                function(j){
                  df<-c()
                  tpm <- j
                  tpm <- sort(tpm,decreasing = TRUE)
                  highTPM <- tpm[1]

                  ####Check for Not Expressed
                  if(highTPM >= expressedGeneThreshold)
                  {
                    secondHighTPM <- tpm[2]
                    foldChangeHigh <- highTPM/secondHighTPM
                    ###Check for Tissue Enriched
                    if(foldChangeHigh >= foldChangeThreshold)
                    {
                      df <- cbind(Gene = gene,Tissue = names(tpm)[1],
                                  Group = "Tissue-Enriched")
                    }else
                    {
                      ####Check for Group Enriched
                      thresholdForGroupTPM <- highTPM/foldChangeThreshold
                      groupTPM <- tpm[(tpm >= thresholdForGroupTPM) &
                                        (tpm >= expressedGeneThreshold),drop=FALSE]
                      isFound <- FALSE
                      if(length(groupTPM) <= maxNumberOfTissues &&
                         length(groupTPM) >= minNumberOfTissues)
                      {
                        # df<-lapply(2:(length(groupTPM)),FUN = function(i){
                        #
                        #   meanTPMForGroup <- mean(groupTPM[seq(1,i)])
                        #   highestTPMOutsideGroup <- tpm[i+1]
                        #   if((meanTPMForGroup/highestTPMOutsideGroup) >= foldChangeThreshold)
                        #   {
                        #     ##To put this genes as group enriched for only the
                        #     ##tissues which satisfy the threhold or all the
                        #     ##tissues??
                        #     x <- lapply(seq(1,i), FUN =
                        #                   function(j){
                        #                     cbind(Gene=gene,
                        #                           Tissue=names(tpm)[j],
                        #                           Group="Group-Enriched")});
                        #     df <- do.call("rbind", x)
                        #     isFound <<- TRUE
                        #     return(df)
                        #     break;
                        #   }
                        #
                        # })

                        for(i in 2:(length(groupTPM)))
                        {
                          meanTPMForGroup <- mean(groupTPM[seq(1,i)])
                          highestTPMOutsideGroup <- tpm[i+1]
                          if((meanTPMForGroup/highestTPMOutsideGroup) >=
                             foldChangeThreshold)
                          {
                            x <- lapply(seq(1,i), FUN =
                                          function(j){
                                            cbind(Gene=gene,
                                                  Tissue=names(tpm)[j],
                                                  Group="Group-Enriched")});
                            df <- do.call("rbind", x)
                            isFound <- TRUE
                            break;
                          }
                        }

                      }
                      if(!isFound)
                      {
                        ####Check for Expressed In All
                        if(all(tpm >= expressedGeneThreshold))
                        {
                          #isFound <- TRUE
                          df <- cbind(Gene=gene,Tissue="All",Group=
                                        "Expressed-In-All")
                        }else
                        {
                          ####Check for Tissue Enhanced
                          tissueEnhancedThreshold <- mean(tpm) * foldChangeThreshold
                          enhancedGene <- tpm[(tpm >= tissueEnhancedThreshold)
                                              & (tpm >= expressedGeneThreshold),
                                              drop=FALSE]
                          if(length(enhancedGene) >=1 )
                          {
                            #isFound <- TRUE
                            x <- lapply(names(enhancedGene), FUN =
                                          function(enhancedTissue){
                                            cbind(Gene=gene,
                                                  Tissue=enhancedTissue,
                                                  Group="Tissue-Enhanced")
                                          });
                            df <- do.call("rbind", x)
                          }else
                          {
                            df <- cbind(Gene=gene,Tissue="All",Group="Mixed")
                          }
                        }
                      }
                    }
                  }else
                  {
                    df <- cbind(Gene=gene,Tissue="All",Group="Not-Expressed")
                  }
                  gene<<-gene + 1
                  return(df)
                });

  TSGenes <- do.call("rbind", x)
  TSGenes[,"Gene"]<-unlist(geneList[as.numeric(TSGenes[,"Gene"])])
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  return(SummarizedExperiment(assays = SimpleList(TSGenes),
                              colData = colnames(TSGenes)))
}

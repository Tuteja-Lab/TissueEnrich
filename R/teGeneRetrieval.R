##To Supress Note
utils::globalVariables(c("dataset", "%>%","Gene","Gene.name","Gene.stable.ID"
                         ,".","Human.gene.name","Human.gene.stable.ID","Group",
                         "Tissue","metadata","geneIds","geneIdType",
                         "ENSEMBLIdentifier","Log10PValue","SimpleList"))

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
                                             is.numeric(.) && (. >=2),err_desc =
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
  notExpressed <- data.frame()
  expressedInAll <- data.frame()
  mixedGenes <- data.frame()
  tissueEnriched <- data.frame()
  groupEnriched <- data.frame()
  tissueEnhanced <- data.frame()

  ###Creating the expression dataframe object from summarized object.
  expressionData <- setNames(data.frame(assay(expressionData),row.names =
                                        rowData(expressionData)[,1]),
                           colData(expressionData)[,1])
  for(gene in row.names(expressionData))
  {
    tpm <- data.frame(t(expressionData[gene, ,drop = FALSE]))
    ##Sorting the genes based on expression
    ##Needed to group the genes based on different categories
    #print(tpm)
    tpm <- tpm[order(-tpm[gene]), ,drop = FALSE]
    highTPM <- tpm[1,]

    ####Check for Not Expressed
    if(highTPM >= expressedGeneThreshold)
    {
      secondHighTPM <- tpm[2,]
      foldChangeHigh <- highTPM/secondHighTPM
      ###Check for Tissue Enriched
      if(foldChangeHigh >= foldChangeThreshold)
      {
        df <- cbind(Gene = gene,Tissue = row.names(tpm)[1],
                  Group = "Tissue-Enriched")
        tissueEnriched <- rbind(tissueEnriched,df)
      }else
      {
        ####Check for Group Enriched
        thresholdForGroupTPM <- highTPM/foldChangeThreshold
        groupTPM <- tpm[(tpm[gene] >= thresholdForGroupTPM) &
                        (tpm[gene] >= expressedGeneThreshold), ,drop=FALSE]
        isFound <- FALSE
        if(nrow(groupTPM) <= maxNumberOfTissues &&
           nrow(groupTPM) >= minNumberOfTissues)
        {
          for(i in 2:(nrow(groupTPM)))
          {
            meanTPMForGroup <- mean(groupTPM[seq(1,i),])
            highestTPMOutsideGroup <- tpm[i+1,]
            if((meanTPMForGroup/highestTPMOutsideGroup) >= foldChangeThreshold)
            {
              ##To put this genes as group enriched for only the
              ##tissues which satisfy the threhold or all the
              ##tissues??
              x <- lapply(seq(1,i), FUN =
                            function(j){
                              cbind(Gene=gene,
                                    Tissue=row.names(tpm)[j],
                                    Group="Group-Enriched")});
              df <- do.call("rbind", x)
              groupEnriched <- rbind(groupEnriched,df)

              # for(j in seq(1,i))
              # {
              #   df<-cbind(Gene=gene,Tissue=row.names(tpm)[j],
              #             Group="Group-Enriched")
              #   groupEnriched<-rbind(groupEnriched,df)
              # }
              isFound <- TRUE
              break;
            }
          }
          if(!isFound)
          {
            ####Check for Expressed In All
            if(tpm[nrow(tpm),] >= expressedGeneThreshold)
            {
              isFound <- TRUE
              df <- cbind(Gene=gene,Tissue="All",Group=
                          "Expressed-In-All")
              expressedInAll <- rbind(expressedInAll,df)
            }else
            {
              ####Check for Tissue Enhanced
              tissueEnhancedThreshold <- mean(tpm[,gene]) * foldChangeThreshold
              enhancedGene <- tpm[(tpm[gene] >= tissueEnhancedThreshold) &
                                (tpm[gene] >=expressedGeneThreshold),
                                ,drop=FALSE]
              if(nrow(enhancedGene) >=1 )
              {
                isFound <- TRUE
                x <- lapply(row.names(enhancedGene), FUN =
                              function(enhancedTissue){
                                cbind(Gene=gene,Tissue=enhancedTissue,
                                      Group="Tissue-Enhanced")
                                });
                df <- do.call("rbind", x)
                tissueEnhanced <- rbind(tissueEnhanced,df)
                # for(enhancedTissue in row.names(enhancedGene))
                # {
                #   df<-cbind(Gene=gene,Tissue=enhancedTissue,
                #             Group="Tissue-Enhanced")
                #   tissueEnhanced<-rbind(tissueEnhanced,df)
                # }
              }
              if(!isFound)
              {
                df <- cbind(Gene=gene,Tissue="All",Group="Mixed")
                mixedGenes <- rbind(mixedGenes,df)
                #mixedGenes<-c(mixedGenes,gene)
              }
            }
          }
        }else
        {
          ###Put the same code here also
          ####Check for Expressed In All
          if(tpm[nrow(tpm),] >= expressedGeneThreshold)
          {
            isFound <- TRUE
            df <- cbind(Gene=gene,Tissue="All",Group="Expressed-In-All")
            expressedInAll <- rbind(expressedInAll,df)
          }else
          {
            ####Check for Tissue Enhanced
            tissueEnhancedThreshold <- mean(tpm[,gene]) * foldChangeThreshold
            enhancedGene <- tpm[(tpm[gene] > tissueEnhancedThreshold)
                          & (tpm[gene] >= expressedGeneThreshold), ,drop=FALSE]
            if(nrow(enhancedGene) >= 1)
            {
              isFound <- TRUE
              x <- lapply(row.names(enhancedGene), FUN =
                            function(enhancedTissue){
                              cbind(Gene=gene,Tissue=enhancedTissue,
                                    Group="Tissue-Enhanced")
                              });
              df <- do.call("rbind", x)
              tissueEnhanced <- rbind(tissueEnhanced,df)
              # for(enhancedTissue in row.names(enhancedGene))
              # {
              #   df<-cbind(Gene=gene,Tissue=enhancedTissue,
              #             Group="Tissue-Enhanced")
              #   tissueEnhanced<-rbind(tissueEnhanced,df)
              # }
            }
            if(!isFound)
            {
              df <- cbind(Gene=gene,Tissue="All",Group="Mixed")
              mixedGenes <- rbind(mixedGenes,df)
            }
          }
        }
      }
    }else
    {
      df <- cbind(Gene=gene,Tissue="All",Group="Not-Expressed")
      notExpressed <- rbind(notExpressed,df)
    }
  }
  TSGenes <- rbind(tissueEnriched,groupEnriched,tissueEnhanced,
                 notExpressed,expressedInAll,mixedGenes)
  return(SummarizedExperiment(assays = SimpleList(as.matrix(TSGenes)),
                              colData = colnames(TSGenes)))
}

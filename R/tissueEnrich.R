#ensure_numeric<-ensure_that(is.numeric(.),"Please enter numeric values.")
library(dplyr)
library(ensurer)
library(utils)
library(SummarizedExperiment)
library(GSEABase)


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
#' @param expressionData A SummarizedExperiment object containing gene expression values.
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


teGeneRetrieval<-function(expressionData,foldChangeThreshold=5,
                          maxNumberOfTissues=7,expressedGeneThreshold=1)
{
    ###Add checks for the conditions
    expressionData<-ensurer::ensure_that(expressionData, !is.null(.) &&
                                           class(.) == "SummarizedExperiment"
                                         && !is.null(assay(.)) && (nrow(assay(.)) > 0)
                                         && (ncol(assay(.)) > 1) && (ncol(rowData(.)) == 1)
                                         && (ncol(colData(.)) == 1),
                                         err_desc = "expressionData should be a
                                          non-empty SummarizedExperiment object with atleast 1 gene and 2 tissues.")
    foldChangeThreshold<-ensurer::ensure_that(foldChangeThreshold, !is.null(.) && is.numeric(.) && (. >=1),err_desc = "foldChangeThreshold should be a numeric value greater than or equal to 1.")
    maxNumberOfTissues<-ensurer::ensure_that(maxNumberOfTissues, !is.null(.) && is.numeric(.) && (. >=2),err_desc = "maxNumberOfTissues should be an integer value greater than or equal to 2.")
    expressedGeneThreshold<-ensurer::ensure_that(expressedGeneThreshold, !is.null(.) && is.numeric(.) && (. >=0),err_desc = "expressedGeneThreshold should be a numeric value greater than or equal to 0.")

    SEExpressionData<-expressionData
    minNumberOfTissues<-2
    notExpressed<-data.frame()
    expressedInAll<-data.frame()
    mixedGenes<-data.frame()
    tissueEnriched<-data.frame()
    groupEnriched<-data.frame()
    tissueEnhanced<-data.frame()

    ###Creating the expression dataframe object from summarized object.
    expressionData<-setNames(data.frame(assay(expressionData),row.names = rowData(expressionData)[,1]),colData(expressionData)[,1])
    for(gene in row.names(expressionData))
    {
        tpm<-data.frame(t(as.matrix(expressionData[gene,])))
        colnames(tpm)<-c(gene)
        tpm<-tpm[with(tpm, order(-tpm[gene])), ,drop = FALSE]
        highTPM<-tpm[1,]
        ####Check for Not Expressed
        if(highTPM>=expressedGeneThreshold)
        {
          secondHighTPM<-tpm[2,]
          foldChangeHigh<-highTPM/secondHighTPM
          ###Check for Tissue Enriched
          if(foldChangeHigh >=foldChangeThreshold)
          {
            df<-cbind(Gene=gene,Tissue=row.names(tpm)[1],Group="Tissue-Enriched")
            tissueEnriched<-rbind(tissueEnriched,df)
          }else
          {
            ####Check for Group Enriched
            thresholdForGroupTPM<-highTPM/foldChangeThreshold
            groupTPM<-tpm[(tpm[gene] >= thresholdForGroupTPM) & (tpm[gene] >=expressedGeneThreshold), ,drop=FALSE]
            isFound<-FALSE
            if(nrow(groupTPM)<=maxNumberOfTissues && nrow(groupTPM) >= minNumberOfTissues)
            {

              for(i in 2:(nrow(groupTPM)))
              {
                meanTPMForGroup<-mean(groupTPM[seq(1,i),])
                highestTPMOutsideGroup<-tpm[i+1,]
                if((meanTPMForGroup/highestTPMOutsideGroup) >= foldChangeThreshold)
                {
                  ##To put this genes as group enriched for only the tissues which satisfy the threhold or all the tissues??
                  for(j in seq(1,i))
                  {
                    df<-cbind(Gene=gene,Tissue=row.names(tpm)[j],Group="Group-Enriched")
                    groupEnriched<-rbind(groupEnriched,df)
                  }
                  isFound<-TRUE
                  break;
                }
              }
              if(!isFound)
              {
                ####Check for Expressed In All
                if(tpm[nrow(tpm),] >= expressedGeneThreshold)
                {
                  isFound<-TRUE
                  df<-cbind(Gene=gene,Tissue="All",Group="Expressed-In-All")
                  expressedInAll<-rbind(expressedInAll,df)
                }else
                {
                  ####Check for Tissue Enhanced
                  tissueEnhancedThreshold<-mean(tpm[,gene]) * foldChangeThreshold
                  enhancedGene<-tpm[(tpm[gene] >= tissueEnhancedThreshold) & (tpm[gene] >=expressedGeneThreshold), ,drop=FALSE]
                  if(nrow(enhancedGene)>=1)
                  {
                    isFound<-TRUE
                    for(enhancedTissue in row.names(enhancedGene))
                    {
                      df<-cbind(Gene=gene,Tissue=enhancedTissue,Group="Tissue-Enhanced")
                      tissueEnhanced<-rbind(tissueEnhanced,df)
                    }
                  }
                  if(!isFound)
                  {
                    df<-cbind(Gene=gene,Tissue="All",Group="Mixed")
                    mixedGenes<-rbind(mixedGenes,df)
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
                isFound<-TRUE
                df<-cbind(Gene=gene,Tissue="All",Group="Expressed-In-All")
                expressedInAll<-rbind(expressedInAll,df)
              }else
              {
                ####Check for Tissue Enhanced
                tissueEnhancedThreshold<-mean(tpm[,gene]) * foldChangeThreshold
                enhancedGene<-tpm[(tpm[gene] > tissueEnhancedThreshold) & (tpm[gene] >=expressedGeneThreshold), ,drop=FALSE]
                if(nrow(enhancedGene)>=1)
                {
                  isFound<-TRUE
                  for(enhancedTissue in row.names(enhancedGene))
                  {
                    df<-cbind(Gene=gene,Tissue=enhancedTissue,Group="Tissue-Enhanced")
                    tissueEnhanced<-rbind(tissueEnhanced,df)
                  }
                }
                if(!isFound)
                {
                  df<-cbind(Gene=gene,Tissue="All",Group="Mixed")
                  mixedGenes<-rbind(mixedGenes,df)
                }
              }
            }
          }
        }else
        {
          df<-cbind(Gene=gene,Tissue="All",Group="Not-Expressed")
          notExpressed<-rbind(notExpressed,df)
        }
    }
    TSGenes<-rbind(tissueEnriched,groupEnriched,tissueEnhanced,notExpressed,expressedInAll,mixedGenes)
    return(SummarizedExperiment(assays = SimpleList(as.matrix(TSGenes)),colData = colnames(TSGenes)))
}

#' Calculate tissue-specific gene enrichment using the hypergeometric test
#'
#' @description The teEnrichment function is used to calculate the enrichment of tissue-specific genes,
#' given an input gene set. It uses tissue-specific genes defined by processing RNA-Seq datasets
#' from human and mouse.
#' @author Ashish Jain, Geetu Tuteja
#' @param inputGenes A GeneSet object containing the input genes, organism type ("Homo Sapiens" or "Mus Musculus"), and gene id identifier (Gene Symbol or ENSEMBL identifier).
#' @param rnaSeqDataset An integer describing the dataset to be used for enrichment analysis. 1 for "Human Protein Atlas" (default), 2 for "GTEx", 3 for "Mouse ENCODE". Default 1.
#' @param tissueSpecificGeneType An integer describing the type of tissue-specific genes to be used.
#' 1 for "All" (default), 2 for "Tissue-Enriched", 3 for "Tissue-Enhanced", and 4 for "Group-Enriched". Default 1.
#' @param multiHypoCorrection Flag to correct P-values for multiple hypothesis using BH method. Default TRUE.
#' @export
#' @return The output is a list with three objects. The first object is the
#' SummarizedExperiment object containing the enrichment results, the second
#' object contains the expression values and tissue-specificity information
#' of the tissue-specific genes for genes from the input gene set, and the
#' third is a GeneSet object containing genes that were not identified in the
#' tissue-specific gene data.
#'
#' @examples
#' library(dplyr)
#' library(ggplot2)
#' genes<-system.file("extdata", "inputGenes.txt", package = "TissueEnrich")
#' inputGenes<-scan(genes,character())
#' gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",geneIdType=SymbolIdentifier())
#' output<-teEnrichment(gs)
#' seEnrichmentOutput<-output[[1]]
#' enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput)),
#' colData(seEnrichmentOutput)[,1])
#' enrichmentOutput$Tissue<-row.names(enrichmentOutput)
#' #Plotting the P-Values
#' ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,
#' label = Tissue.Specific.Genes,fill = Tissue))+
#' geom_bar(stat = 'identity')+
#' labs(x='', y = '-LOG10(P-Value)')+
#' theme_bw()+
#' theme(legend.position="none")+
#' theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
#' theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#' panel.grid.major= element_blank(),panel.grid.minor = element_blank())


teEnrichment<-function(inputGenes = NULL,
                                       rnaSeqDataset=1,#c("Human Protein Atlas","GTEx combine","GTEx sub-tissue","Mouse ENCODE"),
                                       tissueSpecificGeneType=1,multiHypoCorrection=TRUE)
{
  ###Add checks for the conditions
  #inputGenes<-ensurer::ensure_that(inputGenes, !is.null(.) && is.vector(.),err_desc = "Please enter correct inputGenes. It should be a character vector")
  inputGenes<-ensurer::ensure_that(inputGenes, !is.null(.) && (class(.) == "GeneSet") && !is.null(geneIds(.)),err_desc = "Please enter correct inputGenes. It should be a Gene set object.")
  rnaSeqDataset<-ensurer::ensure_that(rnaSeqDataset, !is.null(.) && is.numeric(.) && . <=3 || . >=1,err_desc = "Please enter correct rnaSeqDataset. It should be 1 for Human Protein Atlas, 2 for GTEx combine, 3 for GTEx sub-tissue, 4 for Mouse ENCODE.")
  #organism<-ensurer::ensure_that(organism, !is.null(.) && is.numeric(.) && . <=2 || . >=1,err_desc = "Please enter correct organism. It should be either 1 for Homo Sapiens, 2 for Mus Musculus.")
  #geneFormat<-ensurer::ensure_that(geneFormat,!is.null(.) && is.numeric(.) && . <=2 || . >=1,err_desc = "Please enter correct geneFormat. It should be either 1 for EnsemblId, 2 for Gene Symbol.")
  tissueSpecificGeneType<-ensurer::ensure_that(tissueSpecificGeneType,!is.null(.) && is.numeric(.)  && . <=4 || . >=1,err_desc = "Please enter correct tissueSpecificGeneType. It should be 1 for All, 2 for Tissue-Enriched,3 for Tissue-Enhanced, and 4 for Group-Enriched")
  #isHomolog<-ensurer::ensure_that(isHomolog,!is.null(.) && is.logical(.) ,err_desc = "Please enter correct isHomolog. It should be either TRUE or FALSE")
  multiHypoCorrection<-ensurer::ensure_that(multiHypoCorrection,!is.null(.) && is.logical(.) ,err_desc = "Please enter correct multiHypoCorrection. It should be either TRUE or FALSE")
  ##Enviroment to load datasets
  #env<-loading("data/combine-expression.rda")
  # withCallingHandlers(warning("Column `Gene` joining factors with different levels, coercing to character vector"), warning = function(w) {
  #   print("No")
  # })
  geInputGenes<-inputGenes
  inputGenes<-geneIds(inputGenes)
  ##Load gene mapping and orthologs Data
  isHomolog<-FALSE
  e <- new.env()
  load(file = system.file("extdata", "combine-expression.rda", package = "TissueEnrich"),envir = e)
  if(rnaSeqDataset == 1)
  {
    expressionDataLocal<-e$dataset$`Protein-Atlas`$expressionData
    tissueDetails<-e$dataset$`Protein-Atlas`$tissueDetails
    tissueSpecificGenes<-e$dataset$`Protein-Atlas`$tissueSpecificGenes
  }else if(rnaSeqDataset == 2)
  {
    expressionDataLocal<-e$dataset$`GTEx-Combine`$expressionData
    tissueDetails<-e$dataset$`GTEx-Combine`$tissueDetails
    tissueSpecificGenes<-e$dataset$`GTEx-Combine`$tissueSpecificGenes
  }else if(rnaSeqDataset == 3)
  {
    expressionDataLocal<-e$dataset$`ENCODE Dataset`$expressionData
    tissueDetails<-e$dataset$`ENCODE Dataset`$tissueDetails
    tissueSpecificGenes<-e$dataset$`ENCODE Dataset`$tissueSpecificGenes
  }else
  {
    stop("Please enter correct dataset id.",call. = FALSE)
    #print("Please enter correct dataset")
  }


  ##Check for organism and homolog to update geneMapping Variable
  if(organism(geInputGenes) == "Homo Sapiens")#"Homo Sapiens")
  {
    organism <-1
    geneMapping<-e$dataset$humanGeneMapping
    if(rnaSeqDataset == 3)
    {
      isHomolog<-TRUE
    }
  }else if(organism(geInputGenes) == "Mus Musculus")#"Mus Musculus")
  {
    organism <-2
    geneMapping<-e$dataset$mouseGeneMapping
    if(rnaSeqDataset == 1 || rnaSeqDataset == 2)
    {
      isHomolog<-TRUE
    }
  }else
  {
    stop("Please enter correct organism.",call. = FALSE)
    #print("Please enter correct dataset")
  }

  if(isHomolog)
  {
    geneMapping<-e$dataset$mouseHumanOrthologs
  }

  ###Update tissueSpecificGene variable based on the group
  if(tissueSpecificGeneType == 1)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes
  }else if(tissueSpecificGeneType == 2)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes %>% dplyr::filter(Group == "Tissue-Enriched")
  }else if(tissueSpecificGeneType == 3)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes %>% dplyr::filter(Group == "Tissue-Enhanced")
  }else if(tissueSpecificGeneType == 4)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes %>% dplyr::filter(Group == "Group-Enriched")
  }else
  {
    stop("Tissue-specific gene type is not correct.",call. = FALSE)
    #print("Something is wrong!")
  }

  if(geneIdType(geneIdType(geInputGenes)) == geneIdType(SymbolIdentifier()))
  {
    geneFormat<-2
  }else if(geneIdType(geneIdType(geInputGenes)) == geneIdType(ENSEMBLIdentifier()))
  {
    geneFormat<-1
  }else{
    stop(paste0("Gene Id type is not correct."),call. = FALSE)
  }

  inputGenes<-toupper(inputGenes)
  inputGenes<-unique(inputGenes)
  ###Check if it is ortholog comparison or not#######
  if(isHomolog)
  {
    ####Calculation when working with orthologs
    #####Convert the Gene Id########
    genesNotFound<-c()
    if(organism == 1)
    {
      if(geneFormat == 2)
      {
        inputEnsemblGenes<-geneMapping %>% dplyr::filter(Human.gene.name %in% inputGenes) %>% dplyr::select(Gene.stable.ID)
      }else
      {
        inputEnsemblGenes<-geneMapping %>% dplyr::filter(Human.gene.stable.ID %in% inputGenes) %>% dplyr::select(Gene.stable.ID)
      }
      inputEnsemblGenes<-as.character(inputEnsemblGenes$Gene.stable.ID)

      ##Updating the geneMapping
      geneMappingForCurrentDataset<-geneMapping %>% dplyr::filter(Gene.stable.ID %in% row.names(expressionDataLocal))

      #####Check for genes which are not there in our list########
      if(length(inputGenes) != length(intersect(as.character(geneMappingForCurrentDataset$Gene.stable.ID),inputEnsemblGenes)))
      {
        if(geneFormat == 2)
        {
          genesNotFoundList<-base::setdiff(inputGenes,as.character(geneMappingForCurrentDataset$Human.gene.name))
        }else
        {
          genesNotFoundList<-base::setdiff(inputGenes,as.character(geneMappingForCurrentDataset$Human.gene.stable.ID))
        }
        genesNotFound<-c(length(genesNotFoundList),genesNotFoundList)
        inputEnsemblGenes<-intersect(row.names(expressionDataLocal),inputEnsemblGenes)
      }

      ###retreving mapping based on gene format.
      if(geneFormat == 2)
      {
        geneMappingForCurrentDataset<-geneMappingForCurrentDataset %>% dplyr::select(Gene.stable.ID,Human.gene.name)
      }else if(geneFormat == 1)
      {
        geneMappingForCurrentDataset<-geneMappingForCurrentDataset %>% dplyr::select(Gene.stable.ID,Human.gene.stable.ID)
      }
      colnames(geneMappingForCurrentDataset)<-c("Gene","Gene.name")

    }else
    {
      if(geneFormat == 2)
      {
        inputEnsemblGenes<-geneMapping %>% dplyr::filter(Gene.name %in% inputGenes) %>% dplyr::select(Human.gene.stable.ID)
      }else
      {
        inputEnsemblGenes<-geneMapping %>% dplyr::filter(Gene.stable.ID %in% inputGenes) %>% dplyr::select(Human.gene.stable.ID)
      }
      inputEnsemblGenes<-as.character(inputEnsemblGenes$Human.gene.stable.ID)
      ##Updating the geneMapping
      geneMappingForCurrentDataset<-geneMapping %>% dplyr::filter(Human.gene.stable.ID %in% row.names(expressionDataLocal))
      #####Check for genes which are not there in our list########
      if(length(inputGenes) != length(intersect(as.character(geneMappingForCurrentDataset$Human.gene.stable.ID),inputEnsemblGenes)))
      {
        if(geneFormat == 2)
        {
          genesNotFoundList<-base::setdiff(inputGenes,as.character(geneMappingForCurrentDataset$Gene.name))
        }else
        {
          genesNotFoundList<-base::setdiff(inputGenes,as.character(geneMappingForCurrentDataset$Gene.stable.ID))
        }
        genesNotFound<-c(length(genesNotFoundList),genesNotFoundList)
        inputEnsemblGenes<-intersect(row.names(expressionDataLocal),inputEnsemblGenes)
      }
      ###retreving mapping based on gene format.
      if(geneFormat == 2)
      {
        geneMappingForCurrentDataset<-geneMappingForCurrentDataset %>% dplyr::select(Human.gene.stable.ID,Gene.name)
      }else if(geneFormat == 1)
      {
        geneMappingForCurrentDataset<-geneMappingForCurrentDataset %>% dplyr::select(Human.gene.stable.ID,Gene.stable.ID)
      }
      colnames(geneMappingForCurrentDataset)<-c("Gene","Gene.name")
    }

    finalTissueSpecificGenes<-finalTissueSpecificGenes %>% dplyr::filter(Gene %in% geneMappingForCurrentDataset$Gene)
  }else{

    ####Normal Calculation
    #####Convert the Gene Id########
    if(geneFormat == 2)
    {
      inputEnsemblGenes<-geneMapping %>% dplyr::filter(Gene.name %in% inputGenes) %>% dplyr::select(Gene)
      inputEnsemblGenes<-as.character(inputEnsemblGenes$Gene)
    }else
    {
      inputEnsemblGenes<-inputGenes
    }
    ##Updating the geneMapping
    geneMappingForCurrentDataset<-geneMapping %>% dplyr::filter(Gene %in% row.names(expressionDataLocal))
    #####Check for genes which are not there in our list########
    genesNotFound<-c()
    if(length(inputGenes) != length(intersect(as.character(geneMappingForCurrentDataset$Gene),inputEnsemblGenes)))
    {
      if(geneFormat == 2)
      {
        genesNotFound<-base::setdiff(inputGenes,as.character(geneMappingForCurrentDataset$Gene.name))
      }else
      {
        genesNotFound<-base::setdiff(inputGenes,as.character(geneMappingForCurrentDataset$Gene))
      }
      inputEnsemblGenes<-intersect(row.names(expressionDataLocal),inputEnsemblGenes)
    }
  }

  ####Calculate the Hypergeometric P-Value#########
  tissueNames<-as.character(tissueDetails$RName)
  pValueList<-c()
  overlapGenesList<-c()
  overlapTissueGenesList<-list()
  for(tissue in tissueNames)
  {
    tissueGenes<-finalTissueSpecificGenes %>% dplyr::filter(Tissue==tissue)
    overlapGenes<-length(intersect(tissueGenes$Gene,inputEnsemblGenes))
    if(isHomolog || geneFormat == 2)
    {
      ##code to convert ensembl Id to gene names
      intGenes<- geneMappingForCurrentDataset %>% dplyr::filter(Gene %in% intersect(tissueGenes$Gene,inputEnsemblGenes))
      teExpressionData<-expressionDataLocal[as.character(intGenes$Gene),]
      #genes<-as.factor(row.names(teExpressionData))
      genes<-row.names(teExpressionData)
      teExpressionData$Gene<-genes

      ##Supress warnings due to different levels in factors
      #oldw <- getOption("warn")
      #options(warn = -1)
      teExpressionData<-left_join(teExpressionData,geneMappingForCurrentDataset, by = "Gene")
      #options(warn = oldw)
      if(isHomolog)
      {
        teExpressionData$Gene.name<-as.character(teExpressionData$Gene.name)
      }
      #row.names(teExpressionData)<-teExpressionData[,ncol(teExpressionData)-1]

      if(nrow(teExpressionData) > 0)
      {
        #####Code to take the mean of the genes with multiple ensembl Ids.
        res <- as.data.frame( # sapply returns a list here, so we convert it to a data.frame
          t(sapply(unique(teExpressionData$Gene.name), # for each unique column name
                   FUN=function(col) colMeans(teExpressionData[teExpressionData$Gene.name == col,c(seq(1,(ncol(teExpressionData)-2)))]) # calculate row means
          )
          )
        )
        teExpressionData<-res
      }else
      {
        teExpressionData<-teExpressionData[,seq(1,(ncol(teExpressionData)-2))]
      }

      #teExpressionData<-teExpressionData[,1:(ncol(teExpressionData)-2)]
      ##code to convert ensembl Id for groups
      teInputGeneGroups<-tissueGenes %>% dplyr::filter(Gene %in% intGenes$Gene) %>% select(Gene,Group)

      ##Supress warnings due to different levels in factors
      #oldw <- getOption("warn")
      #options(warn = -1)
      teInputGeneGroups<-left_join(teInputGeneGroups,geneMappingForCurrentDataset, by = "Gene") %>% select(Gene.name,Group)
      #options(warn = oldw)
      ##Code to remove gene names with
      teInputGeneGroups<-unique(teInputGeneGroups)
    }else
    {
      teExpressionData<- expressionDataLocal[intersect(tissueGenes$Gene,inputEnsemblGenes),]
      teInputGeneGroups<-tissueGenes %>% dplyr::filter(Gene %in% inputEnsemblGenes) %>% select(Gene,Group)
    }
    colnames(teInputGeneGroups) <-c("Gene","Group")
    colnames(teExpressionData)<-tissueDetails$TissueName
    teExpressionData<-log2(teExpressionData+1)
    overlapTissueGenesList[[as.character(tissueDetails[tissueDetails$RName == tissue,"TissueName"])]]<-list(SummarizedExperiment(assays = SimpleList(as.matrix(teExpressionData)),rowData=row.names(teExpressionData),colData = colnames(teExpressionData)),SummarizedExperiment(assays = SimpleList(as.matrix(teInputGeneGroups)),colData = colnames(teInputGeneGroups)))
    GenesInTissue<-nrow(tissueGenes)
    pValue<-stats::phyper(overlapGenes-1,GenesInTissue,nrow(geneMappingForCurrentDataset)-GenesInTissue,length(inputEnsemblGenes),lower.tail = FALSE)
    pValueList<-c(pValueList,pValue)
    overlapGenesList<-c(overlapGenesList,overlapGenes)
  }

  if(multiHypoCorrection)
  {
    ##Multiple hypothesis correction
    pValueList<-stats::p.adjust(pValueList,method = "BH")
  }
  pValueList<-(-log10(pValueList))
  output<-data.frame(Log10PValue=pValueList,Tissue.Specific.Genes=overlapGenesList,row.names = tissueDetails$TissueName)
  output<-output[with(output, order(-Log10PValue)), ]
  #colnames(output)<-c("Tissue","-Log10PValue","Tissue-Specific-Genes")
  return(list(SummarizedExperiment(assays = SimpleList(as.matrix(output)),rowData=row.names(output),colData = colnames(output)),overlapTissueGenesList,GeneSet(geneIds=genesNotFound)))
}


#' Calculate tissue-specific gene enrichment using the hypergeometric test for custom datasets
#'
#' @description The teEnrichmentCustom function is used to calculate tissue-specific gene
#' enrichment using tissue-specific genes defined using the teGeneRetrieval function.
#' @author Ashish Jain, Geetu Tuteja
#' @param inputGenes An GeneSet object containing the input genes.
#' @param tissueSpecificGenes A SummarizedExperiment object. Output from `teGeneRetrieval` function. Default NULL.
#' @param tissueSpecificGeneType An integer describing the type of tissue-specific genes to be used. 1 for "All" (default), 2 for "Tissue-Enriched",3 for "Tissue-Enhanced", and 4 for "Group-Enriched". Default 1.
#' @param multiHypoCorrection Flag to correct P-values for multiple hypothesis using BH method. Default TRUE.
#' @export
#' @return The output is a list with three objects. The first object is the
#' SummarizedExperiment object containing the enrichment results, the second
#' object contains the tissue-specificity information of the tissue-specific
#' genes for genes from the input gene set, and the third is a GeneSet object
#' containing genes that were not identified in the tissue-specific gene data.
#' @examples
#' library(dplyr)
#' data<-system.file("extdata", "test.expressiondata.txt", package = "TissueEnrich")
#' expressionData<-read.table(data,header=TRUE,row.names=1,sep='\t')
#' se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),rowData =
#' row.names(expressionData),colData = colnames(expressionData))
#' output<-teGeneRetrieval(se)
#' head(metadata(output)[["TissueSpecificGenes"]])
#' genes<-system.file("extdata", "inputGenesEnsembl.txt", package = "TissueEnrich")
#' inputGenes<-scan(genes,character())
#' gs<-GeneSet(geneIds=inputGenes)
#' output2<-teEnrichmentCustom(gs,output)
#' #Plotting the P-Values
#' enrichmentOutput<-setNames(data.frame(assay(output2[[1]])),colData(output2[[1]])[,1])
#' enrichmentOutput$Tissue<-row.names(enrichmentOutput)
#' ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,
#' label = Tissue.Specific.Genes,fill = Tissue))+
#' geom_bar(stat = 'identity')+
#' labs(x='', y = '-LOG10(P-Value)')+
#' theme_bw()+
#' theme(legend.position="none")+
#' theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title = element_text(size=15))+
#' theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#' panel.grid.major= element_blank(),panel.grid.minor = element_blank())

teEnrichmentCustom<-function(inputGenes=NULL,tissueSpecificGenes=NULL,
                                             tissueSpecificGeneType= 1,multiHypoCorrection = TRUE)
{
  ###Add checks for the conditions
  inputGenes<-ensurer::ensure_that(inputGenes, !is.null(.) && (class(.) == "GeneSet") && !is.null(geneIds(.)),err_desc = "Please enter correct inputGenes. It should be a Gene set object.")
  ##Check for dataset class, rows and should not be NULL
  tissueSpecificGenes<-ensurer::ensure_that(tissueSpecificGenes, !is.null(.) && class(.) == "SummarizedExperiment" && !is.null(assay(.)) && (nrow(assay(.)) > 0) && (ncol(assay(.)) > 1) && (ncol(colData(.)) == 1),err_desc = "expressionData should be a non-empty SummarizedExperiment object with atleast 1 gene and 2 tissues.")
  tissueSpecificGenes<-setNames(data.frame(assay(tissueSpecificGenes)),colData(tissueSpecificGenes)[,1])
  tissueSpecificGenes<-ensurer::ensure_that(tissueSpecificGenes,!is.null(.) && is.data.frame(.) && !is.null(.) && (nrow(.) > 0) && (ncol(.) == 3),err_desc = "Please enter correct tissueSpecificGenes. It should be a non-empty dataframe object.")
  tissueSpecificGeneType<-ensurer::ensure_that(tissueSpecificGeneType,!is.null(.) && is.numeric(.)  && . <=4 || . >=1 ,err_desc = "Please enter correct tissueSpecificGeneType. It should be 1 for All, 2 for Tissue-Enriched,3 for Tissue-Enhanced, and 4 for Group-Enriched")
  multiHypoCorrection<-ensurer::ensure_that(multiHypoCorrection,!is.null(.) && is.logical(.),err_desc = "Please enter correct multiHypoCorrection. It should be either TRUE or FALSE")

  if(tissueSpecificGeneType == 1)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes %>% dplyr::filter(Group %in% c("Tissue-Enriched","Tissue-Enhanced","Group-Enriched"))
  }else if(tissueSpecificGeneType == 2)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes %>% dplyr::filter(Group == "Tissue-Enriched")
  }else if(tissueSpecificGeneType == 3)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes %>% dplyr::filter(Group == "Tissue-Enhanced")
  }else if(tissueSpecificGeneType == 4)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes %>% dplyr::filter(Group == "Group-Enriched")
  }else
  {
    stop("Tissue-specific gene type is not correct.",call. = FALSE)
    #print("Something is wrong!")
  }

  geInputGenes<-inputGenes
  inputGenes<-geneIds(inputGenes)
  inputGenes<-unique(inputGenes)
  totalGenes<-as.character(unique(tissueSpecificGenes$Gene))
  genesNotFound<-GeneSet(geneIds=base::setdiff(inputGenes,totalGenes))
  inputEnsemblGenes<-intersect(inputGenes,totalGenes)

  ####Calculate the Hypergeometric P-Value#########
  tissueNames<-as.character(unique(finalTissueSpecificGenes$Tissue))
  pValueList<-c()
  overlapGenesList<-c()
  overlapTissueGenesList<-list()
  for(tissue in tissueNames)
  {
    tissueGenes<-finalTissueSpecificGenes %>% dplyr::filter(Tissue==tissue)
    overlapGenes<-length(intersect(tissueGenes$Gene,inputEnsemblGenes))
    teInputGeneGroups<-tissueGenes %>% dplyr::filter(Gene %in% inputEnsemblGenes) %>% select(Gene,Group)
    seTeInputGeneGroups<-SummarizedExperiment(assays = SimpleList(as.matrix(teInputGeneGroups)),colData = colnames(teInputGeneGroups))
    overlapTissueGenesList[[tissue]]<-seTeInputGeneGroups
    GenesInTissue<-nrow(tissueGenes)
    pValue<-stats::phyper(overlapGenes-1,GenesInTissue,length(totalGenes)-GenesInTissue,length(inputEnsemblGenes),lower.tail = FALSE)
    pValueList<-c(pValueList,pValue)
    overlapGenesList<-c(overlapGenesList,overlapGenes)
  }
  if(multiHypoCorrection)
  {
    ##Multiple hypothesis correction
    pValueList<-stats::p.adjust(pValueList,method = "BH")
  }
  pValueList<-(-log10(pValueList))
  output<-data.frame(Log10PValue=pValueList,Tissue.Specific.Genes=overlapGenesList,row.names = tissueNames)
  output<-output[with(output, order(-Log10PValue)), ]
  seOutput<-SummarizedExperiment(assays = SimpleList(as.matrix(output)),rowData=row.names(output),colData = colnames(output))
  return(list(seOutput,overlapTissueGenesList,genesNotFound))
}

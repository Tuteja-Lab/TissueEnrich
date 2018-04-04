##To Supress Note
utils::globalVariables(c("dataset", "%>%","Gene","Gene.name","Gene.stable.ID"
                         ,".","Human.gene.name","Human.gene.stable.ID","Group",
                         "Tissue","metadata","geneIds","geneIdType",
                         "ENSEMBLIdentifier","Log10PValue","SimpleList"))

#' Calculate tissue-specific gene enrichment using the hypergeometric test
#'
#' @description The teEnrichment function is used to calculate the enrichment
#' of tissue-specific genes, given an input gene set. It uses tissue-specific
#' genes defined by processing RNA-Seq datasets from human and mouse.
#' @author Ashish Jain, Geetu Tuteja
#' @param inputGenes A GeneSet object containing the input genes, organism type
#' ("Homo Sapiens" or "Mus Musculus"), and gene id identifier (Gene Symbol or
#' ENSEMBL identifier).
#' @param rnaSeqDataset An integer describing the dataset to be used for
#' enrichment analysis. 1 for "Human Protein Atlas" (default), 2 for "GTEx",
#' 3 for "Mouse ENCODE". Default 1.
#' @param tissueSpecificGeneType An integer describing the type of tissue-
#' specific genes to be used. 1 for "All" (default), 2 for "Tissue-Enriched",
#'  3 for "Tissue-Enhanced", and 4 for "Group-Enriched". Default 1.
#' @param multiHypoCorrection Flag to correct P-values for multiple hypothesis
#' using BH method. Default TRUE.
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
#' gs<-GeneSet(geneIds=inputGenes,organism="Homo Sapiens",
#' geneIdType=SymbolIdentifier())
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
#' theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
#' element_text(size=15))+
#' theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#' panel.grid.major= element_blank(),panel.grid.minor = element_blank())


teEnrichment<-function(inputGenes = NULL,
                       rnaSeqDataset=1,
                       tissueSpecificGeneType=1,
                       multiHypoCorrection=TRUE)
{
  ###Add checks for the conditions
  inputGenes<-ensurer::ensure_that(inputGenes, !is.null(.) &&
                                     (class(.) == "GeneSet") &&
                                     !is.null(geneIds(.)),
                                   err_desc = "Please enter correct inputGenes.
                                   It should be a Gene set object.")
  rnaSeqDataset<-ensurer::ensure_that(rnaSeqDataset, !is.null(.) &&
                                        is.numeric(.) && . <=3 || . >=1,
                                      err_desc = "Please enter correct
                                      rnaSeqDataset. It should be 1 for Human
                                      Protein Atlas, 2 for GTEx combine, 3 for
                                      GTEx sub-tissue, 4 for Mouse ENCODE.")
  tissueSpecificGeneType<-ensurer::ensure_that(tissueSpecificGeneType,
                                               !is.null(.) && is.numeric(.)
                                               && . <=4 || . >=1,err_desc =
                                                 "Please enter correct
                                               tissueSpecificGeneType. It should be
                                               1 for All, 2 for Tissue-Enriched,3 for
                                               Tissue-Enhanced, and 4 for Group-Enriched")

  multiHypoCorrection<-ensurer::ensure_that(multiHypoCorrection,!is.null(.)
                                            && is.logical(.) ,err_desc =
                                              "Please enter correct multiHypoCorrection.
                                            It should be either TRUE or FALSE")
  geInputGenes<-inputGenes
  inputGenes<-geneIds(inputGenes)
  ##Load gene mapping and orthologs Data
  isHomolog<-FALSE
  e <- new.env()
  load(file = system.file("extdata", "combine-expression.rda", package =
                            "TissueEnrich"),envir = e)
  rnaSeqDatasetGroups <- list(
    "Protein-Atlas", "GTEx-Combine", "ENCODE Dataset"
  )
  if (rnaSeqDataset > length(rnaSeqDatasetGroups))
    stop("Please enter correct dataset id.",call. = FALSE)

  d <- rnaSeqDatasetGroups[[rnaSeqDataset]]
  expressionDataLocal<-e$dataset[[d]]$expressionData
  tissueDetails<-e$dataset[[d]]$tissueDetails
  tissueSpecificGenes<-e$dataset[[d]]$tissueSpecificGenes

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

  groups <- list(
    c("Tissue-Enriched", "Tissue-Enhanced", "Group-Enriched"),
    "Tissue-Enriched", "Tissue-Enhanced", "Group-Enriched"
  )
  if (tissueSpecificGeneType > length(groups))
    stop("Tissue-specific gene type is not correct.",call. = FALSE)

  group <- groups[[tissueSpecificGeneType]]
  tsGenes <- tissueSpecificGenes$Group %in% group
  finalTissueSpecificGenes <- tissueSpecificGenes[tsGenes, , drop=FALSE]

  symbol<-geneIdType(SymbolIdentifier())
  ensembl<-geneIdType(ENSEMBLIdentifier())
  geneId<-list(1, 2)
  names(geneId)<-c(symbol,ensembl)
  inputGeneId<-geneIdType(geneIdType(geInputGenes))
  if (inputGeneId != symbol &&
      inputGeneId != ensembl)
    stop("Gene Id type is not correct.",call. = FALSE)

  geneFormat <- geneId[[inputGeneId]]
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
        inputEnsemblGenes<-geneMapping %>%
          dplyr::filter(Human.gene.name %in% inputGenes) %>%
          dplyr::select(Gene.stable.ID)
      }else
      {
        inputEnsemblGenes<-geneMapping %>%
          dplyr::filter(Human.gene.stable.ID %in% inputGenes) %>%
          dplyr::select(Gene.stable.ID)
      }
      inputEnsemblGenes<-as.character(inputEnsemblGenes$Gene.stable.ID)

      ##Updating the geneMapping
      geneMappingForCurrentDataset<-geneMapping %>%
        dplyr::filter(Gene.stable.ID %in% row.names(expressionDataLocal))

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
        genesNotFound<-c(genesNotFoundList)
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
    ensemblGenesInCurrentDataset<-as.character(geneMappingForCurrentDataset$Gene)
    geneSymbolInCurrentDataset<-as.character(geneMappingForCurrentDataset$Gene.name)
    if(length(inputGenes) != length(intersect(ensemblGenesInCurrentDataset,inputEnsemblGenes)))
    {
      if(geneFormat == 2)
      {
        genesNotFound<-base::setdiff(inputGenes,geneSymbolInCurrentDataset)
      }else
      {
        genesNotFound<-base::setdiff(inputGenes,ensemblGenesInCurrentDataset)
      }
      inputEnsemblGenes<-intersect(ensemblGenesInCurrentDataset,inputEnsemblGenes)
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
      teExpressionData<-left_join(teExpressionData,geneMappingForCurrentDataset, by = "Gene")
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

      ##code to convert ensembl Id for groups
      teInputGeneGroups<-tissueGenes %>% dplyr::filter(Gene %in% intGenes$Gene) %>% select(Gene,Group)
      teInputGeneGroups<-left_join(teInputGeneGroups,geneMappingForCurrentDataset, by = "Gene") %>% select(Gene.name,Group)
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

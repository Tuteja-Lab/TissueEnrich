ensure_numeric<-ensure_that(is.numeric(.),"Please enter numeric values.")


#' Calculation of tissue specific genes by using the algorithm from human protein atlas project.
#'
#' @param expressionData A dataframe object containing gene expression values (Rows are genes and Tissues are columns) .
#' @param foldChangeThreshold A number. Threshold of fold change, default 5.
#' @param maxNumberOfTissues A number. Maximum number of tissues in a group for group enriched genes, default 7.
#' @param expressedGeneThreshold A number. Minimum gene expression cutoff for the gene to be called as expressed, default 1.
#' @export
#' @return A data frame object with three columns, Gene, Tissues, and the enrichment group of the gene in the given tissue.
#' @examples
#' library(tidyverse)
#' data<-read.table("humanProteinAtlasV18.tsv",header = T,sep = "\t")
#' spreadData<-data %>% spread(key = Sample,value= Value)
#' row.names(spreadData)<-spreadData[,1]
#' humanProteinAtlas<-spreadData[gene,c(4:ncol(spreadData))]
#' TSgenes<-tissueSpecificGenesRetrieval(humanProteinAtlas)
#' head(TSgenes)

tissueSpecificGenesRetrieval<-function(expressionData,foldChangeThreshold=5,maxNumberOfTissues=7,expressedGeneThreshold=1)
{
  ###Add checks for the conditions
  expressionData<-ensure_that(expressionData, is.data.frame(.) && !is.null(.) && (nrow(.) > 0) && (ncol(.) > 1),err_desc = "expressionData should be a non-empty dataframe object with atleast 1 gene and 2 tissues. Rows are treated as genes and columns as tissues.")
  foldChangeThreshold<-ensure_that(foldChangeThreshold, is.numeric(.) && (. >=1),err_desc = "foldChangeThreshold should be a numeric value greater than or equal to 1.")
  maxNumberOfTissues<-ensure_that(maxNumberOfTissues, is.numeric(.) && (. >=2),err_desc = "maxNumberOfTissues should be an integer value greater than or equal to 2.")
  expressedGeneThreshold<-ensure_that(expressedGeneThreshold, is.numeric(.) && (. >=0),err_desc = "expressedGeneThreshold should be a numeric value greater than or equal to 0.")

  minNumberOfTissues<-2
  notExpressed<-data.frame()
  expressedInAll<-data.frame()
  mixedGenes<-data.frame()
  tissueEnriched<-data.frame()
  groupEnriched<-data.frame()
  tissueEnhanced<-data.frame()
  for(gene in row.names(expressionData))
  {
    tpm<-data.frame(t(as.matrix(expressionData[gene,])))
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
            meanTPMForGroup<-mean(groupTPM[1:i,])
            highestTPMOutsideGroup<-tpm[i+1,]
            if((meanTPMForGroup/highestTPMOutsideGroup) >= foldChangeThreshold)
            {
              ##To put this genes as group enriched for only the tissues which satisfy the threhold or all the tissues??
              for(j in 1:i)
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
  return(rbind(tissueEnriched,groupEnriched,tissueEnhanced,notExpressed,expressedInAll,mixedGenes))
}


loading <- function(rdata_file)
{
  e <- new.env()
  load(rdata_file, envir = e)
  e
}


#' Calculation of tissue specific gene enrichment using hypergeometric test
#'
#' @param inputGenes A vector containing the input genes.
#' @param dataset An integer describing the dataset to be used for enrichment, 1 for "Human Protein Atlas", 2 for "GTEx combine", 3 for "GTEx sub-tissue", 4 for "Mouse ENCODE". Default 1.
#' @param organism An integer describing the organism, 1 for "Homo Sapiens", 2 for "Mus Musculus". Default 1.
#' @param tissueSpecificGeneType An integer describing type of tissue specific genes to be used, 1 for "All", 2 for "Tissue-Enriched",3 for "Tissue-Enhanced", and 4 for "Group-Enriched". Default 1.
#' @param multiHypoCorrection Flag to carry out multiple hypothesis correction. Default TRUE.
#' @param geneFormat Type of gene symbol, 1 for "Gene Symbol", 2 for "EnsemblId". Default 1.
#' @param isHomolog Flag for usage of orthologus genes, default FALSE.
#' @export
#' @return A data frame object with three columns, Gene, Tissues, and the enrichment group of the gene in the given tissue.
#' @examples
#' library(tidyverse)
#' data<-read.table("humanProteinAtlasV18.tsv",header = T,sep = "\t")
#' spreadData<-data %>% spread(key = Sample,value= Value)
#' row.names(spreadData)<-spreadData[,1]
#' humanProteinAtlas<-spreadData[gene,c(4:ncol(spreadData))]
#' TSgenes<-tissueSpecificGenesRetrieval(humanProteinAtlas)
#' head(TSgenes)

tissueSpecificGeneEnrichment<-function(inputGenes = NULL,
                                       dataset=1,#c("Human Protein Atlas","GTEx combine","GTEx sub-tissue","Mouse ENCODE"),
                                       organism=c("Homo Sapiens","Mus Musculus"),tissueSpecificGeneType=c("All","Tissue-Enriched","Tissue-Enhanced","Group-Enriched"),
                                        geneFormat=c("EnsemblId","Gene Symbol"),isHomolog=FALSE,multiHypoCorrection=TRUE)
{
  ###Add checks for the conditions
  inputGenes<-ensure_that(inputGenes, !is.null(.) && is.vector(.),err_desc = "inputGenes should be a null vector")
  dataset<-ensure_that(dataset, !is.null(.) && is.integer(.) && ((. >=1) || (. <=4)),err_desc = "dataset should be 1, 2, 3 or 4")
  organism<-ensure_that(organism, !is.null(.) && is.integer(.) && ((. >=1) || (. <=2)),err_desc = "organism should be 1 or 2")
  #tissueSpecificGenes<-ensure_that(tissueSpecificGeneType, is.data.frame(.) && !is.null(.) && (nrow(.) > 0) && (ncol(.) == 3),err_desc = "tissueSpecificGenes should be a non-empty dataframe object.")
  geneFormat<-ensure_that(geneFormat,!is.null(.) && is.integer(.) && ((. >=1) || (. <=2)),err_desc = "geneFormat should be 1 or 2")
  tissueSpecificGeneType<-ensure_that(tissueSpecificGeneType,!is.null(.) && is.integer(.)  && ((. >=1) || (. <=4)),err_desc = "tissueSpecificGeneType should be 1, 2, 3 or 4")
  isHomolog<-ensure_that(isHomolog,!is.null(.) && is.logical(.) ,err_desc = "isHomolog should be a TRUE or FALSE")
  multiHypoCorrection<-ensure_that(multiHypoCorrection,!is.null(.) && is.logical(.) ,err_desc = "multiHypoCorrection should be a TRUE or FALSE")
  #foldChangeThreshold<-ensure_that(foldChangeThreshold, is.numeric(.) && (. >=1),err_desc = "foldChangeThreshold should be a numeric value greater than or equal to 1.")
  #maxNumberOfTissues<-ensure_that(maxNumberOfTissues, is.numeric(.) && (. >=2),err_desc = "maxNumberOfTissues should be an integer value greater than or equal to 2.")
  #expressedGeneThreshold<-ensure_that(expressedGeneThreshold, is.numeric(.) && (. >=0),err_desc = "expressedGeneThreshold should be a numeric value greater than or equal to 0.")

  ##Enviroment to load datasets
  env<-NULL
  ##Load gene mapping and orthologs Data
  if(dataset == 1)
  {
    env<-loading("data/proteinAtlasV18-TSEA.Rdata")
  }else if(dataset == 2)
  {
    env<-loading("data/GTEx-combine-TSEA.Rdata")
  }else if(dataset == 3)
  {
    env<-loading("data/GTEx-combine-TSEA.Rdata")
  }else if(dataset == 4)
  {
    env<-loading("data/ENCODE-TSEA.Rdata")
  }else
  {
    print("Please enter correct dataset")
  }

  expressionData<-env$expressionData
  tissueDetails<-env$tissueDetails
  ##Load gene mapping and orthologs Data
  envMapping<-loading("data/geneMappings.Rdata")

  ##Check for organism and homolog to update geneMapping Variable
  if(organism == 1)#"Homo Sapiens")
  {
    geneMapping<-envMapping$humanGeneMapping
  }else if(organism == 2)#"Mus Musculus")
  {
    geneMapping<-envMapping$mouseGeneMapping
  }else
  {
    print("Please enter correct dataset")
  }

  if(isHomolog)
  {
    geneMapping<-envMapping$mouseHumanOrthologs
  }

  ###Update tissueSpecificGene variable based on the group
  if(tissueSpecificGenesType == 1)
  {
    finalTissueSpecificGenes<-env$tissueSpecificGenes
  }else if(tissueSpecificGenesType == 2)
  {
    finalTissueSpecificGenes<-env$tissueSpecificGenes %>% filter(Group == "Tissue-Enriched")
  }else if(tissueSpecificGenesType == 3)
  {
    finalTissueSpecificGenes<-env$tissueSpecificGenes %>% filter(Group == "Tissue-Enhanced")
  }else if(tissueSpecificGenesType == 4)
  {
    finalTissueSpecificGenes<-env$tissueSpecificGenes %>% filter(Group == "Group-Enriched")
  }else
  {
    print("Something is wrong!")
  }

  inputGenes<-toupper(inputGenes)

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
        inputEnsemblGenes<-geneMapping %>% filter(Human.gene.name %in% inputGenes) %>% select(Gene.stable.ID)
      }else
      {
        inputEnsemblGenes<-geneMapping %>% filter(Human.gene.stable.ID %in% inputGenes) %>% select(Gene.stable.ID)
      }
      inputEnsemblGenes<-as.character(inputEnsemblGenes$Gene.stable.ID)

      ##Updating the geneMapping
      geneMappingForCurrentDataset<-geneMapping %>% filter(Gene.stable.ID %in% row.names(expressionData))

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
        inputEnsemblGenes<-intersect(row.names(expressionData),inputEnsemblGenes)
      }
      geneMappingForCurrentDataset<-geneMappingForCurrentDataset %>% select(Gene.stable.ID,Gene.name)
      colnames(geneMappingForCurrentDataset)<-c("Gene","Gene.name")

    }else
    {
      if(geneFormat == 2)
      {
        inputEnsemblGenes<-geneMapping %>% filter(Gene.name %in% inputGenes) %>% select(Human.gene.stable.ID)
      }else
      {
        inputEnsemblGenes<-geneMapping %>% filter(Gene.stable.ID %in% inputGenes) %>% select(Human.gene.stable.ID)
      }
      inputEnsemblGenes<-as.character(inputEnsemblGenes$Human.gene.stable.ID)
      ##Updating the geneMapping
      geneMappingForCurrentDataset<-geneMapping %>% filter(Human.gene.stable.ID %in% row.names(expressionData))
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
        inputEnsemblGenes<-intersect(row.names(expressionData),inputEnsemblGenes)
      }
      geneMappingForCurrentDataset<-geneMappingForCurrentDataset %>% select(Human.gene.stable.ID,Human.gene.name)
      colnames(geneMappingForCurrentDataset)<-c("Gene","Gene.name")
    }

    finalTissueSpecificGenes<-finalTissueSpecificGenes %>% filter(Gene %in% geneMappingForCurrentDataset$Gene)
  }else{

    ####Normal Calculation
    #####Convert the Gene Id########
    print(geneFormat)
    if(geneFormat == 2)
    {
      inputEnsemblGenes<-geneMapping %>% filter(Gene.name %in% inputGenes) %>% select(Gene)
      inputEnsemblGenes<-as.character(inputEnsemblGenes$Gene)
    }else
    {
      inputEnsemblGenes<-inputGenes
    }
    ##Updating the geneMapping
    geneMappingForCurrentDataset<-geneMapping %>% filter(Gene %in% row.names(expressionData))
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
      inputEnsemblGenes<-intersect(row.names(expressionData),inputEnsemblGenes)
    }
  }

  ####Calculate the Hypergeometric P-Value#########
  tissueNames<-as.character(tissueDetails$RName)
  pValueList<-c()
  overlapGenesList<-c()
  overlapTissueGenesList<-list()
  for(tissue in tissueNames)
  {
    tissueGenes<-finalTissueSpecificGenes %>% filter(Tissue==tissue)
    overlapGenes<-length(intersect(tissueGenes$Gene,inputEnsemblGenes))
    overlapTissueGenesList[[tissue]]<-intersect(tissueGenes$Gene,inputEnsemblGenes)
    GenesInTissue<-nrow(tissueGenes)
    pValue<-phyper(overlapGenes-1,GenesInTissue,nrow(geneMappingForCurrentDataset)-GenesInTissue,length(inputEnsemblGenes),lower.tail = F)
    pValueList<-c(pValueList,pValue)
    overlapGenesList<-c(overlapGenesList,overlapGenes)
  }

  if(multiHypoCorrection)
  {
    ##Multiple hypothesis correction
    pValueList<-p.adjust(pValueList,method = "BH")
  }
  pValueList<-(-log10(pValueList))

  output<-data.frame(Tissue=tissueNames,Log10PValue=pValueList,Tissue.Specific.Genes=overlapGenesList)
  output<-output[with(output, order(-Log10PValue)), ]
  return(list(output,overlapTissueGenesList,genesNotFound))
}


#' Calculation of tissue specific gene enrichment using hypergeometric test for custom datasets
#'
#' @param inputGenes A vector containing the input genes.
#' @param tissueSpecificGenes A dataframe object. Output from `tissueSpecificGenesRetrieval` function. Default NULL.
#' @param tissueSpecificGeneType An integer describing type of tissue specific genes to be used, 1 for "All", 2 for "Tissue-Enriched",3 for "Tissue-Enhanced", and 4 for "Group-Enriched". Default 1.
#' @param multiHypoCorrection Flag to carry out multiple hypothesis correction. Default TRUE.
#' @export
#' @return A data frame object with three columns, Gene, Tissues, and the enrichment group of the gene in the given tissue.
#' @examples
#' library(tidyverse)
#' data<-read.table("humanProteinAtlasV18.tsv",header = T,sep = "\t")
#' spreadData<-data %>% spread(key = Sample,value= Value)
#' row.names(spreadData)<-spreadData[,1]
#' humanProteinAtlas<-spreadData[gene,c(4:ncol(spreadData))]
#' TSgenes<-tissueSpecificGenesRetrieval(humanProteinAtlas)
#' head(TSgenes)

tissueSpecificGeneEnrichmentCustom<-function(inputGenes=NULL,tissueSpecificGenes=NULL,
                                       tissueSpecificGeneType= 1,multiHypoCorrection = TRUE)
{
  ###Add checks for the conditions
  inputGenes<-ensure_that(inputGenes, !is.null(.) && is.vector(.) && !is.null(.),err_desc = "inputGenes should be a null vector")
  ##Check for dataset class, rows and should not be NULL
  tissueSpecificGenes<-ensure_that(tissueSpecificGenes,!is.null(.) && is.data.frame(.) && !is.null(.) && (nrow(.) > 0) && (ncol(.) == 3),err_desc = "tissueSpecificGenes should be a non-empty dataframe object.")
  tissueSpecificGeneType<-ensure_that(tissueSpecificGeneType,!is.null(.) && is.integer(.) && ((. >=1) || (. <=4)),err_desc = "tissueSpecificGeneType should be 1, 2, 3 or 4")
  multiHypoCorrection<-ensure_that(multiHypoCorrection,!is.null(.) && is.logical(.),err_desc = "multiHypoCorrection should be 0 or 1")

  if(tissueSpecificGenesType == 1)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes %>% filter(Group %in% c("Tissue-Enriched","Tissue-Enhanced","Group-Enriched"))
  }else if(tissueSpecificGenesType == 2)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes %>% filter(Group == "Tissue-Enriched")
  }else if(tissueSpecificGenesType == 3)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes %>% filter(Group == "Tissue-Enhanced")
  }else if(tissueSpecificGenesType == 4)
  {
    finalTissueSpecificGenes<-tissueSpecificGenes %>% filter(Group == "Group-Enriched")
  }else
  {
    print("Something is wrong!")
  }

  totalGenes<-as.character(unique(tissueSpecificGenes$Gene))
  genesNotFound<-c()
  genesNotFound<-base::setdiff(inputGenes,totalGenes)
  inputEnsemblGenes<-intersect(inputGenes,totalGenes)

  ####Calculate the Hypergeometric P-Value#########
  tissueNames<-as.character(tissueDetails$RName)
  pValueList<-c()
  overlapGenesList<-c()
  overlapTissueGenesList<-list()
  for(tissue in tissueNames)
  {
    tissueGenes<-finalTissueSpecificGenes %>% filter(Tissue==tissue)
    overlapGenes<-length(intersect(tissueGenes$Gene,inputEnsemblGenes))
    overlapTissueGenesList[[tissue]]<-intersect(tissueGenes$Gene,inputEnsemblGenes)
    GenesInTissue<-nrow(tissueGenes)
    pValue<-phyper(overlapGenes-1,GenesInTissue,length(totalGenes)-GenesInTissue,length(inputEnsemblGenes),lower.tail = F)
    pValueList<-c(pValueList,pValue)
    overlapGenesList<-c(overlapGenesList,overlapGenes)
  }
  if(multiHypoCorrection)
  {
    ##Multiple hypothesis correction
    pValueList<-p.adjust(pValueList,method = "BH")
  }
  pValueList<-(-log10(pValueList))
  output<-data.frame(Tissue=tissueNames,Log10PValue=pValueList,Tissue.Specific.Genes=overlapGenesList)
  output<-output[with(output, order(-Log10PValue)), ]
  return(list(output,overlapTissueGenesList,genesNotFound))
}

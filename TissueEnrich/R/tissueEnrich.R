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
  # foldChangeThreshold<-5
  # maxNumberOfTissues<-7
  minNumberOfTissues<-2
  # expressedGeneThreshold<-1
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

#' Calculation of tissue specific gene enrichment using hypergeometric test
#'
#' @param inputGenes A vector containing the input genes.
#' @param tissueSpecificGenes A dataframe object. Output from `tissueSpecificGenesRetrieval` function, default NULL.
#' @param dataset A String describing the dataset to be used for enrichment, Default "Human Protein Atlas". Other sources "GTEx combine","GTEx sub-tissue","Mouse ENCODE". Not required when providing \code{tissueSpecificGenes}.
#' @param organism A String describing the organism, default "Homo Sapiens", other "Mus Musculus". Not required when providing \code{tissueSpecificGenes}.
#' @param tissueSpecificGeneType A String describing type of tissue specific genes to be used, default "All", others "Tissue-Enriched","Tissue-Enhanced","Group-Enriched".
#' @param geneSymbol Type of gene symbol, default "Gene Symbol", other "EnsemblId".
#' @param isHomolog Flag for usage of orthologus genes, default FALSE. Not required when providing \code{tissueSpecificGenes}.
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

tissueSpecificGeneEnrichment<-function(inputGenes,tissueSpecificGenes=NULL,
                                       dataset=c("Human Protein Atlas","GTEx combine","GTEx sub-tissue","Mouse ENCODE"),
                                       organism=c("Homo Sapiens","Mus Musculus"),tissueSpecificGeneType=c("All","Tissue-Enriched","Tissue-Enhanced","Group-Enriched"),
                                        geneSymbol=c("EnsemblId","Gene Symbol"),isHomolog=FALSE)
{

  ###Add checks for the conditions
  inputGenes<-toupper(inputGenes)
  ###Check if it is ortholog comparison or not#######
  if(orthlog)
  {
    ####Calculation when working with orthologs
    #####Convert the Gene Id########
    genesNotFound<-c()
    if(organism == 1)
    {
      if(input$geneFormatId == 2)
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
        if(input$geneFormatId == 2)
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
      if(input$geneFormatId == 2)
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
        if(input$geneFormatId == 2)
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
    print(input$geneFormatId)
    if(input$geneFormatId == 2)
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
      if(input$geneFormatId == 2)
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
  tissueNames<-as.character(tissueDetails$TissueName)
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
    pValueList<-c(pValueList,-log10(pValue))
    overlapGenesList<-c(overlapGenesList,overlapGenes)
  }

  ###Update the global variable
  overlapTissueGenesListGlobal<<-overlapTissueGenesList
  print(length(overlapGenesList))

  output<-data.frame(Tissue=tissueNames,Log10PValue=pValueList,Samples=as.numeric(tissueDetails$Sample),Tissue.Specific.Genes=overlapGenesList)
  output<-output[with(output, order(-Log10PValue)), ]
  return(list(output,genesNotFound))
}

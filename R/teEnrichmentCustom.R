## To Supress Note
utils::globalVariables(c("%>%", ".", "geneIds", "Log10PValue", "SimpleList"))

#' Calculate tissue-specific gene enrichment using the hypergeometric test for
#' custom datasets
#'
#' @description The teEnrichmentCustom function is used to calculate
#' tissue-specific gene enrichment using tissue-specific genes defined using
#' the teGeneRetrieval function.
#' @author Ashish Jain, Geetu Tuteja
#' @param inputGenes An GeneSet object containing the input genes.
#' @param tissueSpecificGenes A SummarizedExperiment object. Output from
#' `teGeneRetrieval` function. Default NULL.
#' @param tissueSpecificGeneType An integer describing the type of
#' tissue-specific genes to be used. 1 for 'All' (default), 2 for
#' 'Tissue-Enriched',3 for 'Tissue-Enhanced', and 4 for 'Group-Enriched'.
#' Default 1.
#' @param multiHypoCorrection Flag to correct P-values for multiple hypothesis
#' using BH method. Default TRUE.
#' @param backgroundGenes A GeneSet object containing the background gene
#' list, organism type ('Homo Sapiens' or 'Mus Musculus'), and gene id
#' identifier (Gene Symbol or ENSEMBL identifier). The input genes must
#' be present in the background gene list. If not provided all the genes
#' will be used as background.
#' @export
#' @return The output is a list with three objects. The first object is the
#' SummarizedExperiment object containing the enrichment results, the second
#' object contains the tissue-specificity information of the tissue-specific
#' genes for genes from the input gene set, and the third is a GeneSet object
#' containing genes that were not identified in the tissue-specific gene data.
#' @examples
#' library(dplyr)
#' data<-system.file('extdata', 'test.expressiondata.txt', package =
#' 'TissueEnrich')
#' expressionData<-read.table(data,header=TRUE,row.names=1,sep='\t')
#' se<-SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),
#' rowData = row.names(expressionData),colData = colnames(expressionData))
#' output<-teGeneRetrieval(se)
#' head(metadata(output)[['TissueSpecificGenes']])
#' genes<-system.file('extdata', 'inputGenesEnsembl.txt', package =
#' 'TissueEnrich')
#' inputGenes<-scan(genes,character())
#' gs<-GeneSet(geneIds=inputGenes)
#' output2<-teEnrichmentCustom(gs,output)
#' #Plotting the P-Values
#' enrichmentOutput<-setNames(data.frame(assay(output2[[1]]),
#' row.names = rowData(output2[[1]])[,1]),
#' colData(output2[[1]])[,1])
#' enrichmentOutput$Tissue<-row.names(enrichmentOutput)
#' ggplot(enrichmentOutput,aes(x=reorder(Tissue,-Log10PValue),y=Log10PValue,
#' label = Tissue.Specific.Genes,fill = Tissue))+
#' geom_bar(stat = 'identity')+
#' labs(x='', y = '-LOG10(P-Value)')+
#' theme_bw()+
#' theme(legend.position='none')+
#' theme(plot.title = element_text(hjust = 0.5,size = 20),axis.title =
#' element_text(size=15))+
#' theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
#' panel.grid.major= element_blank(),panel.grid.minor = element_blank())

teEnrichmentCustom <- function(inputGenes = NULL, tissueSpecificGenes = NULL,
    tissueSpecificGeneType = 1, multiHypoCorrection = TRUE,
    backgroundGenes = NULL) {
    # start.time <- Sys.time()
    ### Add checks for the conditions
    inputGenes <- ensurer::ensure_that(inputGenes, !is.null(.) && (class(.) ==
        "GeneSet") && !is.null(geneIds(.)),
        err_desc = "Please enter correct inputGenes.
                    It should be a Gene set object.")
    ## Check for dataset class, rows and should not be NULL
    tissueSpecificGenes <- ensurer::ensure_that(tissueSpecificGenes,
        !is.null(.) && class(.) == "SummarizedExperiment" &&
        !is.null(assay(.)) &&
        (nrow(assay(.)) > 0) &&
        (ncol(assay(.)) > 1) &&
        (ncol(colData(.)) == 1),
            err_desc = "expressionData should
            be a non-empty SummarizedExperiment
            object with atleast 1 gene and 2
            tissues.")

    tissueSpecificGenes <- setNames(data.frame(assay(tissueSpecificGenes)),
        colData(tissueSpecificGenes)[,1])
    tissueSpecificGenes <- ensurer::ensure_that(tissueSpecificGenes,
        !is.null(.) &&
        is.data.frame(.) &&
        !is.null(.) &&
        (nrow(.) > 0) &&
        (ncol(.) == 3),
            err_desc = "Please enter correct tissueSpecificGenes. It should be
                        a non-empty dataframe object.")
    tissueSpecificGeneType <- ensurer::ensure_that(tissueSpecificGeneType,
        !is.null(.) &&
        is.numeric(.) && . <= 4 || . >= 1, err_desc = "Please enter correct
                                                tissueSpecificGeneType. It
                                                should be 1 for All,
                                                2 for Tissue-Enriched,
                                                3 for Tissue-Enhanced,
                                                and 4 for Group-Enriched")
    multiHypoCorrection <- ensurer::ensure_that(multiHypoCorrection,
        !is.null(.) &&
        is.logical(.), err_desc = "Please enter correct
                                            multiHypoCorrection. It should be
                                            either TRUE or FALSE")
    groups <- list(c("Tissue-Enriched", "Tissue-Enhanced", "Group-Enriched"),
        "Tissue-Enriched", "Tissue-Enhanced", "Group-Enriched")

    if (tissueSpecificGeneType > length(groups))
        stop("Tissue-specific gene type is not correct.", call. = FALSE)

    group <- groups[[tissueSpecificGeneType]]
    ridx <- tissueSpecificGenes$Group %in% group
    finalTissueSpecificGenes <- tissueSpecificGenes[ridx, , drop = FALSE]

    geInputGenes <- inputGenes
    inputGenes <- geneIds(inputGenes)
    inputGenes <- unique(inputGenes)

    backgroundGenes <- ensurer::ensure_that(backgroundGenes,
        is.null(.) || (class(.) == "GeneSet" && !is.null(geneIds(.)) &&
        length(intersect(inputGenes,geneIds(.))) == length(inputGenes)),
        err_desc = "Please enter correct backgroundGenes.
        It should be a Gene set object.
        All input genes must be present in the background list")

    totalGenes <- as.character(unique(tissueSpecificGenes$Gene))
    genesNotFound <- GeneSet(geneIds = base::setdiff(inputGenes, totalGenes))
    inputEnsemblGenes <- intersect(inputGenes, totalGenes)

    bgGenesFlag <- FALSE
    if(!is.null(backgroundGenes))
    {
        bgGenesFlag <- TRUE
        totalGenes <- geneIds(backgroundGenes)
        finalTissueSpecificGenes <- finalTissueSpecificGenes %>%
            dplyr::filter(Gene %in% totalGenes)
    }else
    {
        message("No background list provided. Using all the
                genes as background.")
    }
    #### Calculate the Hypergeometric P-Value#########
    tissueNames <- as.character(unique(finalTissueSpecificGenes$Tissue))
    pValueList <- c()
    x <- lapply(tissueNames, FUN =
                function(tissue){
                    tissueGenes <- finalTissueSpecificGenes %>%
                        dplyr::filter(Tissue == tissue)
                    overlapGenes <- length(intersect(tissueGenes$Gene,
                                                        inputEnsemblGenes))
                    teInputGeneGroups <- tissueGenes %>%
                        dplyr::filter(Gene %in% inputEnsemblGenes) %>%
                            select(Gene, Group)

                    seTeInputGeneGroups <- SummarizedExperiment(
                        assays = SimpleList(as.matrix(teInputGeneGroups)),
                                    rowData = row.names(teInputGeneGroups),
                                    colData = colnames(teInputGeneGroups))

                    nTeGenesInTissue <- nrow(tissueGenes)
                    nTotalGenes <- length(totalGenes)
                    nTotalInputGenes <- length(inputEnsemblGenes)
                    pValue <- stats::phyper(overlapGenes - 1, nTeGenesInTissue,
                                            nTotalGenes - nTeGenesInTissue,
                                            nTotalInputGenes,
                                            lower.tail = FALSE)
                    foldChange <- (overlapGenes/nTotalInputGenes)/
                        (nTeGenesInTissue/nTotalGenes)
                    return(c(pValue,overlapGenes,tissue,
                                seTeInputGeneGroups,foldChange))
                })
    df <- do.call("rbind", x)
    pValueList <- unlist(df[,1])
    ## Multiple hypothesis correction
    if (multiHypoCorrection) {
        pValueList <- stats::p.adjust(pValueList, method = "BH")
    }
    pValueList <- (-log10(pValueList))
    output <- matrix(c(pValueList,unlist(df[,2]),unlist(df[,5])),
        nrow=length(pValueList))
    seOutput <- SummarizedExperiment(assays = SimpleList(output),
                                    rowData = unlist(df[,3]),
                                    colData = c("Log10PValue",
                                                "Tissue.Specific.Genes",
                                                "fold.change"))
    overlapTissueGenesList<-df[,4]
    names(overlapTissueGenesList)<-df[,3]
    # end.time <- Sys.time()
    # time.taken <- end.time - start.time
    # print(time.taken)
    return(list(seOutput, overlapTissueGenesList, genesNotFound))
}

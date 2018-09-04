## To Supress Note
utils::globalVariables(c("dataset", "%>%", "Gene",
    "Gene.name", "Gene.stable.ID", ".", "Human.gene.name",
    "Human.gene.stable.ID", "Group", "Tissue", "geneIds",
    "geneIdType", "ENSEMBLIdentifier", "Log10PValue",
    "SimpleList"))

#' Calculate tissue-specific gene enrichment using the hypergeometric test
#'
#' @description The teEnrichment function is used to calculate the enrichment
#' of tissue-specific genes, given an input gene set. It uses tissue-specific
#' genes defined by processing RNA-Seq datasets from human and mouse.
#' @author Ashish Jain, Geetu Tuteja
#' @param inputGenes A GeneSet object containing the input genes, organism type
#' ('Homo Sapiens' or 'Mus Musculus'), and gene id identifier (Gene Symbol or
#' ENSEMBL identifier).
#' @param rnaSeqDataset An integer describing the dataset to be used for
#' enrichment analysis. 1 for 'Human Protein Atlas' (default), 2 for 'GTEx',
#' 3 for 'Mouse ENCODE'. Default 1.
#' @param tissueSpecificGeneType An integer describing the type of tissue-
#' specific genes to be used. 1 for 'All' (default), 2 for 'Tissue-Enriched',
#'  3 for 'Tissue-Enhanced', and 4 for 'Group-Enriched'. Default 1.
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
#' and the third object contains the expression values and tissue-specificity
#' information of the tissue-specific genes for genes from the input
#' gene set, and the fourth is a GeneSet object containing genes that
#' were not identified in the tissue-specific gene data.
#'
#' @examples
#' library(dplyr)
#' library(ggplot2)
#' genes<-system.file('extdata', 'inputGenes.txt', package = 'TissueEnrich')
#' inputGenes<-scan(genes,character())
#' gs<-GeneSet(geneIds=inputGenes,organism='Homo Sapiens',
#' geneIdType=SymbolIdentifier())
#' output<-teEnrichment(gs)
#' seEnrichmentOutput<-output[[1]]
#' enrichmentOutput<-setNames(data.frame(assay(seEnrichmentOutput),
#' row.names = rowData(seEnrichmentOutput)[,1]),
#' colData(seEnrichmentOutput)[,1])
#' enrichmentOutput$Tissue<-row.names(enrichmentOutput)
#' #Plotting the P-Values
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


teEnrichment <- function(inputGenes = NULL, rnaSeqDataset = 1,
    tissueSpecificGeneType = 1, multiHypoCorrection = TRUE,
    backgroundGenes = NULL) {
    ### Add checks for the conditions
    inputGenes <- ensurer::ensure_that(inputGenes,
        !is.null(.) && (class(.) == "GeneSet") && !is.null(geneIds(.)),
        err_desc = "Please enter correct inputGenes.
                    It should be a Gene set object.")
    rnaSeqDataset <- ensurer::ensure_that(rnaSeqDataset,
        !is.null(.) && is.numeric(.) && . <= 3 || . >=
            1, err_desc = "Please enter correct
                            rnaSeqDataset. It should be 1 for Human
                            Protein Atlas, 2 for GTEx combine, 3 for
                            GTEx sub-tissue, 4 for Mouse ENCODE.")
    tissueSpecificGeneType <- ensurer::ensure_that(tissueSpecificGeneType,
        !is.null(.) && is.numeric(.) && . <= 4 || . >=
            1, err_desc = "Please enter correct
                            tissueSpecificGeneType. It should be
                            1 for All, 2 for Tissue-Enriched,3 for
                            Tissue-Enhanced, and 4 for Group-Enriched")

    multiHypoCorrection <- ensurer::ensure_that(multiHypoCorrection,
        !is.null(.) && is.logical(.),
        err_desc = "Please enter correct multiHypoCorrection.
                    It should be either TRUE or FALSE")

    geInputGenes <- inputGenes
    inputGenes <- geneIds(inputGenes)

    backgroundGenes <- ensurer::ensure_that(backgroundGenes,
        is.null(.) || (class(.) == "GeneSet" && !is.null(geneIds(.)) &&
        length(intersect(inputGenes,geneIds(.))) == length(unique(inputGenes))),
        err_desc = "Please enter correct backgroundGenes.
                    It should be a Gene set object.
                    All input genes must be present in the background list")

    bgGenesFlag <- FALSE
    if(!is.null(backgroundGenes))
    {
        bgGenesFlag <- TRUE
        geBackgroundGenes <- backgroundGenes
        backgroundGenes <- geneIds(geBackgroundGenes)
    }else
    {
        message("No background list provided. Using all the
                genes as background.")
    }
    ## Load the RNA-Seq dataset with mappings
    e <- new.env()
    load(file = system.file("extdata", "combine-expression.rda",
        package = "TissueEnrich"), envir = e)
    rnaSeqDatasetGroups <- list("Protein-Atlas", "GTEx-Combine",
        "ENCODE Dataset")
    if (rnaSeqDataset > length(rnaSeqDatasetGroups))
        stop("Please enter correct dataset id.", call. = FALSE)

    d <- rnaSeqDatasetGroups[[rnaSeqDataset]]
    expressionDataLocal <- e$dataset[[d]]$expressionData
    tissueDetails <- e$dataset[[d]]$tissueDetails
    tissueSpecificGenes <- e$dataset[[d]]$tissueSpecificGenes

    ## Map the organism type
    organismGroups <- list(`Homo Sapiens` = c(1, "humanGeneMapping"),
        `Mus Musculus` = c(2, "mouseGeneMapping"))
    if (is.null(organismGroups[[organism(geInputGenes)]]))
        stop("Please enter correct organism.", call. = FALSE)

    organism <- as.numeric(organismGroups[[organism(geInputGenes)]][1])
    geneMapping <- e$dataset[[organismGroups[[organism(geInputGenes)]][2]]]

    ## Check homolog functionality
    isHomologMapping <- list(c(FALSE, FALSE, TRUE),
        c(TRUE, TRUE, FALSE))
    isHomolog <- isHomologMapping[[organism]][rnaSeqDataset]
    if (isHomolog) {
        geneMapping <- e$dataset$mouseHumanOrthologs
    }
    groups <- list(c("Tissue-Enriched", "Tissue-Enhanced",
        "Group-Enriched"), "Tissue-Enriched", "Tissue-Enhanced",
        "Group-Enriched")
    if (tissueSpecificGeneType > length(groups))
        stop("Tissue-specific gene type is not correct.",
            call. = FALSE)

    group <- groups[[tissueSpecificGeneType]]
    tsGenes <- tissueSpecificGenes$Group %in% group
    finalTissueSpecificGenes <- tissueSpecificGenes[tsGenes,
        , drop = FALSE]

    ## Map gene format Id
    symbol <- geneIdType(SymbolIdentifier())
    ensembl <- geneIdType(ENSEMBLIdentifier())
    geneIdMap <- list(1, 2)
    names(geneIdMap) <- c(ensembl, symbol)
    inputGeneId <- geneIdType(geneIdType(geInputGenes))
    if (is.null(geneIdMap[[inputGeneId]]))
        stop("Gene Id type is not correct.", call. = FALSE)

    geneFormat <- geneIdMap[[inputGeneId]]
    inputGenes <- toupper(inputGenes)
    inputGenes <- unique(inputGenes)
    if(bgGenesFlag)
    {
        backgroundGenes <- toupper(backgroundGenes)
        backgroundGenes <- unique(backgroundGenes)
    }
    # print(isHomolog) Check if it is ortholog
    # comparison or not#######
    if (isHomolog) {
        #### Calculation when working with orthologs Convert
        #### the Gene Id########
        genesNotFound <- c()
        homologMappingHeader <- list(list(c("Human.gene.stable.ID",
            "Gene.stable.ID"), c("Human.gene.name",
            "Gene.stable.ID")), list(c("Gene.stable.ID",
            "Human.gene.stable.ID"), c("Gene.name",
            "Human.gene.stable.ID")))
        mappingList <- homologMappingHeader[[organism]][[geneFormat]]
        geneId <- geneMapping[, mappingList[2]] %in%
            row.names(expressionDataLocal)
        geneMappingForCurrentDataset <- geneMapping[geneId,]
        genesNotFound <- base::setdiff(inputGenes,
                            geneMappingForCurrentDataset[, mappingList[1]])
        ##Change the mapping genes based on background genes
        if(bgGenesFlag)
        {
            geneId <- geneMappingForCurrentDataset[, mappingList[1]] %in%
                backgroundGenes
            geneMappingForCurrentDataset <- geneMappingForCurrentDataset[
                                            geneId, ]
        }
        geneId <- geneMappingForCurrentDataset[, mappingList[1]] %in%
            inputGenes
        inputEnsemblGenes <- geneMappingForCurrentDataset[geneId,
            mappingList[2]]
        geneMappingForCurrentDataset <- geneMappingForCurrentDataset[,
            rev(mappingList)]
        colnames(geneMappingForCurrentDataset) <- c("Gene", "Gene.name")
        finalTissueSpecificGenes <- finalTissueSpecificGenes %>%
            dplyr::filter(Gene %in% geneMappingForCurrentDataset$Gene)
    } else {

        geneIdHeaderMapping <- list(c("Gene", "Gene"),
            c("Gene", "Gene.name"))
        #### Normal Calculation
        mappingHeader <- geneIdHeaderMapping[[geneFormat]]
        ## Updating the geneMapping
        geneMappingForCurrentDataset <- geneMapping %>%
            dplyr::filter(Gene %in% row.names(expressionDataLocal))
        ##### Check for genes which are not there in our
        ##### list########
        genesNotFound <- c()
        genesNotFound <- base::setdiff(inputGenes,
            geneMappingForCurrentDataset[, mappingHeader[2]])

        ##Change the mapping genes based on background genes
        if(bgGenesFlag)
        {
            geneId <- geneMappingForCurrentDataset[, mappingHeader[2]] %in%
                backgroundGenes
            geneMappingForCurrentDataset <- geneMappingForCurrentDataset[
                geneId, ]
            finalTissueSpecificGenes <- finalTissueSpecificGenes %>%
                dplyr::filter(Gene %in% geneMappingForCurrentDataset$Gene)
        }
        ##### Convert the Gene Id########
        geneId <- geneMappingForCurrentDataset[, mappingHeader[2]] %in%
            inputGenes
        inputEnsemblGenes <- geneMappingForCurrentDataset[geneId,
            mappingHeader[1]]
    }

    #### Calculate the Hypergeometric P-Value#########
    tissueNames <- as.character(tissueDetails$RName)
    print(paste0("Total background genes:",nrow(geneMappingForCurrentDataset)))
    print(paste0("Total input genes:",length(inputEnsemblGenes)))
    x <- lapply(tissueNames, FUN = function(tissue) {

        tissueGenes <- finalTissueSpecificGenes %>%
            dplyr::filter(Tissue == tissue)
        overlapGenes <- length(intersect(tissueGenes$Gene,
            inputEnsemblGenes))
        ## Code block to take the mean of the genes with
        ## multiple ensembl Ids.
        if (geneFormat == 2 || isHomolog) {
            ## code to convert ensembl Id to gene names
            intGenes <- geneMappingForCurrentDataset %>%
                dplyr::filter(Gene %in% intersect(tissueGenes$Gene,
                inputEnsemblGenes))
            teExpressionData <- expressionDataLocal[
                as.character(intGenes$Gene), ]
            genes <- row.names(teExpressionData)
            teExpressionData$Gene <- genes
            teExpressionData <- dplyr::left_join(teExpressionData,
                geneMappingForCurrentDataset, by = "Gene")
            if (nrow(teExpressionData) > 0) {
                ##### Code to take the mean of the genes with multiple
                ##### ensembl Ids.
                teExpressionData <- as.data.frame(
                    t(vapply(unique(teExpressionData$Gene.name),
                    FUN = function(col)
                        colMeans(teExpressionData[
                            teExpressionData$Gene.name == col,
                            c(seq(1, (ncol(teExpressionData) - 2)))]),
                                numeric(length(tissueNames)))))
            } else {
                teExpressionData <- teExpressionData[,
                    seq(1, (ncol(teExpressionData) - 2))]
            }

            ## code to convert ensembl Id for groups
            teInputGeneGroups <- tissueGenes %>% dplyr::filter(Gene %in%
                intGenes$Gene) %>% select(Gene, Group)
            teInputGeneGroups <- dplyr::left_join(teInputGeneGroups,
                geneMappingForCurrentDataset, by = "Gene") %>%
                select(Gene.name, Group)
            colnames(teInputGeneGroups) <- c("Gene",
                "Group")
            ## Code to remove Gene Symbol with two or more
            ## Ensmbel Ids
            teInputGeneGroups <- unique(teInputGeneGroups)
        } else {
            teExpressionData <- expressionDataLocal[intersect(tissueGenes$Gene,
                inputEnsemblGenes), ]
            teInputGeneGroups <- tissueGenes %>% dplyr::filter(Gene %in%
                inputEnsemblGenes) %>% select(Gene,
                Group)
        }

        ## Update the names of tissues with special
        ## characters
        colnames(teExpressionData) <- tissueDetails$TissueName
        teExpressionData <- log2(teExpressionData +
            1)
        seTeExpressionData <- SummarizedExperiment(
            assays = SimpleList(as.matrix(teExpressionData)),
            rowData = row.names(teExpressionData),
            colData = colnames(teExpressionData))

        seTeInputGeneGroups <- SummarizedExperiment(
            assays = SimpleList(as.matrix(teInputGeneGroups)),
            rowData = row.names(teInputGeneGroups),
            colData = colnames(teInputGeneGroups))
        tissueName <- tissueDetails[tissueDetails$RName ==
            tissue, "TissueName"]
        nTeGenesInTissue <- nrow(tissueGenes)
        nTotalGenes <- nrow(geneMappingForCurrentDataset)
        nTotalInputGenes <- length(inputEnsemblGenes)
        print(paste0("Total TS genes in ",tissue," ",nTeGenesInTissue))
        pValue <- stats::phyper(overlapGenes - 1, nTeGenesInTissue,
            nTotalGenes - nTeGenesInTissue,
            nTotalInputGenes, lower.tail = FALSE)
        foldChange <- (overlapGenes/nTotalInputGenes)/
            (nTeGenesInTissue/nTotalGenes)
        return(c(pValue, overlapGenes, tissueName,
            seTeExpressionData, seTeInputGeneGroups,foldChange))
    })

    df <- do.call("rbind", x)
    pValueList <- unlist(df[, 1])
    if (multiHypoCorrection) {
        ## Multiple hypothesis correction
        pValueList <- stats::p.adjust(pValueList, method = "BH")
    }
    pValueList <- (-log10(pValueList))
    output <- data.frame(Log10PValue=pValueList,
                            Tissue.Specific.Genes=unlist(df[, 2]),
                                fold.change=unlist(df[, 6]),
                                row.names = unlist(df[, 3]))
    seOutput <- SummarizedExperiment(assays = SimpleList(as.matrix(output)),
        rowData = row.names(output), colData = colnames(output))
    seTeExpressionData <- df[, 4]
    names(seTeExpressionData) <- df[, 3]
    seTeInputGeneGroups <- df[, 5]
    names(seTeInputGeneGroups) <- df[, 3]
    ## Not working with R 3.5
    # overlapTissueGenesList <- do.call(function(...) mapply(c,
    #     ..., SIMPLIFY = FALSE), args = list(df[, 4],
    #     df[, 5]))
    # overlapTissueGenesList <- apply(cbind(df[,4], df[,5]),1,
    #                                 function(x) unname(unlist(x)))
    # names(overlapTissueGenesList) <- df[, 3]
    return(list(seOutput, seTeExpressionData, seTeInputGeneGroups,
                GeneSet(geneIds = genesNotFound)))
}

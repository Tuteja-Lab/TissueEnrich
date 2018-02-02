#' Processed data for GTEx, Human Protein Atlas, and Mouse ENCODE
#'
#' A list object containing the processed data from GTEx, Human Protein Atlas, and Mouse ENCODE
#' containing the expression values, calculated tissue specific genes, and information about
#' the tissues.
#' @author Ashish Jain, Geetu Tuteja
#' @format A list object containing 7 objects:
#' \describe{
#'   \item{Protein-Atlas}{Human Protein atlas processed data with normalized expression values (in 35 tissues), tissue details, and tissue specific genes.}
#'   \item{GTEx-Combine}{GTEx-combine processed data with normalized expression values (in 29 tissues), tissue details, and tissue specific genes.}
#'   \item{GTEx-Subtissues}{GTEx-subtissues processed data with normalized expression values (in 50 tissues), tissue details, and tissue specific genes.}
#'   \item{ENCODE Dataset}{Human Protein atlas processed data with normalized expression values (in 17 tissues), tissue details, and tissue specific genes.}
#'   \item{mouseGeneMapping}{Ensembl Id and Gene symbol mapping dataframe for mouse protein coding genes.}
#'   \item{mouseHumanOrthologs}{One to one human and mouse orthologs mapping.}
#'   \item{humanGeneMapping}{Ensembl Id and Gene symbol mapping dataframe for human protein coding genes.}
#' }
#'
#' @references
"dataset"

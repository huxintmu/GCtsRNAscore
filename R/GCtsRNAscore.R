#' Calculate Gastric Cancer tsRNA Score (GCtsRNAscore)
#'
#' This function calculates the GC tsRNA scores based on a given gene expression matrix.
#' The expression matrix should have rows as samples and columns as genes.
#' If the matrix is missing required genes, an error will be thrown.
#'
#' @param my_expr A numeric matrix with rows as samples and columns as genes.
#' @return A data frame with sample names and their tsRNA scores.
#' @export
predict_tsRNA <- function(my_expr) {
  required_genes <- c("AKAP12","APOD","AXL","CDH11","CMTM3",
                      "CPE","CTGF","CXCR4","EDNRA","EFEMP1",
                      "FBN1","FNDC1","HTRA3","IGFBP7","LHFP",
                      "LPPR4","LUM","MEOX2","OLFML2B","PRSS23",
                      "SERPINF1","SLC22A17","THBS1","TIMP2","VGLL3")
  missing_genes <- required_genes[!required_genes %in% colnames(my_expr)]
  if (length(missing_genes) > 0) {
    stop("Error: The following genes are missing in the expression matrix and are required for GCtsRNAscore calculation: ",
         paste(missing_genes, collapse = ", "), ". Please provide a complete expression matrix.")
  }
  load(system.file("data/fit.rdata", package = "GCtsRNAscore"))
  GCtsRNAscore <- cbind(rownames(my_expr), tsRNAscore=predict(fit, newdata = my_expr)$predicted)
  return(GCtsRNAscore)
}

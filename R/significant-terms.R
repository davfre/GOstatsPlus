#' @title significant_terms
#'
#' @description
#' \code{significant_terms} extracts the significant GO terms from a GOstat
#' results object
#'
#' @details
#' Details of function go here
significant_terms<-function(GO_OBJ, cutoff){
  terms_pvals = pvalues(GO_OBJ)
  return(terms_pvals[terms_pvals<cutoff])
}

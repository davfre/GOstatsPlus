#' significant terms
#'
#' Extract the significant GO terms from a GOstat results object

significant_terms<-function(GO_OBJ, cutoff){
  terms_pvals = pvalues(GO_OBJ)
  return(terms_pvals[terms_pvals<cutoff])
}

#' @title test_GO
#'
#' @description
#' \code{test_GO} returns GO terms overrepresented in a gene set
#'
#' @details
#' Use GO annotations and GOstats in R to quickly test for over-represented GO terms

test_GO<-function(genes, ontology, gsc, direction = "over", universe, pval = 0.01){
  Test<-NA

  # default to test against all genes with any annotation
  if(missing(universe)){
    universe<-unique(unlist(geneIds(gsc)))
  }
  
  if(length(genes[genes %in% universe])>0){
		params <- GSEAGOHyperGParams(name = "GSEA based annotation Parameters",
                                 geneSetCollection = gsc,
                                 geneIds = genes,
                                 universeGeneIds = universe,
                                 ontology = ontology,
                                 pvalueCutoff = pval,
                                 conditional = FALSE,
                                 testDirection = direction)
 		Test <- hyperGTest(params)
 	}
	return(Test)
}


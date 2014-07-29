## ------------------------------------------------------------------------
library(GOstats)
library(GSEABase)
library(GOstatsPlus)

## ------------------------------------------------------------------------
if(file.exists("organism_gsc.rda")){
  load(file="organism_gsc.rda")
}else{
  fn_annot = "example_GO_b2g.annot"
  GO_gsc = b2g_to_gsc(file = fn_annot, organism = "myOrganism")
  save(GO_gsc, file = "organism_gsc.rda")
}


## ------------------------------------------------------------------------
gene_IDs_of_interest = head(unique(unlist(geneIds(GO_gsc))),500)

GO_BP = test_GO(gene_IDs_of_interest, ontology = "BP", gsc=GO_gsc, pval = 0.01)
GO_MF = test_GO(gene_IDs_of_interest, ontology = "MF", gsc=GO_gsc, pval = 0.01)
GO_CC = test_GO(gene_IDs_of_interest, ontology = "CC", gsc=GO_gsc, pval = 0.01)

## ------------------------------------------------------------------------
head(summary(GO_BP))
head(summary(GO_MF))
head(summary(GO_CC))

## ------------------------------------------------------------------------
go_terms = unlist(lapply(c(GO_BP, GO_MF, GO_CC),
                         function(d){significant_terms(d, cutoff = 0.01)}
                  ))
write.table(go_terms, file = "go_terms.significant.txt",quote=FALSE)


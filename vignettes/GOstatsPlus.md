Use GOstats to test for over-represented GO terms
========================================================


```r
library(GOstats)
```

```
## Loading required package: Biobase
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Loading required package: Category
## Loading required package: AnnotationDbi
## Loading required package: Matrix
## Loading required package: GO.db
## Loading required package: DBI
## 
## Loading required package: graph
## 
## Attaching package: 'GOstats'
## 
## The following object is masked from 'package:AnnotationDbi':
## 
##     makeGOGraph
```

```r
library(GSEABase)
```

```
## Loading required package: annotate
```

```r
library(GOstatsPlus)
```

### Create gene set
Create a GOstat gene set collection (gsc) object based on a text file having
the format gene_id\tGO_id\tDescription\n (e.g .annot file export from Blast2GO)

This is a bit time consuming, so do save the R object for later reuse


```r
if(file.exists("data/organism_gsc.rda")){
  load(file="organism_gsc.rda")
}else{
  fn_annot <- system.file("extdata", "example_GO_b2g.annot", package="GOstatPlus")
  GO_gsc = b2g_to_gsc(file=fn_annot,organism="myOrganism")
  save(GO_gsc,file="organism_gsc.rda")
}
```

### Perform test
Genes of interest vs all genes that have any GO annotation (default)
Perform tests for over-represented (default) biological process, molecular function and cellular component terms:


```r
gene_IDs_of_interest = head(unique(unlist(geneIds(GO_gsc))),500)

GO_BP<-test_GO(gene_IDs_of_interest, ontology="BP", gsc=GO_gsc, pval = 0.01)
GO_MF<-test_GO(gene_IDs_of_interest, ontology="MF", gsc=GO_gsc, pval = 0.01)
GO_CC<-test_GO(gene_IDs_of_interest, ontology="CC", gsc=GO_gsc, pval = 0.01)
```

### View results

```r
head(summary(GO_BP))
```

```
##       GOBPID     Pvalue OddsRatio ExpCount Count Size
## 1 GO:0044702  0.000e+00     249.7   11.283   266  330
## 2 GO:0022414  0.000e+00     784.1   17.027   413  498
## 3 GO:0032504  0.000e+00     275.2   11.283   270  330
## 4 GO:0000003  0.000e+00    4329.5   19.899   484  582
## 5 GO:0048609  0.000e+00     273.0   11.249   269  329
## 6 GO:0019953 1.157e-306     240.8    9.471   228  277
##                                            Term
## 1          single organism reproductive process
## 2                          reproductive process
## 3           multicellular organism reproduction
## 4                                  reproduction
## 5 multicellular organismal reproductive process
## 6                           sexual reproduction
```

```r
head(summary(GO_MF))
```

```
##       GOMFID    Pvalue OddsRatio ExpCount Count Size
## 1 GO:0005515 3.462e-43     3.860  206.757   353 7269
## 2 GO:0043566 4.896e-30    14.336    4.295    42  151
## 3 GO:0043168 1.394e-28     3.190   73.953   171 2600
## 4 GO:0032559 3.814e-28     3.671   45.737   128 1608
## 5 GO:0005524 3.995e-28     3.686   45.140   127 1587
## 6 GO:0030554 4.299e-28     3.665   45.794   128 1610
##                             Term
## 1                protein binding
## 2 structure-specific DNA binding
## 3                  anion binding
## 4  adenyl ribonucleotide binding
## 5                    ATP binding
## 6      adenyl nucleotide binding
```

```r
head(summary(GO_CC))
```

```
##       GOCCID    Pvalue OddsRatio ExpCount Count Size
## 1 GO:0044464 8.677e-31     9.927    382.2   464 9516
## 2 GO:0005623 1.038e-30     9.906    382.3   464 9520
## 3 GO:0044424 7.116e-28     3.937    322.9   423 8041
## 4 GO:0043231 5.843e-26     2.758    214.5   326 5340
## 5 GO:0043227 8.360e-26     2.749    214.8   326 5349
## 6 GO:0043226 1.859e-24     2.819    256.6   362 6388
##                                       Term
## 1                                cell part
## 2                                     cell
## 3                       intracellular part
## 4 intracellular membrane-bounded organelle
## 5               membrane-bounded organelle
## 6                                organelle
```

###Export list of significant term, e.g. for import/visualization in REViGO (http://revigo.irb.hr)

```r
go_terms = unlist(lapply(c(GO_BP,GO_MF,GO_CC),function(d){significant_terms(d,cutoff=0.01)}))
write.table(go_terms,file="go_terms.significant.txt",quote=FALSE)
```

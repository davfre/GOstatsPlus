#' @title b2g_to_gsc
#'
#' @description
#' \code{b2g_to_gsc} converts a tab-delimited text file to a GOstat gsc object
#'
#' @details
#' Input file format: geneId GOId geneDesc (tab separated)

b2g_to_gsc = function(file, organism="my-organism"){
  blast2GOannot = read.table(file=file, sep="\t", fill=TRUE)
  colnames(blast2GOannot)<-c("Seq","GO","desc")

  # Create GO package
  GOdata<-data.frame(GO.id=blast2GOannot$GO, evidence.code=rep("ISA",dim(blast2GOannot)[1]), gene.id=as.character(blast2GOannot$Seq))
  goFrame = GOFrame(GOdata, organism = organism)
  goAllFrame = GOAllFrame(goFrame)

  GO_gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  #  GO_mappings<-getGOFrameData(goAllFrame)

  return(GO_gsc)
}

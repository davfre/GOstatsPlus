#' Blast2GO to GSC
#'
#' converts a tab-delimited text file to a GOstat gsc object
#'
#' @param file Tab separated text: gene_id GO_accession GO_description
#' @param organism The name of your organism
#' 
#' @return A GSEA GeneSetCollection object
#'
#' @keywords GO, GSEA, Blast2GO
#'   
#' @export
#' 
#' @examples
#' GO_gsc = b2g_to_gsc(file="myBlast2GO.annot",organism="myOrganism")

b2g_to_gsc = function(file, organism="my-organism"){
  blast2GOannot = read.table(file=file, sep="\t", fill=TRUE)
  colnames(blast2GOannot)<-c("Seq", "GO", "desc")

  # Create GO package
  GOdata<-data.frame(GO.id = blast2GOannot$GO,
                     evidence.code = rep("ISA", dim(blast2GOannot)[1]),
                     gene.id = as.character(blast2GOannot$Seq))
  goFrame = GOFrame(GOdata, organism = organism)
  goAllFrame = GOAllFrame(goFrame)

  GO_gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  #  GO_mappings<-getGOFrameData(goAllFrame)

  return(GO_gsc)
}

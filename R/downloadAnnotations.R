#' Download transcriptome annotatin from Ensemblin a format that customProDB can use
#' This Uses the import function from CustomProDB in combination with BioMart and can take a few hours.
#' @param annotation_path Path to a folder where Annotation files should be stored.
#' @import customProDB
#' @export


downloadAnnotations <- function(annotation_path, archive =FALSE, host ="dec2016.archive.ensembl.org"){
  # If you want to reproduce the results set archive to true, otherwhise it will download the current ensembl version
  ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
                     host=host, path="/biomart/martservice",
                     archive=archive)


  PrepareAnnotationEnsembl(mart=ensembl, annotation_path=annotation_path,
                           splice_matrix=TRUE, dbsnp=NULL,
                           transcript_ids=NULL, COSMIC=FALSE)

}


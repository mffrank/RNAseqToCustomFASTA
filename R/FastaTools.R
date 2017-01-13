
#' Download uniprot ids from BioMart
#' @details internal function
#' @param ids Id mapping from annotation folder
#' @import biomaRt
#' @import UniProt.ws

getBioMart <- function(ids){

  # ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
  #                    host="oct2016.archive.ensembl.org", path="/biomart/martservice",
  #                    archive=FALSE)
  # dict_c <- getBM(mart = ensembl, values = ids$pro_name, filters = "ensembl_peptide_id",
  #               attributes = c("ensembl_peptide_id","uniprot_swissprot","uniprot_sptrembl","uniprot_genename"))

  up <- UniProt.ws(taxId=9606)
  dict <-UniProt.ws::select(up, keys = ids$pro_name,
                            keytype = "ENSEMBL_PROTEIN", columns = c("ENSEMBL_PROTEIN", "UNIPROTKB"))

  dict_2 <- read.delim("~/Master_project_local_copy/uniprot_all_ENSEMBL.tab")
  dict <- merge(ids[ids$pro_name != "",c("tx_name","pro_name")], dict2[,c("")])
  dict <- sapply(ids[ids$pro_name != "","tx_name"], function(x){
    c(x, as.character(dict2[grep(x, dict2$Cross.reference..Ensembl.),"Entry"]))
  })
}
convertExpressionHeaders <- function(headers){

}

#' @title Concatenate all fastas in a Folder.
#' @details Can be used to produce a single Fasta ready for library generation.
#' File is saved in the path directory.
#' @param path Path to file containing folder
#' @param pattern Pattern used to grep filenames (default: ".fa")
#' @param outfile Name of the combined file.
#' @export

concatenateFastas <- function(path, pattern = "\\.fa", outfile = "combined.fasta ..."){
  files <- list.files(path = path, pattern = pattern)
  files <- files[! files == outfile]
  message("Combining ", paste("\n",files), "\n into ", outfile)
  allfasta <-  do.call("c",lapply(files, function(x) readLines(file(paste0(path,"/",x)))))
  message("Writing output...")
  write(allfasta, file = paste0(path,outfile))
}


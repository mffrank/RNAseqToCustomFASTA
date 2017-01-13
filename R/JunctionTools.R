
#' Generate a Junction Summary Report
#' @details Saves a table of the types of the Splice junctions to working directory
#' @param sharedjun GRanges object containing all splice junctions
#' @param splicemax Splice Matrix from annotation folder
#' @param ids Id mapping from annotation folder
#' @param txdb Txdb object from annotation folder
#' @param print boolean, should the result be printed in working directory as pdf
#' or be returned as a list.
#' @export

junctionSummaryReport <- function(junction_type, save = TRUE){

  output <- list(junction_type = table(junction_type[, 'jun_type']))
  if(save){
    saveRDS(output,"VariantSummary.rda")
    # pdf("Variant_summary_report.pdf", width = 6, height = 4)
    # for(i in 1:length(output)){
    #   print(output[[i]])
    # }
    # dev.off()
  } else{
    return(output)
  }
}

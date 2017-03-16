
#' Generate a Junction Summary Report
#' @details Saves a table of the types of the Splice junctions to working directory
#' @param sharedjun GRanges object containing all splice junctions
#' @param splicemax Splice Matrix from annotation folder
#' @param ids Id mapping from annotation folder
#' @param txdb Txdb object from annotation folder
#' @param print boolean, should the result be printed in working directory as pdf
#' or be returned as a list.
#' @export

filterSpacePlot <- function(covfilter_unique = c(1,2,3,4,5,10), share_num = c(1,2), bed_directory, splicemax, txdb, ids){

  combinations <- expand.grid(covfilter_unique, share_num)
  colnames(combinations) <- c("covfilter_unique", "share_num")
  combinations$remaining_junctions <- 0
  juncs_list <- list()
  message("Importing Junctions with different filters \n \n")
  for(i in 1:nrow(combinations)){
    juncs <- getStarJunctions(bed_directory, splicemax = splicemax, txdb = txdb, ids = ids,
                              pattern = "\\.tab", share_num = combinations[i,2], extend = 100, skip = 0,
                              covfilter_unique = combinations[i,1], covfilter_multi = 0, overhang = "max_overhang")
    combinations[i,3] <- length(juncs)
    juncs_list[[i]] <- juncs
  }

  m <- ggplot(combinations, aes(x=covfilter_unique, y = remaining_junctions, colour = share_num)) +
    geom_point()

  return(m)
}
#' Generate a Junction Summary Report
#' @details Saves a table of the types of the Splice junctions to working directory
#' @param sharedjun GRanges object containing all splice junctions
#' @param splicemax Splice Matrix from annotation folder
#' @param ids Id mapping from annotation folder
#' @param txdb Txdb object from annotation folder
#' @param print boolean, should the result be printed in working directory as pdf
#' or be returned as a list.
#' @export

junctionTypeDistribution <- function(covfilter_unique = c(1,2,5), share_num = c(2), bed_directory, splicemax, txdb, ids){

  combinations <- expand.grid(covfilter_unique, share_num)
  colnames(combinations) <- c("covfilter_unique", "share_num")

  juncs_list <- list()
  message("Importing Junctions with different filters \n \n")
  for(i in 1:nrow(combinations)){
    juncs <- getStarJunctions(bed_directory, splicemax = splicemax, txdb = txdb, ids = ids,
                              pattern = "\\.tab", share_num = combinations[i,2], extend = 100, skip = 0,
                              covfilter_unique = combinations[i,1], covfilter_multi = 0, overhang = "max_overhang")
    juncs_list[[i]] <- juncs
  }
  summary_list <- list()
  for(i in 1: length(juncs_list)){
    junction_type <- JunctionType(juncs_list[[i]], splicemax, txdb, ids)
    junc_summary <- junctionSummaryReport(junction_type, save = F)
    summary_list[[i]] <- junc_summary
  }


  m <- ggplot(combinations, aes(x=covfilter_unique, y = remaining_junctions, colour = share_num)) +
    geom_point()

  return(m)
}

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

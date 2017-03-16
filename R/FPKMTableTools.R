#' internal, makes named vector dict
#' @param ids query ids
#' @param mapping_table table with ids of both types
#' @param query_id_type name of the query column
#' @param map_id_type name of the map_id column


makeDict <- function(ids, mapping_table, query_id_type, map_id_type){
  if(!(any(map_id_type %in% colnames(mapping_table))) | !(any(map_id_type %in% colnames(mapping_table)))){
    stop("id types must be in colnames of the mapping table")
  }
  dict_tab <- subset(mapping_table, mapping_table[,query_id_type] %in% ids)
  dict <- dict_tab[,map_id_type]
  names(dict) <- dict_tab[,query_id_type]
  return(dict)
}
#' Output Histogram of FPKM values
#' @details Plots a Histogram with the 30% quantile as a vertical line
#' @param fpkm_table Table with FPKM vaules of (multiple) samples
#' @param mapping_table table with ids of both types
#' @import ggplot2
#' @export

FPKMhist <- function(fpkm_table, qt){
  fpkm_table$average_FPKM <- apply(fpkm_table, 1,mean)
  fpkm_table_measured <- fpkm_table[fpkm_table$average > 0,]
  # qt <- c(0.2,0.22,0.25,0.3)
  cutoffs <- data.frame(Quantile = factor(qt), val = quantile(fpkm_table_measured$average_FPKM,qt))

  # Plot FPKM histogram
  m <- ggplot2::ggplot(fpkm_table[fpkm_table$average_FPKM > 10^-5,], aes(x = average_FPKM))
  m <- m + geom_histogram(aes(y = ..density..),binwidth = 0.2) + scale_x_log10()
  m <- m + geom_density(colour = "black",lwd = 1.0)
  m <- m + geom_vline(data = cutoffs,aes(xintercept = val, lty = Quantile))
  return(m)
}

#' Output Histogram of FPKM values for protein-coding vs. non protein coding
#' @details Plots a Histogram with the 30% quantile as a vertical line
#' @param fpkm_table Table with FPKM vaules of (multiple) samples
#' @param id_table table with ids of both types
#' @import ggplot2
#' @export

FPKMhistProteincoding <- function(fpkm_table, id_table, qt){
  fpkm_table$average_FPKM <- apply(fpkm_table, 1,mean)
  fpkm_table_measured <- fpkm_table[fpkm_table$average > 0,]
  dict <- makeDict(rownames(fpkm_table), id_table, "tx_name", "pro_name")
  fpkm_table$protein_id <- dict[rownames(fpkm_table)]
  fpkm_table$protein_coding <- !fpkm_table$protein_id == ""
  # qt <- c(0.2,0.22,0.25,0.3)
  cutoffs <- data.frame(Quantile = factor(qt), val = quantile(fpkm_table_measured$average_FPKM,qt))
  # Plot FPKM histogram
  m <- ggplot2::ggplot(fpkm_table[fpkm_table$average_FPKM > 10^-5,], aes(average_FPKM))
  m <- m + scale_x_log10()
  m <- m + geom_histogram(position='dodge',aes(fill = protein_coding, colour = protein_coding))
  m <- m + geom_vline(data = cutoffs,aes(xintercept = val, lty = Quantile))
  return(m)
}

#' Scatterplot of sd vs log transformed expression
#' @details Log transforms expression, calculates sd and plots it over expression dimension
#' @param fpkm_table Table with FPKM vaules of (multiple) samples
#' @import ggplot2
#' @import vsn
#' @export

plotSDvsMeanExpression <- function(fpkm_table){
  fpkm_table_f <- filterDetectedTranscripts(fpkm_table, quiet = TRUE)
  # fpkm_table_f$log_average_FPKM <- apply(fpkm_table_f, 1,mean)
  # fpkm_table_f <- log10(fpkm_table_f)
  m <- meanSdPlot(log2(as.matrix(fpkm_table_f)),plot = F)
  # fpkm_table_f$sd_log_FPKM <- apply(fpkm_table_f,1,function(x) sd(x[-length(x)])) # /x["average_FPKM"]
  # m <- ggplot2::ggplot(data = fpkm_table_f[fpkm_table_f$log_average_FPKM > -5,])
  # m <- m + geom_point(aes(x = log_average_FPKM, y= sd_log_FPKM),alpha = 0.2)
  return(m)
}


#' Print Pdf with plots to help evaluate FPKM distribution
#' @details Plots a Histograms and sd vs. FPKM
#' @param fpkm_table Table with FPKM vaules of (multiple) samples
#' @param id_table table with ids of both types
#' @param print boolean, should the result be printed in working directory as pdf
#'  or returned as list of ggplot objects
#' @export

FPKMsummaryReport <- function(fpkm_table, id_table, qt = c(0.2,0.22,0.25,0.3), print = TRUE){
  a <- FPKMhist(fpkm_table = fpkm_table, qt)
  b <- FPKMhistProteincoding(fpkm_table = fpkm_table, id_table = id_table, qt)
  c <- plotSDvsMeanExpression(fpkm_table = fpkm_table)
  output <- list(p1 = a, p2 = b, p3 = c)
  if(print){
    pdf("FPKM_summary_plots.pdf", width = 6, height = 4)
    for(i in 1:length(output)){
      print(output[[i]])
    }
    dev.off()
  } else{
    return(output)
  }
}

#' Convert FPKM Table with transcript ids to protein ids
#' @details Matches rownames of fpkm table to protein ids with a mapping table (has to be supplied)
#' @details Rows with no matching protein id are discarded
#' @param fpkm_table Table with FPKM vaules of (multiple) samples
#' @param id_table table with ids of both types
#' @param transcript_id_name colname of the transcript id column in mapping table
#' @param protein_id_name colname of the protein id column in mapping table
#' @export

getProteinCoding <- function(fpkm_table, id_table, transcript_id_name = "tx_name", protein_id_name = "pro_name"){
  dict <- makeDict(rownames(fpkm_table), id_table, transcript_id_name, protein_id_name)
  # fpkm_table$protein_id <- dict[rownames(fpkm_table)]
  fpkm_table_protc <- subset(fpkm_table, dict[rownames(fpkm_table)] != "")
  rownames(fpkm_table_protc) <- dict[rownames(fpkm_table_protc)]
  return(fpkm_table_protc)
}

#' Filter FPKM Table
#' @details Filters out all rows where transcripts are detected in every sample
#' @param fpkm_table Table with FPKM vaules of (multiple) samples
#' @param quiet print how many rows are filtered, default FALSE
#' @export

filterDetectedTranscripts <- function(fpkm_table, quiet = FALSE){
  filter <- apply(fpkm_table,1, function(x) all(x > 0))
  fpkm_table_filtered <- subset(fpkm_table, filter)
  if(!quiet) message(paste0("Filtered "),length(which(!filter)), " rows out of ", nrow(fpkm_table))
  return(fpkm_table_filtered)
}






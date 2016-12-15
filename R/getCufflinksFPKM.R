#' Takes Cufflinks gtf output files as input and outputs a named numeric vector (names = Transcript/gene id, value = FPKM)
#' @param gtf_file_path Path to GTF file that is output by cufflinks
#' @import rtracklayer
#' @export



getCufflinksgtfFPKM<-function(gtf_file_path){

  gtf<-readGFF(gtf_file_path)
  FPKM<-as.numeric(gtf$FPKM[gtf$type=="transcript"])
  tr_names<-gsub("\\..*","",gtf$transcript_id[gtf$type=="transcript"])
  # dict<-ids$pro_name
  # names(dict)<-ids$tx_name
  # pro_names<-dict[tr_names]
  names(FPKM)<-tr_names
  FPKM<-FPKM[!(names(FPKM)=="" | is.na(names(FPKM)))]
  return(FPKM)
}



#' Takes Cufflinks genes.fpkm_tracking output files as input and outputs a named numeric vector (names = Transcript/gene id, value = FPKM)
#' @param genes.fpkm_tracking_file_path Path to GTF file that is output by cufflinks
#' @export

getCufflinksgeneFPKM<-function(genes.fpkm_tracking_file_path){

  quant_table<-read.delim(genes.fpkm_tracking_file_path)
  FPKM<-as.numeric(quant_table$FPKM)
  tr_names<-gsub("\\..*","",quant_table$gene_id)
  # dict<-ids$pro_name
  # names(dict)<-ids$tx_name
  # pro_names<-dict[tr_names]
  names(FPKM)<-tr_names
  FPKM<-FPKM[!(names(FPKM)=="" | is.na(names(FPKM)))]
  return(FPKM)
}


#' Takes Cufflinks isoforms.fpkm_tracking output files as input and outputs a named numeric vector (names = Transcript/gene id, value = FPKM)
#' @param genes.fpkm_tracking_file_path Path to GTF file that is output by cufflinks
#' @export

getCufflinkstranscriptFPKM<-function(isoforms.fpkm_tracking_file_path){

  quant_table<-read.delim(isoforms.fpkm_tracking_file_path)
  FPKM<-as.numeric(quant_table$FPKM)
  tr_names<-gsub("\\..*","",quant_table$tracking_id)
  # dict<-ids$pro_name
  # names(dict)<-ids$tx_name
  # pro_names<-dict[tr_names]
  names(FPKM)<-tr_names
  FPKM<-FPKM[!(names(FPKM)=="" | is.na(names(FPKM)))]
  return(FPKM)
}

#' Import FPKM from Cufflinks
#' @description  Takes multiple Cufflinks output files as input and outputs a data.frame with FPKM values
#' @param files Path to Cufflinks files to import. Can be isoforms.fpkm_tracking, genes.fpkm_tracking, or GTF
#' @param filetype possible values : "isoforms.fpkm_tracking", "genes.fpkm_tracking", "GTF"
#' @export


getCufflinksFPKM <- function(files, filetype = c("isoforms.fpkm_tracking", "genes.fpkm_tracking", "GTF")){

  fpkm <- list()
  for(file in files){
    sample_name <- paste0(strsplit(file, "_")[[1]][1:2],collapse = "_")
    if(filetype == "isoforms.fpkm_tracking"){
      fpkm[[sample_name]] <- getCufflinkstranscriptFPKM(isoforms.fpkm_tracking_file_path = file)
    } else if(filetype == "genes.fpkm_tracking"){
      fpkm[[sample_name]] <- getCufflinksgeneFPKM(genes.fpkm_tracking_file_path = file)
    }else if(filetype == "GTF"){
      fpkm[[sample_name]] <- getCufflinksgtfFPKM(gtf_file_path = file)
    } else{
      stop(paste("Please specify valid filetype"))
    }
    fpkm[[sample_name]] <- hela_fpkm[[sample_name]][sort(names(hela_fpkm[[sample_name]]))]
  }
  fpkm_table <- as.data.frame(do.call(cbind, fpkm))
  fpkm_table <- fpkm_table[-which(duplicated(rownames(fpkm_table))),]
  return(fpkm_table)
}


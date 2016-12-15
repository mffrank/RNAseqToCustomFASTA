#' Takes Cufflinks gtf output files as input and outputs a named numeric vector (names = Transcript/gene id, value = FPKM)
#' @param gtf_file_path Path to GTF file that is output by cufflinks
#' @import rtracklayer
#' @export



getCufflinksFPKM<-function(gtf_file_path){

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

getCufflinksgeneFPKM<-function(isoforms.fpkm_tracking_file_path){

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

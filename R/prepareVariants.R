#' Read Vcf files to GRanges
#' @details Takes VCF output files from variant calling pipeline and combines
#'  a GRanges object ready for fasta generation.
#' @param vcf_dir Path to directory with Vcf files of all replicate samples
#' @param pattern Pattern to grep the specific Vcf files in the directory
#' @param share_num Number of samples that must contain the variant.
#' Two options, percentage format or sample number.
#' @import customProDB
#' @import GenomeInfoDb
#' @export


getVcfs <- function(vcf_dir, pattern =".vcf", share_num = "100%"){

  vcfFile_path <- paste(vcf_dir,list.files(vcf_dir, pattern = pattern),sep="/")
  vcfs <- lapply(vcfFile_path, function(x) InputVcf(x))

  for (i in 1:length(vcfs)) {
    seqlevelsStyle(vcfs[[i]][[1]])<-"NCBI"
  }

  shared <- Multiple_VCF(vcfs, share_num=share_num)

  # shared_out <- getSNVandINDEL(shared, exon_annotation, procodingseq)

  return(shared)
}

#' Annotate Variants
#' @details Annotates the variants in a granges object produced by getVcfs
#' @param granges_variants GRanges object produced by getVcfs
#' @param exon_annotation Exon annotation from annotation folder
#' @param procodingseq Coding sequences from annotation folder
#' @import customProDB
#' @import rtracklayer
#' @export

getSNVandINDEL <- function(granges_variants, exon_annotation, procodingseq){

  variants_df <- values(granges_variants)
  indel_index <- sapply(variants_df$REF,nchar) != sapply(variants_df$ALT,nchar)

  indelvcf <- granges_variants[indel_index]
  snvvcf <- granges_variants[!indel_index]

  postable_snv <- Positionincoding(snvvcf, exon_annotation)
  postable_indel <- Positionincoding(indelvcf, exon_annotation)

  txlist <- unique(postable_snv[, 'txid'])
  codingseq <- procodingseq[procodingseq[, 'tx_id'] %in% txlist,]
  snvtab <- aaVariation (postable_snv, codingseq)

  txlist_indel <- unique(postable_indel[, 'txid'])
  codingseq_indel <- procodingseq[procodingseq[, 'tx_id'] %in% txlist_indel, ]

  return(list(postable_snv = postable_snv, SNV_tab = snvtab,
              postable_indel = postable_indel, codingseq_indel = codingseq_indel,
              snvvcf = snvvcf, indelvcf = indelvcf))
}

#' Generate a Variant Summary Report
#' @details Prints a table of the location of the variants. Caution:
#' Because the Varlocation function of CustomProDB is very slow, this takes a long time
#' to compute.
#' @param snv_and_indel List object produced by getSNVandINDEL
#' @param ids Id mapping from annotation folder
#' @param txdb Txdb object from annotation folder
#' @param print boolean, should the result be printed in working directory as pdf
#' or be returned as a list.
#' @import customProDB
#' @import ggplot2
#' @export


VariantSummaryReport <- function(snv_and_indel, ids, txdb, print = TRUE){

  SNVloc <- Varlocation(snv_and_indel$snvvcf,txdb,ids)
  indelloc <- Varlocation(snv_and_indel$indelvcf,txdb,ids)
  output <- list(SNV = table(SNVloc[,'location']),
                 INDEL = table(indelloc[,'location']))

  if(print){
    pdf("Variant_summary_report.pdf", width = 6, height = 4)
    # SNV's
    plot_df <- data.table(Position = names(output[[1]][-3]), SNVs = as.numeric(output[[1]][-3]), AA_Impact = FALSE)
    coding_df <- data.table(Position = "Coding", 
                            SNVs =c(nrow(unique(snv_and_indel$SNV_tab[snv_and_indel$SNV_tab$vartype == "synonymous",c("pos", "chr")])),
                                    nrow(unique(snv_and_indel$SNV_tab[snv_and_indel$SNV_tab$vartype == "non-synonymous",c("pos", "chr")]))),
                            AA_Impact = c(FALSE, TRUE))
    plot_df <- rbind(plot_df, coding_df)
    p <- ggplot(plot_df) + 
      geom_bar(aes(x=Position, y=SNVs, fill = AA_Impact), stat = "identity") + 
      ggtitle(paste0("Location of ", names(output)[1],"s")) +
      theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
      ylab(paste0("# of ", names(output[1])))
    print(p)  
    # Indels
    plot_df2 <- data.table(Position = names(output[[2]]), SNVs = as.numeric(output[[2]]),
                          AA_Impact = c(F,F,T,F,F,F,F))
    p2 <- ggplot(plot_df2) + 
      geom_bar(aes(x=Position, y=SNVs, fill = AA_Impact), stat = "identity") + 
      ggtitle(paste0("Location of ", names(output)[2],"s")) +
      theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
      ylab(paste0("# of ", names(output[2])))
    
    # pl <- rbind(plot_df[, Variation_type := "SNV"], plot_df2[, Variation_type := "INDEL"])
    # p <- ggplot(pl) + 
    #   geom_bar(aes(x=Position, y=SNVs, fill = AA_Impact), stat = "identity") + 
    #   ggtitle(paste0("Location of ", names(output)[2],"s")) +
    #   theme(axis.text.x = element_text(angle=30, hjust=1, vjust=1)) +
    #   ylab(paste0("# of ", names(output[2]))) +
    #   facet_grid(. ~ Variation_type)
    print(p2)  
    dev.off()
  } else{
    return(output)
  }

  table(SNVloc[,'location'])
  table(indelloc[,'location'])

}

#' Output Splice Junction FASTA
#' @description  Takes STAR output files (SJ.out.tab as input and outputs a
#' named numeric vector (names = Transcript/gene id, value = FPKM)
#' @param jun Junction file produced by the JunctionType function
#' @param exons Exon table from annotation folder
#' @param procodingseq Protein coding sequences from annotation folder
#' @import GenomicRanges


#' @export


OutputJunctions <- function(outfile, juncs, exons, procodingseq, proteinseq, genome, alt_spl = TRUE, all_junc_peptides = TRUE, trimTry = TRUE){

  if(alt_spl){
    alt_spl_out <- getAltSpliceProteins(juncs, exon, procodingseq)

    juncs_noalt <- subset(juncs, jun_type != "novel alternative splicing junction")
  }else{
    alt_spl_out <- numeric(0)
  }
  if(all_junc_peptides){
    junpeps_out <-OutputNovelJun(junction_type = juncs_noalt, genome = genome, trimTry = trimTry, proteinseq = proteinseq)
  }else{
    junpeps_out <- numeric(0)
  }
  all_out <- rbind(alt_spl_out,junpeps_out)

  tmp <- paste('>', all_out[, 1], '\n', all_out[, 2], sep='')
  write(tmp, file=outfile)


}

#' Get alternatively spliced proteins
#' @description  Uses junction table to produce spliced sequences of known transcripts. A protein is produced when a
#' junction connects two of its exons that are not adjacent.
#' @param jun Junction file produced by the JunctionType function
#' @param exons Exon table from annotation folder
#' @param procodingseq Protein coding sequences from annotation folder
#' @import GenomicRanges
#' @import Biostrings

#' @export


getAltSpliceProteins <- function(jun, exon, procodingseq){

  # tr_part1 <- strsplit(juncs$tx_name_part1, split = ",")
  # tr_part2 <- strsplit(juncs$tx_name_part2, split = ",")
  #
  # transcripts <- lapply(1:nrow(juncs), function(i){
  #   intersect(tr_part1[[i]], tr_part2[[i]])
  # })

  exon_cds <- subset(exon, !is.na(cds_start) & !is.na(cds_end))
  cdsanno <- GRanges(seqnames = exon_cds$chromosome_name,
                     ranges = IRanges(start = exon_cds$cds_chr_start, end = exon_cds$cds_chr_end),
                     strand = exon_cds$strand, trans_id = exon_cds$tx_id,
                     trans_name = exon_cds$tx_name, pro_name = exon_cds$pro_name,
                     gene_name = exon_cds$gene_name, cds_s2coding = exon_cds$cds_start,
                     cds_e2coding = exon_cds$cds_end, ex_no = exon_cds$rank)
  jun <- subset(jun, jun_type == "novel alternative splicing junction")
  junRange1 <- GRanges(seqnames = jun$seqnames, ranges = IRanges(start = jun$part1_sta,
                                                                 end = jun$part1_end),
                       strand = jun$strand, junction_id = jun$id)
  junRange2 <- GRanges(seqnames = jun$seqnames, ranges = IRanges(start = jun$part2_sta,
                                                                 end = jun$part2_end),
                       strand = jun$strand, junction_id = jun$id, type = jun$jun_type,
                       cov = jun$cov)

  match_part1 <- findOverlaps(junRange1, cdsanno, type = "end")
  match_part2 <- findOverlaps(junRange2, cdsanno, type = "start")

  tr_part1 <- data.frame(Jun_id = queryHits(match_part1), EX_id_1 = subjectHits(match_part1),
                         tx_name = cdsanno[subjectHits(match_part1)]$trans_name,
                         ex_no_1 = cdsanno[subjectHits(match_part1)]$ex_no,
                         cds_e1coding = cdsanno[subjectHits(match_part1)]$cds_e2coding,
                         cds_s1coding = cdsanno[subjectHits(match_part1)]$cds_s2coding)

  tr_part2 <- data.frame(Jun_id = queryHits(match_part2), EX_id_2 = subjectHits(match_part2),
                         tx_name = cdsanno[subjectHits(match_part2)]$trans_name,
                         ex_no_2 = cdsanno[subjectHits(match_part2)]$ex_no,
                         cds_e2coding = cdsanno[subjectHits(match_part2)]$cds_e2coding,
                         cds_s2coding = cdsanno[subjectHits(match_part2)]$cds_s2coding,
                         strand = strand(cdsanno[subjectHits(match_part2)]),
                         type =junRange2[queryHits(match_part2)]$type,
                         cov = junRange2[queryHits(match_part2)]$cov)

  junc_tr <- merge(tr_part1,tr_part2)
  # Alternatively, if junctions are not filtered for alternative splicing
  # canoncal_idx <- subset(junc_tr, abs(junc_tr$ex_no_1 - junc_tr$ex_no_2) == 1)$Jun_id
  # junc_tr <- subset(junc_tr, !(Jun_id %in% canoncal_idx))

  junc_seq <- merge(junc_tr, procodingseq)
  junc_seq$spliced_seq <- sapply(1:nrow(junc_seq), function(i){
    if(junc_seq$strand[i] == "-"){
      splice_mar <- junc_seq[i,c("cds_s1coding","cds_e2coding")]
    }else{
      splice_mar <- junc_seq[i,c("cds_s2coding","cds_e1coding")]
    }
    paste(substr(junc_seq$coding[i],1,min(splice_mar)),
          substr(junc_seq$coding[i],max(splice_mar),nchar(junc_seq$coding[i])),
          sep = "")
  })

  junc_seq$spliced_aa <- Biostrings::translate(DNAStringSet(junc_seq$spliced_seq),if.fuzzy.codon = "solve")
  # Fuzzy translations are skipped
  junc_seq$spliced_aa <- gsub("X", "", junc_seq$spliced_aa)
  # Proteins with stop codons are cleaved
  junc_seq$spliced_aa <- gsub("\\*.*", "", junc_seq$spliced_aa)
  # Junctions that generate the same coding sequence as another junction are removed
  junc_seq <- junc_seq[!duplicated(junc_seq$spliced_aa),]
  # Resulting sequences with less than 7 AA are deleted
  junc_seq <- subset(junc_seq, nchar(junc_seq$spliced_aa) >=7)
  # Add gene name
  dict <- unique(exon_cds[,c("tx_name", "gene_name")])
  junc_seq <- merge(junc_seq, dict, all.x = T)

  header <- sapply(1:nrow(junc_seq), function(i){
    if(junc_seq$strand[i] == "-"){
      splice_mar <- junc_seq[i,c("cds_s1coding","cds_e2coding")]
    }else{
      splice_mar <- junc_seq[i,c("cds_s2coding","cds_e1coding")]
    }
    paste0(junc_seq$pro_name[i],"_", min(splice_mar), ":", max(splice_mar), " ALTSPLICEJUNC", i, "|",
                     junc_seq$tx_name[i], "|", junc_seq$gene_name[i], "|detected in ", junc_seq$cov[i], " samples")
  })

  output <- cbind(header, junc_seq$spliced_aa)
  return(output)
}


### Custom Function for trimming splice variants to second Lys/Arg after junction

##' Three-frame translation of novel junctions. And remove those could be found in normal protein sequences.
##' This function requires a genome built by BSgenome package.
##'
##' @title Generate peptide FASTA file that contains novel junctions.
##' @param junction_type a data frame which is the output of function JunctionType()
##' @param genome a BSgenome object. (e.g. Hsapiens)
##' @param outfile output file name
##' @param proteinseq a data frame cotaining amino acid sequence for each protein.
##' @param ... Additional arguments.
##' @return FASTA file that contains novel junction peptides.
##' @author Xiaojing Wang
##' @author Max Frank

OutputNovelJun <- function(junction_type, genome,
                           trimTry = F,
                           proteinseq, ...)
{
  options(stringsAsFactors=FALSE)
  #ids <- subset(ids,pro_name!='')

  #trans <- transcripts(txdb)
  #index <- which(values(trans)[['tx_name']] %in% ids[,'tx_name'])
  #pro_trans <- trans[index]
  novel_junc <- subset(junction_type, jun_type != 'known junction')
  if(!length(grep('chr', novel_junc[, 'seqnames'], fixed=T))>0) {
    novel_junc[, 'seqnames'] <- paste('chr', novel_junc[, 'seqnames'], sep='')
    idx <- which(novel_junc[, 'seqnames'] %in% seqnames(genome))
    novel_junc <- novel_junc[idx, ]
  }

  ###remove abnormal junctions
  idx_abn <- union(which(novel_junc[, 'start'] < 0),
                   which(novel_junc[, 'end'] < 0))
  if(length(idx_abn > 0)) novel_junc <- novel_junc[-idx_abn, ]

  junRange1 <- GRanges(seqnames=novel_junc$seqnames,
                       ranges=IRanges(start=novel_junc$part1_sta,
                                      end=novel_junc$part1_end),
                       strand=novel_junc$strand,
                       junction_id=novel_junc$id)

  junRange2 <- GRanges(seqnames=novel_junc$seqnames,
                       ranges=IRanges(start=novel_junc$part2_sta,
                                      end=novel_junc$part2_end),
                       strand=novel_junc$strand,
                       junction_id=novel_junc$id)

  #match1_protx <- findOverlaps(junRange1,pro_trans)
  #match2_protx <- findOverlaps(junRange2,pro_trans)

  #juntransRange1 <- junRange1[unique(queryHits(match1_protx))]
  #juntransRange2 <- junRange2[unique(queryHits(match2_protx))]

  #junseq1 <- getSeq(genome,'chr1',start=1000,end=2000,as.character=TRUE)
  ###already did reverseComplement
  junseq1 <- getSeq(genome, junRange1)
  junseq2 <- getSeq(genome, junRange2)

  junseq_cat <- DNAStringSet(mapply(function(x, y, z)
    ifelse(z == '+', paste(x, y, sep=''), paste(y, x, sep='')),
    as.data.frame(junseq1)[, 1],
    as.data.frame(junseq2)[, 1], as.character(strand(junRange1))))

  #index_plus <- which(strand(junRange1) == '+')
  #index_minus <- which(strand(junRange1) == '-')
  #seqs_plus <- junseq_cat[index_plus]
  #seqs_minus <- reverseComplement(junseq_cat[index_minus])
  #seqs <- c(seqs_plus, seqs_minus)

  #novel_junc_new <- rbind(novel_junc[index_plus, ],
  #                                novel_junc[index_minus, ])
  novel_junc_new <- novel_junc
  seqs <- junseq_cat

  ##Remove sequences contains NNN
  Nindx <- grep('N', seqs)
  if(length(Nindx) > 0){
    seqs <- seqs[-Nindx]
    novel_junc_new <- novel_junc_new[-Nindx, ]
  }

  seqs_name <- paste(paste(novel_junc_new[, 'id'], '_',
                           novel_junc_new[, 'seqnames'], ':',
                           novel_junc_new[, 'start'], '-', novel_junc_new[, 'end'],
                           sep=''), novel_junc_new[, 'cov'],sep='|')

  junpepcoding <- data.frame('pro_name'=seqs_name,
                             'coding'=as.data.frame(seqs)[, 1])
  #save(junpepcoding, file=outfile_c)

  peptides_r1 <- Biostrings::translate(seqs)
  peptides_r2 <- Biostrings::translate(subseq(seqs, start=2))
  peptides_r3 <- Biostrings::translate(subseq(seqs, start=3))

  junpos_rna_p1 <- ifelse(novel_junc_new[, 'strand'] == '+',
                          as.numeric(novel_junc_new[, 'part1_len']),
                          as.numeric(novel_junc_new[, 'part2_len']))
  junpos_rna_p2 <- ifelse(novel_junc_new[, 'strand'] == '+',
                          as.numeric(novel_junc_new[, 'part1_len'])+1,
                          as.numeric(novel_junc_new[, 'part2_len'])+1)

  junpos_r1_p1 <- ceiling(junpos_rna_p1/3)
  junpos_r1_p2 <- ceiling(junpos_rna_p2/3)

  junpos_r2_p1 <- ceiling((junpos_rna_p1-1)/3)
  junpos_r2_p2 <- ceiling((junpos_rna_p2-1)/3)

  junpos_r3_p1 <- ceiling((junpos_rna_p1-2)/3)
  junpos_r3_p2 <- ceiling((junpos_rna_p2-2)/3)

  name_r1 <- paste(paste(novel_junc_new[, 'id'], '_', novel_junc_new[, 'seqnames'],
                         ':', novel_junc_new[, 'start'], '-', novel_junc_new[, 'end'],
                         sep=''), novel_junc_new[, 'cov'], 'ORF1 ',
                   paste('Junpos:', junpos_r1_p1, '-', junpos_r1_p2, sep=''),
                   novel_junc_new[, 'strand'], novel_junc_new[, 'tx_name_part1'],
                   novel_junc_new[, 'tx_name_part2'],
                   novel_junc_new[, 'jun_type'], sep='|')
  name_r2 <- paste(paste(novel_junc_new[, 'id'],'_', novel_junc_new[, 'seqnames'],
                         ':',novel_junc_new[, 'start'], '-', novel_junc_new[, 'end'],
                         sep=''), novel_junc_new[, 'cov'], 'ORF2 ',
                   paste('Junpos:', junpos_r2_p1, '-', junpos_r2_p2, sep=''),
                   novel_junc_new[, 'strand'], novel_junc_new[, 'tx_name_part1'],
                   novel_junc_new[, 'tx_name_part2'],
                   novel_junc_new[, 'jun_type'], sep='|')
  name_r3 <- paste(paste(novel_junc_new[, 'id'],'_', novel_junc_new[, 'seqnames'],
                         ':',novel_junc_new[, 'start'], '-', novel_junc_new[, 'end'],
                         sep=''), novel_junc_new[, 'cov'], 'ORF3 ',
                   paste('Junpos:', junpos_r3_p1, '-', junpos_r3_p2, sep=''),
                   novel_junc_new[, 'strand'], novel_junc_new[, 'tx_name_part1'],
                   novel_junc_new[, 'tx_name_part2'],
                   novel_junc_new[, 'jun_type'], sep='|')




  all_pep<- rbind(cbind(name_r1, as.data.frame(peptides_r1)[, 1]),
                  cbind(name_r2, as.data.frame(peptides_r2)[, 1]),
                  cbind(name_r3, as.data.frame(peptides_r3)[, 1]))

  ### remove peptide contain stop codon
  index_stop <- grep('*', all_pep[, 2], fixed=T)
  if(length(index_stop) > 0){
    all_pep_rmstop <- all_pep[-index_stop, ]
  }else all_pep_rmstop <- all_pep
  ### check if any peptides can be found in the normal database, remove those

  index_nor <- lapply(all_pep_rmstop[, 2], function(x)
    grep(x, proteinseq[, 'peptide'], fixed=T))
  index_nor <- which(unlist(lapply(index_nor, length)) > 0)

  if(length(index_nor) > 0){
    all_pep_new <- all_pep_rmstop[-index_nor, ]
  }else all_pep_new <- all_pep_rmstop

  ### Trim bases to second Lys/Arg after/before junction (for star input)
  if(trimTry){

    trim_right <- apply(all_pep_new, 1,function(x){
      .instr( str1 = x[2], str2 = "[KR]",
              startpos = as.numeric(sub("\\|.*","",strsplit(x[1], "Junpos:[0-9]*-")[[1]][2])),n = 2)
    })

    # for (i in 1:length(trim_right)){
    #   if (trim_right[i] == 0){
    #     trim_right[i] <- nchar(all_pep_new[i,2])
    #   }
    # }

    trim_left <- apply(all_pep_new, 1,function(x){
      nchar(x[2]) -
        .instr( str1 = paste(rev(strsplit(x[2],"")[[1]]),collapse=""),
                str2 = "[KR]",
                startpos = nchar(x[2]) - (as.numeric(sub("\\|.*","",strsplit(x[1], "Junpos:[0-9]*-")[[1]][2])) - 2),n = 2)
    })
    # for (i in 1:length(trim_right)){
    #   if (trim_left[i] == nchar(all_pep_new[i,2])){
    #     trim_left[i] <- 0
    #   }
    # }

    trimmed_pep <- unlist(lapply(1:nrow(all_pep_new), function(i){
      substr(all_pep_new[i,2],trim_left[i]+1,trim_right[i])
    }))

    all_pep_new[,2] <- trimmed_pep
    all_pep_new <- all_pep_new[-which(nchar(all_pep_new[,2]) < 7),,drop = F]
  }
  return(all_pep_new)

}


## Helper function to detect nth ocuurence of query string
.instr <- function(str1,str2,startpos=1,n=2){
  aa=unlist(strsplit(substring(str1,startpos),str2))
  if(length(aa) < n+1 ){
    return(0)
  }else{
    return(sum(nchar(aa[1:n])) + startpos+(n-1) )
  }
}

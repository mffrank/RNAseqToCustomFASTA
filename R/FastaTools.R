
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
#' @param addcRAP logical, Should cRAP and affinity tag sequences be appended
#' @param addiRT logical, Should iRT peptides be appended
#' @import seqinr
#' @import data.table
#' @export

concatenateFastas <- function(path, pattern = "\\.fa", outfile = "combined.fasta", addcRAP = TRUE, addiRT = TRUE, saveMappingTable = TRUE){
  files <- list.files(path = path, pattern = pattern)
  files <- files[! files == outfile]

  fasta_list <- lapply(files, function(x) readLines(file(paste0(path,"/",x))))
  closeAllConnections()
  fasta_list <- lapply(fasta_list, function(x){
    data.table(header = x[seq(1,length(x),by = 2)],
               sequence = x[seq(2,length(x),by = 2)],
               stringsAsFactors = F)
  })
  names(fasta_list) <- gsub("\\..*","", files)

  # Remove stars at end positions
  for(i in 1:length(fasta_list)){
    star_end_pos <- which(substr(fasta_list[[i]]$sequence,nchar(fasta_list[[i]]$sequence),
                                 nchar(fasta_list[[i]]$sequence)) == "*")
    fasta_list[[i]][star_end_pos,"sequence"] <- gsub("\\*","",fasta_list[[i]]$sequence[star_end_pos])
  }

  # message("Mapping headers to gene level...")
  # fasta_list <- .mapToGeneLevel(fasta_list)
  # message("Merging duplicate sequences...")
  # # Merge canonical duplicates
  # fasta_list <- .mergeIdenticalSeq(fasta_list)

  message("Combining ", paste("\n",files), "\n into ", outfile)
  fasta_df <- do.call("rbind", fasta_list)

  message("Mapping headers to gene level...")
  fasta_df$gene <- .getGene(fasta_df)
  fasta_df$Isoform <- .getIsoform(fasta_df)
  fasta_df <- .changeHeaders(fasta_df)

  n_positions <- grep("[n]", fasta_df$sequence)
  # fasta_df[n_positions,"sequence"] <- gsub("n","X",fasta_df[n_positions,"sequence"])
  fasta_df <- fasta_df[-n_positions,]
  message(paste0("\nWarning: removed ", fasta_df[n_positions,header], " from db. Ambigous(n) AAs \n\n"),
          paste0("Removed ", length(n_positions), " Sequences in total"))

  if(saveMappingTable){
    fasta_df$protein <- .getProtein(fasta_df)
    saveRDS(fasta_df, file = paste0(path,gsub("\\..*","",outfile),"_header_mapping_v2.rda"))
    fasta_df[, protein := NULL]
  }

  fasta_df[, c("gene", "Isoform", "header_old") := NULL]
  # @TODO get this into data loading
  iRT <- readRDS("Y:/Master_Project/src/iRT.rda")
  cRap_tags <- readRDS("Y:/Master_Project/src/cRAP_tags.rda")
  if(addiRT){
    message("Adding iRT peptides...")
    fasta_df <- rbind(fasta_df, iRT)
  }
  if(addcRAP){
    message("Adding cRAP peptides...")
    fasta_df <- rbind(fasta_df, cRAP_tags)
  }
  message("Writing output...")
  fasta_df$header <- substr(fasta_df$header,2,nchar(fasta_df$header))
  write.fasta(sequences = as.list(fasta_df$sequence), names = fasta_df$header,
              as.string = TRUE, file.out = paste0(path,outfile))
  # write(allfasta, file = paste0(path,outfile))
}


.mapToGeneLevel <- function(fasta_list){
  # Canonical Proteins
  can_idx <- which(sapply(fasta_list, function(x) grepl("ENSP[0-9]* ",x[1,1])))
  fasta_list[[can_idx]]$gene <- sapply(strsplit(fasta_list[[can_idx]]$header, split = "\\|"), "[",4)
  # fasta_list[[can_idx]] <- .changeHeaders(fasta = fasta_list[[can_idx]], can = TRUE)
  # SNVs
  can_idx <- which(sapply(fasta_list, function(x) grepl("ENSP[0-9]*_[[:upper:]]",x[1,1])))
  fasta_list[[can_idx]]$gene <- sapply(strsplit(fasta_list[[can_idx]]$header, split = "\\|"), "[",3)
  # fasta_list[[can_idx]] <- .changeHeaders(fasta = fasta_list[[can_idx]], can = FALSE)

  # INDELs
  can_idx <- which(sapply(fasta_list, function(x) grepl("ENSP[0-9]*_[0-9]*\\:[[:upper:]]*>",x[1,1])))
  fasta_list[[can_idx]]$gene <- sapply(strsplit(fasta_list[[can_idx]]$header, split = "\\|"), "[",3)
  # fasta_list[[can_idx]] <- .changeHeaders(fasta = fasta_list[[can_idx]], can = FALSE)

  # ALTSPLICEJUNCS
  can_idx <- which(sapply(fasta_list, function(x) grepl("ALTSPLICEJUNC",x[1,1])))
  altsplidx <- grep("ALTSPLICEJUNC", fasta_list[[can_idx]]$header)
  fasta_list[[can_idx]]$gene <- sapply(strsplit(fasta_list[[can_idx]]$header, split = "\\|"), "[",3)
  # fasta_list[[can_idx]][altsplidx] <- .changeHeaders(fasta = fasta_list[[can_idx]][altsplidx], can = FALSE)


  return(fasta_list)
  }

.getGene <- function(fasta){
  canidx <- grep("ENSP[0-9]* ", fasta$header)
  juncidx <- grep(">JUNC", fasta$header)
  rest <- setdiff(1:nrow(fasta), c(canidx, juncidx))
  can_gene <- sapply(strsplit(fasta[canidx]$header, split = "\\|"), "[",4)
  nocan_gene <- sapply(strsplit(fasta[rest]$header, split = "\\|"), "[",3)
  out <- vector()
  out[canidx] <- can_gene
  out[juncidx] <- NaN
  out[rest] <- nocan_gene
  return(out)
}

.getIsoform <- function(fasta, can = F){

  setkey(fasta, gene, sequence)
  # not pretty
  iso <- 1
  is <- 1
  for(i in 2:nrow(fasta)){
    if(fasta[i,gene] == fasta[i-1,gene]){
      if(fasta[i,sequence] != fasta[i-1,sequence]){
        is <- is + 1
      }
    }else{
      is <- 1
    }
    iso <- c(iso,is)
  }
 return(iso)
}

.changeHeaders <- function(fasta){
  fasta$header_old <- fasta$header
  canidx <- grep("ENSP[0-9]* ", fasta$header)
  juncidx <- grep(">JUNC", fasta$header)
  rest <- setdiff(1:nrow(fasta), c(canidx, juncidx))
  # @TODO fix greps so they have 3 entries separated by | (swissprot norm)
  fasta[canidx]$header <-paste0(">cf|",
                          fasta[canidx]$gene, "-",fasta[canidx]$Isoform,
                          gsub("(>)(.*?(?= ))( \\|)(.*?(?=\\|))(\\|)(.*?(?=\\|))(.*)","|\\2-\\4 \\6\\7",
                               fasta[canidx]$header, perl = T))
  fasta[rest]$header <- paste0(">cf|",
                           fasta[rest]$gene, "-",fasta[rest]$Isoform,
                           gsub("(>)(.*?(?=_))(.*?(?= ))( \\|| )(.*)","\\3|\\2 \\5",
                                fasta[rest]$header, perl = T))

  # fasta[, c("gene","Isoform") := NULL]
  return(fasta)
}


saveMappingTable <- function(path, pattern = "\\.fa", outfile = "header_mapping_table.rda"){
  files <- list.files(path = path, pattern = pattern)
  files <- files[! files == outfile]

  fasta_list <- lapply(files, function(x) readLines(file(paste0(path,"/",x))))
  closeAllConnections()
  fasta_list <- lapply(fasta_list, function(x){
    data.table(header = x[seq(1,length(x),by = 2)],
               sequence = x[seq(2,length(x),by = 2)],
               stringsAsFactors = F)
  })
  names(fasta_list) <- gsub("\\..*","", files)

  # Remove stars at end positions
  for(i in 1:length(fasta_list)){
    star_end_pos <- which(substr(fasta_list[[i]]$sequence,nchar(fasta_list[[i]]$sequence),
                                 nchar(fasta_list[[i]]$sequence)) == "*")
    fasta_list[[i]][star_end_pos,"sequence"] <- gsub("\\*","",fasta_list[[i]]$sequence[star_end_pos])
    removeStars <- function(sequence){
      
    
    star_end_pos <- which(substr(sequence,nchar(sequence),
                                 nchar(sequence)) == "*")
    sequence[star_end_pos] <- gsub("\\*","",sequence[star_end_pos])
    return(sequence)
    }
  }

  # message("Mapping headers to gene level...")
  # fasta_list <- .mapToGeneLevel(fasta_list)
  # message("Merging duplicate sequences...")
  # # Merge canonical duplicates
  # fasta_list <- .mergeIdenticalSeq(fasta_list)

  message("Combining ", paste("\n",files), "\n into ", outfile)
  fasta_df <- do.call("rbind", fasta_list)

  message("Mapping headers to gene level...")
  fasta_df$gene <- .getGene(fasta_df)
  fasta_df$Isoform <- .getIsoform(fasta_df)
}

.getProtein <- function(fasta){
  if("header_old" %in% names(fasta)){
    prot <- gsub(">| .*", "", fasta$header_old)
  }else{
    prot <- gsub(">| .*", "", fasta$header)
  }
  return(prot)
}

.mergeIdenticalSeq <- function(fasta_list){
  can_idx <- which(sapply(fasta_list, function(x) grepl("ENSP[0-9]* ",x[1,1])))
  dupseq <- unique(fasta_list[[can_idx]][duplicated(fasta_list[[can_idx]]$sequence),sequence])
  dup <- lapply(as.character(dupseq), function(x) fasta_list[[can_idx]][fasta_list[[can_idx]]$sequence == x,header])
  fasta_list[[can_idx]] <- fasta_list[[can_idx]][!duplicated(fasta_list[[can_idx]]$sequence),]
  dup <- lapply(dup, strsplit, split = "\\||>| \\|")
  new_head <- sapply(1:length(dup), function(i){
    meanexpr <- mean(as.numeric(sapply(dup[[i]], "[",3)))
    names <- sapply(dup[[i]], "[",2)
    paste0(">", paste0(names, collapse = "/"), " |", meanexpr, "|", paste0(dup[[i]][[1]][4:7], collapse = "|"))
  })
  fasta_list[[can_idx]][fasta_list[[can_idx]]$sequence %in% dupseq, "header"] <- new_head
  # Merge SNV duplicates
  can_idx <- which(sapply(fasta_list, function(x) grepl("ENSP[0-9]*_[[:upper:]]",x[1,1])))
  dupseq <- unique(fasta_list[[can_idx]][duplicated(fasta_list[[can_idx]]$sequence),sequence])
  dup <- lapply(dupseq, function(x) fasta_list[[can_idx]][fasta_list[[can_idx]]$sequence == x,header])
  fasta_list[[can_idx]] <- fasta_list[[can_idx]][!duplicated(fasta_list[[can_idx]]$sequence),]
  dup <- lapply(dup, strsplit, split = "\\||>| \\|")
  new_head <- sapply(1:length(dup), function(i){
    names <- sapply(dup[[i]], function(x) gsub("_.*","",x[2]))
    aa_change <- gsub("^.*_","",dup[[i]][[1]][2])
    paste0(">", paste0(names, collapse = "/"), "_", aa_change, " |", paste0(dup[[i]][[1]][3:6], collapse = "|"))
  })
  fasta_list[[can_idx]][fasta_list[[can_idx]]$sequence %in% dupseq, "header"] <- new_head
  # Merge INDEL duplicates
  can_idx <- which(sapply(fasta_list, function(x) grepl("ENSP[0-9]*_[0-9]*\\:[[:upper:]]*>",x[1,1])))
  dupseq <- unique(fasta_list[[can_idx]][duplicated(fasta_list[[can_idx]]$sequence),sequence])
  dup <- lapply(dupseq, function(x) fasta_list[[can_idx]][fasta_list[[can_idx]]$sequence == x,header])
  fasta_list[[can_idx]] <- fasta_list[[can_idx]][!duplicated(fasta_list[[can_idx]]$sequence),]
  dup <- lapply(dup, strsplit, split = "\\|| \\|")
  new_head <- sapply(1:length(dup), function(i){
    names <- sapply(dup[[i]], function(x) gsub("(^>)|(_.*)","",x[1]))
    base_change <- gsub("^.*_","",dup[[i]][[1]][1])
    paste0(">", paste0(names, collapse = "/"), "_", base_change, " |", paste0(dup[[i]][[1]][2:5], collapse = "|"))
  })
  fasta_list[[can_idx]][fasta_list[[can_idx]]$sequence %in% dupseq, "header"] <- new_head
  return(fasta_list)
}






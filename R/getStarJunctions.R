
#' Read in STAR junctions
#' @description  Takes STAR output files (SJ.out.tab as input and outputs a
#' named numeric vector (names = Transcript/gene id, value = FPKM)
#' @param SJout_dir Path to Star output files
#' @import GenomicRanges
#' @import customProDB
#' @export

getStarJunctions <- function(SJout_dir, splicemax, txdb, ids, pattern = "\\.tab", share_num = 1, extend = 100, skip = 0, covfilter_unique = 1, covfilter_multi = 0, overhang = "max_overhang"){

  starfiles <- paste(SJout_dir, list.files(SJout_dir, pattern),sep="")
  message("Converting Star output to Range format...")
  jun <- lapply(starfiles, function(x) StarSJ2Range(x, skip, covfilter_unique, covfilter_multi, overhang))

  for (i in 1:length(jun)) {
    seqlevelsStyle(jun[[i]])<-"NCBI"
  }
  message("Collapsing replicates...")
  sharedjun <- SharedJunc(jun, share_num=share_num, ext_up=extend, ext_down=extend)

  return(sharedjun)
}




StarSJ2Range <- function(SJout_tab, skip = 0, covfilter_unique = 1, covfilter_multi = 0, overhang = "max_overhang", ...){

  jun <- read.table(SJout_tab, sep = "\t", header = F, quote = "\"",
                    stringsAsFactors = F, skip = skip)
  names(jun) <- c("chr","start","end","strand_star","motif","annotated", "cov_unique","cov_multi","max_overhang")
  conv_strand <- c("*","+","-")
  names(conv_strand) <- c("0","1","2")
  jun$strand <- conv_strand[as.character(jun$strand_star)]

  # Coverage Filter and remove Mitochondrial SJs
  jun_filter <- subset(jun, cov_unique >= covfilter_unique & cov_multi >= covfilter_multi)
  jun_filter <- subset(jun_filter, !(chr == "chrM" | chr == "MT"))

  if(overhang == "max_overhang"){
    max_overhang <- as.numeric(jun_filter$max_overhang)
  }else{
    max_overhang <- rep(overhang,nrow(jun_filter))
  }


  part1_sta <- as.numeric(jun_filter$start - max_overhang)
  part1_end <- as.numeric(jun_filter$start - 1)
  part2_sta <- as.numeric(jun_filter$end + 1)
  part2_end <- as.numeric(jun_filter$end + max_overhang)

  junRange <- GRanges(seqnames = jun_filter$chr, ranges = IRanges(start = part1_end,
                                                                  end = part2_sta), strand = jun_filter$strand,
                      id = paste("jun",1:nrow(jun_filter),sep = ""),
                      cov = jun_filter$cov_unique, part1_len = max_overhang,
                      part2_len = max_overhang, part1_sta = part1_sta,
                      part1_end = part1_end, part2_sta = part2_sta,
                      part2_end = part2_end)
  return(junRange)
}

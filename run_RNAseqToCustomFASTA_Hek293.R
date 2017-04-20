# Produce a custom Fasta file from RNAseq data

library(customProDB)
library(devtools)
# install_github("mffrank/RNAseqToCustomFASTA")
install("Y:/Master_Project/src/RNAseqToCustomFASTA")
library(RNAseqToCustomFASTA)
# load_all("Y:/Master_Project/src/RNAseqToCustomFASTA")

#--------------------------
## Execution Parameters
annotation_path <- "Y:/Master_Project/data/Human_genome/GRCh38/CustomProDB_annotation"
output_directory <- "Y:/Master_Project/results/2017-01-18/Hek293_RNAseq_to_Fasta/Custom_fastas/"
rpkm_directory <- "Y:/Master_Project/data/HEK293/RNA_seq/Assembled_Transcriptome/"
vcf_directory <-"Y:/Master_Project/data/HEK293/RNA_seq/Varaints/"
bed_directory <- "Y:/Master_Project/data/HEK293/RNA_seq/Splice_Junctions/"
sample_name <- "Hek293"
wd <- "Y:/Master_Project/results/2017-01-18/Hek293_RNAseq_to_Fasta/"

## Download all nescessary annotations from ENSEMBL - only has to be carried out once
# downloadAnnotations("Y:/Master_Project/data/Human_genome/GRCh38/CustomProDB_annotation")

# Load the annotation objects into the working env.

load(paste(annotation_path, "/exon_anno.RData", sep = ""))
load(paste(annotation_path,"/proseq.RData", sep = ""))
load(paste(annotation_path,"/ids.RData", sep = ""))
load(paste(annotation_path, "/procodingseq.RData", sep = ""))
load(paste(annotation_path,"/splicemax.RData", sep = ""))
txdb <- loadDb(paste(annotation_path,"/txdb.sqlite", sep = ""))


## Step one: Generate Fasta according to expression levels
# -------------------------------------------------------

# Load FPKM values into data frame
setwd(rpkm_directory)
files <- list.files(pattern = paste0(sample_name,".*_isoforms.*"))

fpkm_table <- getCufflinksFPKM(files = files, filetype = "isoforms.fpkm_tracking")
setwd(wd)
# saveRDS(fpkm_table, "Y:/Master_Project/results/2016-12-15/HeLa_custom_fasta/Hela_fpkm_table")

# Look for protein coding transcripts and assign Protein id

FPKMsummaryReport(fpkm_table = fpkm_table, id_table = ids, qt = c(0.2, 0.25, 0.28, 0.3))

# Remove non-proteincoding transcripts
fpkm_table_prot <- getProteinCoding(fpkm_table = fpkm_table,id_table = ids,
                                    transcript_id_name = "tx_name", protein_id_name = "pro_name")
# Filter out all trancripts with 0 FPKM in any replicate, leaving only consistently measured transcrips
fpkm_table_prot_measured <- filterDetectedTranscripts(fpkm_table = fpkm_table_prot)

# Set cutoff at 25% of all consistently measured transcripts
cutoff <- quantile(apply(fpkm_table_prot_measured,1, mean),0.28)
nrow(fpkm_table_prot_measured[])

# Build the filtered Fasta library
setwd(output_directory)
OutputsharedPro(fpkm_table_prot_measured, cutoff=cutoff, share_sample=2, 
                outfile = "expressed_canonical_proteins.fasta", ids = ids, proteinseq =proteinseq)
#39675 proteins were generated

## Step two: Generate Fasta with SNPs and INDELs
# -------------------------------------------------------

shared_variants <- getVcfs(vcf_dir = vcf_directory, pattern = ".vcf$", share_num = 2)

snv_indels <- getSNVandINDEL(granges_variants = shared_variants, exon_annotation = exon, procodingseq = procodingseq)
saveRDS(snv_indels, "snv_indels.rda")
setwd(wd)
variant_summary <- VariantSummaryReport(snv_and_indel = snv_indels, ids = ids, txdb = txdb, print = F)
saveRDS(variant_summary, "Variant_location_table.rda")

setwd(output_directory)
OutputVarproseq_single(vartable = snv_indels$SNV_tab, proteinseq = proteinseq,
                outfile = "SNV_containing_proteins.fasta",ids = ids)

setwd(output_directory)
Outputaberrant(positiontab = snv_indels$postable_indel, coding = snv_indels$codingseq_indel, proteinseq=proteinseq, 
               outfile="INDEL_containing_proteins.fasta", ids=ids)

## Step three: Generate Fasta with alternative Splicing Proteins
# -------------------------------------------------------
setwd(wd)
pdf("Remaining_Junctions_different_filters.pdf")
filterSpacePlot(covfilter_unique = c(1,2,3,4,5,10), share_num = c(1,2), bed_directory, splicemax, txdb, ids)
dev.off()

juncs <- getStarJunctions(bed_directory, splicemax = splicemax, txdb = txdb, ids = ids,
                          pattern = "\\.tab", share_num = 2, extend = 100, skip = 0,
                          covfilter_unique = 3, covfilter_multi = 0, overhang = "max_overhang")

saveRDS(juncs,"Junctions_cov3_share2.rda")
junction_type <- JunctionType(juncs, splicemax, txdb, ids)
saveRDS(junction_type,"Junction_type_cov3_share2.rda")

library("BSgenome.Hsapiens.UCSC.hg38")
setwd(output_directory)
OutputJunctions(juncs = junction_type, procodingseq = procodingseq, proteinseq = proteinseq,
                genome = Hsapiens, exons = exon, outfile="Splice_Junction_containing_proteins.fasta")

 # OutputNovelJun(junction_type =  juncs_alt_spl ,genome =   Hsapiens, outfile =  "Splice_junction_proteins.fasta", 
#                proteinseq = proteinseq)


## Step four: Concatenate all fasta files
# -------------------------------------------------------

concatenateFastas(path = output_directory, pattern = "\\.fa", outfile = "combined_test.fasta")




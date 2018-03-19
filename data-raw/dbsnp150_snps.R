# Download snp150Common.txt.gz from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
# Common SNPs(150) - SNPs with >= 1% minor allele frequency (MAF), mapping only once to reference assembly.
# Double click file to unzip

#------
# bash
#------------------------------------------------------------------------------                          
awk -F "\t" '$12 == "single" && \
$9 ~ /^A$|^T$|^C$|^G$/ && \
$2 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ \
{ print $2"\t"$3"\t"$4"\t"$7"\t"$5"\t"$9"\t"$10 }' snp150Common.txt > \
snp150Common.prefiltered.full.txt
#------------------------------------------------------------------------------                                                                


#---
# R
#-----------------------------------------------------------------------------------------------------                               
filtered <- read.table("/dcl01/scharpf1/data/dbruhm/svpipeline/data/snp150Common.prefiltered.full.txt",
                       header = FALSE,
                       stringsAsFactors = FALSE)
colnames(filtered) <- c("chrom", "chromStart", "chromEnd", "strand", "name", "refUCSC", "observed")

# Some SNP entries contain characters not in {A, T, C, G} -- removing these (see below)

# > unique(filtered$observed)
# [1] "C/G"       "A/G"       "G/T"       "C/T"       "A/T"       "A/C"      
# [7] "A/G/T"     "A/C/T"     "A/C/G"     "C/G/T"     "A/C/G/T"   "A/C/G/T/W"
# [13] "A/C/T/W"   "A/C/G/M/T" "A/C/G/N/T"


n.ind <- grep("N", filtered$observed)
m.ind <- grep("M", filtered$observed)
w.ind <- grep("W", filtered$observed)
ind <- unique(c(n.ind, m.ind, w.ind))

# length(ind) is only 4... just drop them.
filtered <- filtered[-(ind),]


# Making a column for only the alt allele (difference between refUCSC and observed)

library(stringr)

filtered$altAllele <- str_replace(string = filtered$observed, 
                                  pattern = filtered$refUCSC, 
                                  replacement = "")

filtered$altAllele <- gsub("//", "/", filtered$altAllele)
filtered$altAllele <- gsub("^/", "", filtered$altAllele)
filtered$altAllele <- gsub("/$", "", filtered$altAllele)

# Convert to GRanges

library(GenomicRanges)

gr <- GRanges(seqnames = filtered$chrom,
              ranges = IRanges(start = filtered$chromEnd, 
                               end = filtered$chromEnd),
              strand = filtered$strand,
              refUCSC = filtered$refUCSC,
              altAllele = filtered$altAllele)
strand(gr) <- "+"

# Making seqinfo match the seqinfo from svfilters.hg19
data(bins1kb, package = "svfilters.hg19")

gr <- keepSeqlevels(gr, seqlevels(bins1kb), pruning.mode = "coarse")
gr <- sortSeqlevels(gr)
genome(gr) <- "hg19"
seqlengths(gr) <- seqlengths(svfilters.hg19::bins1kb)

# Remove any overlapping positions
gr <- gr[which(countOverlaps(gr, gr, ignore.strand = TRUE) == 1)]

# liftOver to hg18

library(rtracklayer)

chain <- import.chain("/dcl01/scharpf1/data/dbruhm/svpipeline/data/hg19ToHg18.over.chain")
gr18 <- unlist(liftOver(gr, chain))
gr18 <- sortSeqlevels(gr18)

# Remove any overlapping positions
gr18 <- gr18[which(countOverlaps(gr18, gr18, ignore.strand = TRUE) == 1)]

# Replacing seqinfo with that from bins1kb in svfilters.hg18
genome(gr18) <- "hg18"
seqlengths(gr18) <- seqlengths(svfilters.hg18::bins1kb)

# Sort
gr18 <- sort(gr18)

# Save object under the name snps

dbsnp150_snps <- gr18
save(dbsnp150_snps, file = "/dcl01/scharpf1/data/dbruhm/svpipeline/data/hg18/dbsnp150_snps.rda")
#-----------------------------------------------------------------------------------------------------
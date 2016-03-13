##library(EnsDb.Hsapiens.v75)  ## what version for hg18?
library(TxDb.Hsapiens.UCSC.hg18.refGene)
txdb <- TxDb.Hsapiens.UCSC.hg18.refGene
tx <- transcripts(TxDb.Hsapiens.UCSC.hg18.refGene)

library(AnnotationHub)
ah <- AnnotationHub()
query(ah, "OrgDb")
hsapiens <- ah[["AH49582"]]
keytypes(hsapiens)
refseqid <- keys(hsapiens, "REFSEQ")
map <- select(hsapiens, refseqid, "SYMBOL", "REFSEQ")
table(tx$tx_name %in% map$REFSEQ)
## 66 of the transcripts are not in the map
tx$tx_name[! tx$tx_name %in% map$REFSEQ]
tx <- tx[ tx$tx_name %in% map$REFSEQ ]
symbols <- setNames(map$SYMBOL, map$REFSEQ)
symbols <- symbols[tx$tx_name]
tx$gene_name <- symbols


library(GenomeInfoDb)
seqlevelsStyle(tx) <- "UCSC"
tx <- keepSeqlevels(tx, paste0("chr", c(1:22, "X", "Y", "M")))
tx <- sort(tx)
si <- seqinfo(tx)
genome(si) <- "hg18"
seqinfo(tx) <- si

df <- read.csv("~/Software/svpackages/svfilters/inst/extdata/cancer_genes_2016-03-05.csv",
               stringsAsFactors=FALSE)
df <- df[!is.na(df$biol_sign), ]
cancer.genes <- df$biol_sign[df$is_clin_sign]
tx$cancer_connection <- tx$gene_name %in% cancer.genes
tx$biol_sign <- tx$gene_name %in% df$biol_sign
transcripts <- tx
save(transcripts, file="~/Software/svpackages/svfilters.hg18/data/transcripts.rda")

library(svfilters.hg18)
library(svcnvs)
germline_filters <- listGenomeFilters()
germline_filters <- germline_filters[-match("transcripts", names(germline_filters))]
save(germline_filters, file="svfilters.hg18/data/germline_filters.rda")

library(rtracklayer)
chain <- import.chain("~/Downloads/hg19ToHg18.over.chain")
load("svfilters.hg19/data/germline_rear.rda")
## lift over to hg18
gr.hg19 <- germline_rear
germline_rear <- liftOver(gr.hg19, chain)

library(BSgenome.Hsapiens.UCSC.hg18)
genome <- BSgenome.Hsapiens.UCSC.hg18
si <- seqinfo(genome)
si2 <- keepStandardChromosomes(si)
germline_rear <- keepSeqlevels(germline_rear, seqlevels(si2))
germline_rear <- unlist(germline_rear)
seqinfo(germline_rear) <- si2
save(germline_rear, file="svfilters.hg18/data/germline_rear.rda")

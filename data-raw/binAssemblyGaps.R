library(GenomicRanges)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg18)
slevels <- c(paste0("chr", 1:22), "chrX", "chrY")
bins <- tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg18)[slevels, ],
                   tilewidth=1000,
                   cut.last.tile.in.chrom=TRUE)
##
## Assembly gaps / GC extreme
##
gc <- GCcontent(Hsapiens, bins)
binAssemblyGaps_hg18 <- bins[as.numeric(gc) < 0.1]
sum(width(binAssemblyGaps_hg18))/1e6 ## 234o Mb
save(binAssemblyGaps_hg18, file="../data/binAssemblyGaps_hg18.rda")

bins$gc <- gc

bins10kb <- tileGenome(seqinfo(BSgenome.Hsapiens.UCSC.hg18)[slevels, ],
                       tilewidth=10000,
                       cut.last.tile.in.chrom=TRUE)a
hits <- findOverlaps(bins10kb, bins)
hit.list <- split(bins, queryHits(hits))
hit.list <- hit.list[elementNROWS(hit.list) == 10]
tmp <- unlist(hit.list)
gc <- tmp$gc
f <- factor(rep(seq_along(hit.list), each=10))
gclist <- split(gc, f)
gc.mns <- sapply(gclist, sum)/10

bins10kb <- bins10kb[as.integer(names(hit.list))]
bins10kb$gc <- as.integer(round(gc.mns*1000, 0))
bins10kb <- bins10kb[bins10kb$gc >= 0.1]
save(bins10kb, file="/dcl01/scharpf/data/svpackages/svfilters.hg18/data/bins10kb.rda")

mappability <- function(map, bins){
  ## summarize mappability by bin
  score <- rep(NA, length(bins))
  o <- findOverlaps(bins, map)
  if(length(o)==0) return(score)
  j <- subjectHits(o)
  i <- queryHits(o)
  subjectIndex <- split(j, i)
  w <- width(map)
  s <- map$score
  ## weight the mappability score by the width of the mappability interval
  avg <- sapply(subjectIndex, function(i, score, weight){
    (sum(score[i]*weight[i],na.rm=TRUE))/sum(weight[i], na.rm=TRUE)
  }, score=s, weight=w)
  score[as.integer(names(avg))] <- avg
  as.integer(score*1000)
}
path <- "/dcl01/scharpf/data/reference/no_alternates/hg18"
file <- file.path(path, "wgEncodeCrgMapabilityAlign100mer.bw")
map.bw <- import.bw(file)
bins10kb$map <- mappability(map.bw, bins10kb)
save(bins10kb, file="/dcl01/scharpf/data/svpackages/bins10kb.rda")

genomeGaps <- function() {
  library(BSgenome.Hsapiens.UCSC.hg18)
  genome <- "hg18"
  mySession <- browserSession()
  genome(mySession) <- genome
  gaps <- getTable(ucscTableQuery(mySession, track="gap"))
  gaps.gr <- GRanges(gaps$chrom, IRanges(gaps$chromStart,
                                         gaps$chromEnd),
                     type=gaps$type)
  gaps_hg18 <- gaps.gr
  ## no gaps available in mitochondria
  gaps_hg18 <- keepSeqlevels(gaps_hg18, paste0("chr", c(1:22, "X", "Y")), pruning.mode = "coarse")
  seqinfo(gaps_hg18) <- seqinfo(Hsapiens)[seqlevels(gaps_hg18), ]
  gaps_hg18 <- trim(gaps_hg18)
  
  # Reduce the two chrY centromere ranges
  ycen <- gaps_hg18[which((gaps_hg18$type == "centromere") & (seqnames(gaps_hg18) == "chrY"))]
  combined <- range(ycen)
  combined$type <- "centromere"
  gaps_hg18 <- gaps_hg18[-which((gaps_hg18$type == "centromere") & (seqnames(gaps_hg18) == "chrY"))]
  gaps_hg18 <- c(gaps_hg18, combined)
  gaps <- sort(gaps_hg18)
  save(gaps_hg18, file="/Users/dbruhm/Desktop/cancer-genomics/svfilters.hg18/data/gaps.rda")
}


## Add the OncoKB-derived cancer_gene and clinically_significant columns to
## data/transcripts.rda.
##
## Run from the package root, after data-raw/transcripts.R (which sources this
## file at its end), or on its own to re-annotate the committed object without
## rebuilding the transcript set:
##
##   Rscript data-raw/cancer_gene_columns.R
##
## The gene lists come from CancerGenes (cancer-genomics/CancerGenes), the lab's
## single definition of a cancer gene, shared with svfilters.hg19. CancerGenes
## is a private repository -- it holds the parsed OncoKB tables, which MSK's
## permission does not cover -- so this step can be rerun only inside the lab.
## It is a build-time dependency only: the result ships in
## data/transcripts.rda, and the installed package never loads CancerGenes.
##
## OncoKB (https://www.oncokb.org), maintained by Memorial Sloan Kettering
## Cancer Center, is the source of the OncoKB content in these two columns,
## redistributed here with MSK's permission (2026-09-02) under the OncoKB Terms
## of Use (https://www.oncokb.org/terms). The raw OncoKB downloads are never
## part of this repository. See LICENSE.note.

stopifnot("run this script from the package root" = file.exists("DESCRIPTION"))
suppressMessages({
  library(GenomicRanges)
  library(CancerGenes)
})
stopifnot(packageVersion("CancerGenes") >= "0.2.0")

path <- "data/transcripts.rda"
env <- new.env(parent = emptyenv())
load(path, envir = env)
tx <- env[["transcripts"]]

## Idempotent: drop any previous annotation before adding it again.
mcols(tx)$clinically_significant <- NULL
mcols(tx)$cancer_gene <- NULL
stopifnot(identical(names(mcols(tx)),
                    c("tx_id", "tx_name", "gene_name",
                      "cancer_connection", "biol_sign")))

## clinical_aliases = TRUE: this transcript set uses 2016-era OrgDb symbols,
## so H3F3A, WHSC1, HIST1H3B, HIST1H3C and MRE11A would otherwise miss the
## levels OncoKB records under H3-3A, NSD2, H3C2, H3C3 and MRE11.
transcripts <- add_cancer_gene_columns(tx, clinical_aliases = TRUE)

## Pin the result, so a changed snapshot or gene definition cannot slip in
## without these counts being updated deliberately. The OncoKB snapshot counts
## themselves (825 ONCOGENE/TSG symbols, 2,259 aliases) are asserted where the
## snapshot is parsed, in CancerGenes/data-raw/01_oncokb.R.
n_symbols <- function(x) length(unique(transcripts$gene_name[x]))
clinical <- !is.na(transcripts$clinically_significant)
stopifnot(
  length(transcripts) == 51491L,
  sum(transcripts$cancer_gene) == 2919L,
  n_symbols(transcripts$cancer_gene) == 1078L,
  sum(clinical) == 505L,
  n_symbols(clinical) == 170L,
  all(transcripts$cancer_gene[clinical]),
  identical(unique(genome(transcripts)), "hg18"),
  identical(format(metadata(transcripts)$oncokb$snapshots$date),
            c("2026-01-06", "2026-01-12"))
)

save(transcripts, file = path, compress = "xz")

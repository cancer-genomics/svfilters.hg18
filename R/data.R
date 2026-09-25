#' 1kb bins with metadata for GC and mappability
#'
#'  1kb bins with percent GC and average mappability stored in the
#' \code{GRanges} metadata.
#'
#' @docType data
#' @keywords datasets
#' @name bins
#' @usage data(bins1kb)
#' @aliases bins1kb
#'
#' @format a \code{GRanges} object containing 1kb intervals with
#'   metadata columns 'gc' and 'map'.  See details
#'
#' @details
#'
#' The metadata column 'gc' contains the percent GC for the 1kb
#'   intervals.  The percent GC was multiplied by 1000 and rounded to
#'   the nearest integer.
#'
#' The metadata column 'map' contains the average mappability for the
#'   1kb bins.  The average mappability was multiplied by 1000 and
#'   rounded to the nearest integer.
#'
#' @seealso See
#'   \url{ftp://hgdownload1.cse.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig}
#'   for 100mer mappability tracks.
#'
NULL

#' Genome assembly gaps 
#'
#' A \code{GRanges} object representing assembly gaps.  Assembly gaps
#' were defined as 1kb bins with low GC (less than 10
#' percent), including 'N' sequences.  The function \code{reduce} was
#' applied to the 1kb bins with low GC.
#'
#' @docType data
#' @name assembly_gaps
#' @usage data(assembly_gaps)
#' @aliases assembly_gaps
#'
#' @format a \code{GRanges} object with metadata element 'gc'
#'   providing the corresponding GC content for the interval.
#'
NULL

## Low mappability regions
##
## A \code{GRanges} object containing intervals of mappability less
## than 0.75.
##
## @docType data
## @name low_mappability
## @usage data(low_mappability)
## @aliases low_mappability
##
## @format a \code{GRanges} object
## NULL

#' Transcripts with approved HGNC symbols and cancer-gene annotation
#'
#' RefSeq transcripts from \code{TxDb.Hsapiens.UCSC.hg18.refGene}, restricted to
#' chr1-22, X, Y, and M, with \code{gene_name} set to the approved HGNC symbol
#' mapped from the RefSeq accession. NCBI build 36.
#'
#' Two sets of cancer-gene annotation are carried, each assigned by gene
#' symbol:
#'
#' \describe{
#'   \item{\code{cancer_connection}, \code{biol_sign}}{From a literature-based
#'     gene list compiled 2016-03-05
#'     (\code{inst/extdata/cancer_genes_2016-03-05.csv}), with no OncoKB
#'     content. \code{biol_sign} marks genes of biological interest (1,407
#'     symbols) and \code{cancer_connection} the clinically significant subset
#'     (195 symbols).}
#'   \item{\code{clinically_significant}}{For 170 symbols, the highest OncoKB
#'     evidence level within each class (Tx, Dx, Px, R), comma-separated, for
#'     example \code{"Dx1,Px1,R1,Tx1"}; \code{NA} otherwise. Derived entirely
#'     from OncoKB.}
#'   \item{\code{cancer_gene}}{\code{TRUE} for 1,078 symbols: the union of the
#'     OncoKB ONCOGENE/TSG gene list and its aliases, the OncoKB biomarker genes,
#'     and six published driver-gene sets (PMIDs 23539594, 24132290, 24390350,
#'     29056346, 29625053, 32015527).}
#' }
#'
#' Neither OncoKB-derived column carries variant or alteration content, drugs,
#' tumor types, or descriptive text. Five symbols (\code{H3F3A}, \code{WHSC1},
#' \code{HIST1H3B}, \code{HIST1H3C}, \code{MRE11A}) carry the older names of
#' OncoKB genes and are matched to them through OncoKB aliases.
#'
#' @section OncoKB snapshot:
#' The OncoKB content is a \strong{frozen snapshot}: the cancer gene list of
#' 2026-01-06 and the biomarker levels (Tx/Dx/Px/R) of 2026-01-12. OncoKB is
#' updated continuously, so these annotations may be out of date. Anyone
#' needing current OncoKB content should obtain it directly from
#' \url{https://www.oncokb.org} under their own registration and terms. The
#' snapshot dates, sources, and citations are also stored with the object, in
#' \code{metadata(transcripts)$oncokb}.
#'
#' @section Attribution:
#' OncoKB (\url{https://www.oncokb.org}), a precision oncology knowledge base
#' maintained by Memorial Sloan Kettering Cancer Center (MSK), is the source of
#' the OncoKB content in \code{clinically_significant} and \code{cancer_gene}.
#' It is redistributed here with the permission of MSK and used under the
#' OncoKB Terms of Use (\url{https://www.oncokb.org/terms}). MSK makes no
#' warranties or representations with respect to the OncoKB content, and it is
#' not a substitute for professional medical judgment or advice. See
#' \code{LICENSE.note} and \code{citation("svfilters.hg18")}.
#'
#' @references
#' Chakravarty D, Gao J, Phillips SM, et al. OncoKB: A Precision Oncology
#' Knowledge Base. \emph{JCO Precis Oncol.} 2017;2017:PO.17.00011.
#' doi:10.1200/PO.17.00011. PMID: 28890946.
#'
#' Suehnholz SP, Nissan MH, Zhang H, et al. Quantifying the Expanding Landscape
#' of Clinical Actionability for Patients with Cancer. \emph{Cancer Discov.}
#' 2024;14(1):49-65. doi:10.1158/2159-8290.CD-23-0467. PMID: 37849038.
#'
#' @docType data
#' @keywords datasets
#' @name transcripts
#' @usage data(transcripts)
#' @aliases transcripts
#' @format a \code{GRanges} object with metadata columns \code{tx_id},
#'   \code{tx_name}, \code{gene_name}, \code{cancer_connection},
#'   \code{biol_sign}, \code{clinically_significant}, and \code{cancer_gene},
#'   and the OncoKB provenance record in \code{metadata(transcripts)$oncokb}
#'
#' @examples
#' data(transcripts)
#' metadata(transcripts)$oncokb$snapshots
#' unique(transcripts$gene_name[!is.na(transcripts$clinically_significant)])
NULL

#' Genome gaps downloaded from UCSC
#'
#' Includes heterochromatin, centromeres, and telomeres.
#'
#' @docType data
#' @keywords datasets
#' @name gaps
#' @usage data(gaps)
#' @aliases gaps
#' @format a \code{GRanges} object
#'
#' @examples
#' data(gaps)
NULL

#' Coverage filters
#'
#' A \code{GRangesList} of germline deletions, amplified regions, and copy
#' number outliers identified in 10 lymphoblast cell lines and 8 blood samples
#' from ovarian cancer patients (matched-normal).
#'
#' @docType data
#' @keywords datasets
#' @name coverage_filters
#' @usage data(coverage_filters)
#' @aliases coverage_filters
#' @format a \code{GRangesList} object
#'
#' @examples
#' data(coverage_filters)
#' sum(width(coverage_filters))/1e6
NULL

#' Germline rearrangements
#'
#' A reduced \code{GRanges} object of rearrangements identified in 10
#' lymphoblast cell lines and 8 blood samples from ovarian cancer patients.
#'
#' @docType data
#' @keywords datasets
#' @name germline_rear
#' @usage data(germline_rear)
#' @aliases germline_rear
#' @format a reduced \code{GRanges} object
#'
#' @examples
#' data(germline_rear)
NULL


#' A set of 1 million frequently heterozygous SNPs from dbSNP build 150
#'
#' A \code{GRanges} object of SNP positions from dbSNP build 150 lifted over
#' from hg19. Metadata columns for the reference allele and alternate allele are included.  
#'
#' @docType data
#' @keywords datasets
#' @name snps
#' @usage data(snps)
#' @aliases snps
#' @format a \code{GRanges} object
#'
#' @examples
#' data(snps)
NULL


#' A comprehensive set of SNPs from dbSNP build 150 
#'
#' A \code{GRanges} object of SNP positions from dbSNP build 150.
#' Metadata columns for the reference allele and alternate allele are included.
#' 
#' @docType data
#' @keywords datasets
#' @name dbsnp150_snps
#' @usage data(dbsnp150_snps)
#' @aliases dbsnp150_snps
#' @format a \code{GRanges} object
#'
#' @examples
#' data(dbsnp150_snps)
NULL

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

#' Low mappability regions
#'
#' A \code{GRanges} object containing intervals of mappability less
#' than 0.75.
#'
#' @docType data
#' @name low_mappability
#' @usage data(low_mappability)
#' @aliases low_mappability
#'
#' @format a \code{GRanges} object
NULL


#' Lymphoblast coverage-based filters
#'
#' Genomic intervals for outliers, deletions, and amplifications
#' identified in lymphoblast cell lines.
#'
#' @details
#'
#' Coverage estimates for 10 lymphoblast cell lines were developed
#'   using the package \code{svpreprocess} in non-overlapping 1kb bins
#'   along the genome.  The autosomal median absolute deviation (MAD)
#'   was used as a robust measure of variance for log2-transformed
#'   counts. Log2-transformed counts were adjusted for GC content, as
#'   well as a loess smoother of background coverage estimates.  Bins
#'   for which two or more of the 10 lymphoblast cell lines had a
#'   preprocessed coverage estimate 5 or more MADs from zero were
#'   categorized as outliers.  Segmentation of the preprocessed
#'   coverage estimates was performed using circular binary
#'   segmentation with default settings (see \code{svcnvs} package).
#'   For purposes of an ad-hoc germline filter, segments with means
#'   less than -1 or greater than 1 were considered deletions or
#'   amplicons, respectively.  The genomic intervals for outliers,
#'   deletions, and amplicons were combined in a single object.
#'
#' @seealso \code{\link{bins1kb}}
#'
#' @docType data
#' @name lymphoblast_filters
#' @usage data(lymphoblast_filters)
#' @aliases lymphoblast_filters
#' @format a \code{GRanges} object
NULL


#' Transcripts with approved HGNC symbols 
#'
#' Build hg19 corresonds to Ensembl build 75 and NCBI build 37.
#'
#' @docType data
#' @keywords datasets
#' @name transcripts
#' @usage data(transcripts)
#' @aliases transcripts
#' @format a \code{GRanges} object
#'
#' @examples
#' data(transcripts)
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

#' Lymphoblast rearrangement intervals
#'
#' Genomic intervals demarcating clusters of reads involved in
#' improper pairs were identified in a set of 10 lymphoblast cell
#' lines.  The genomic intervals that could could be linked by 5 or
#' more read pairs were reduced.  These reduced intervals are used to
#' remove potential germline rearrangements.
#'
#' @docType data
#' @keywords datasets
#' @name lymphoblast_rear
#' @usage data(lymphoblast_rear)
#' @aliases lymphoblast_rear
#' @format a \code{GRanges} object
#'
#' @examples
#' data(lymphoblast_rear)
#' sum(width(lymphoblast_rear))/1e6
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

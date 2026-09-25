# svfilters.hg18

Genomic filters and annotation for structural-variant analysis on hg18 (NCBI36),
used by [trellis](https://github.com/cancer-genomics/trellis). The package ships
`GRanges` objects of regions to exclude from analysis (assembly gaps, germline
deletions and amplifications, germline and normal-blood rearrangements, dbSNP
positions). It also ships a `transcripts` object used to annotate the deletions
and amplifications that trellis calls. The cancer-gene annotation is the same as
in [svfilters.hg19](https://github.com/cancer-genomics/svfilters.hg19).

```r
# install.packages("remotes")
remotes::install_github("cancer-genomics/svfilters.hg18")
library(svfilters.hg18)
data(transcripts)
?transcripts
```

## Cancer-gene annotation on `transcripts`

`transcripts` has RefSeq transcripts keyed by HGNC symbol and two sets of
cancer-gene annotation:

| Column | Content | Source |
|---|---|---|
| `cancer_connection`, `biol_sign` | clinically significant (195 symbols) / biologically significant (1,407 symbols) genes | literature-based list compiled 2016-03-05 |
| `clinically_significant` | highest OncoKB evidence level per class (Tx, Dx, Px, R) for 170 symbols; `NA` otherwise | OncoKB |
| `cancer_gene` | 1,078 symbols: OncoKB cancer genes and aliases, OncoKB biomarker genes, and six published driver-gene sets | OncoKB and literature |

### OncoKB

**OncoKB is the source of the OncoKB content in this package.** OncoKB
(<https://www.oncokb.org>) is a precision oncology knowledge base maintained by
Memorial Sloan Kettering Cancer Center (MSK). The `clinically_significant`
column and the OncoKB part of `cancer_gene` are redistributed here with MSK's
permission. They are used under the OncoKB Terms of Use:
<https://www.oncokb.org/terms>.

**This is a frozen snapshot.** It uses the OncoKB cancer gene list of
2026-01-06 and the biomarker levels of 2026-01-12. OncoKB is updated
continuously, so these annotations may be out of date. For current OncoKB
content, get it directly from <https://www.oncokb.org> under your own
registration and terms. The snapshot dates are:

- stored with the object, in `metadata(transcripts)$oncokb`
- printed when the package is attached.

Only gene symbols and level labels are included. There are no variants, drugs,
tumor types or descriptive text, and the raw OncoKB downloads are not part of
this package.

MSK makes no warranties or representations with respect to the OncoKB content.
It is not a substitute for professional medical judgment or advice.

If you use these annotations, please cite OncoKB (also available from
`citation("svfilters.hg18")`):

- Chakravarty D, Gao J, Phillips SM, et al. OncoKB: A Precision Oncology
  Knowledge Base. *JCO Precis Oncol.* 2017;2017:PO.17.00011.
  doi:10.1200/PO.17.00011. PMID: 28890946.
- Suehnholz SP, Nissan MH, Zhang H, et al. Quantifying the Expanding Landscape
  of Clinical Actionability for Patients with Cancer. *Cancer Discov.*
  2024;14(1):49-65. doi:10.1158/2159-8290.CD-23-0467. PMID: 37849038.

### Literature driver-gene sets

`cancer_gene` also includes genes from:

- Vogelstein et al. 2013 (PMID 23539594)
- Kandoth et al. 2013 (PMID 24132290)
- Lawrence et al. 2014 (PMID 24390350)
- Martincorena et al. 2017 (PMID 29056346)
- Bailey et al. 2018 (PMID 29625053)
- Dietlein et al. 2020 (PMID 32015527)

## Licensing

The package code is Artistic-2.0. The bundled annotation includes third-party
content with its own terms; see [`LICENSE.note`](LICENSE.note).

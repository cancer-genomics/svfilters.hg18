# svfilters.hg18 0.0.20

- `transcripts` now carries the `clinically_significant` and `cancer_gene`
  columns, matching `svfilters.hg19` (>= 0.0.27). They sit alongside
  `cancer_connection` and `biol_sign`, which are unchanged, so code written
  against either schema keeps working.
- `clinically_significant` is derived from OncoKB, and `cancer_gene` in part
  from OncoKB. Both are redistributed with the permission of Memorial Sloan
  Kettering Cancer Center. They are a frozen snapshot: the OncoKB cancer gene
  list of 2026-01-06 and the biomarker levels of 2026-01-12.
- The snapshot dates are recorded in three places: in `?transcripts`, in
  `metadata(transcripts)$oncokb`, and in a message printed when the package is
  attached. See `LICENSE.note` and `citation("svfilters.hg18")`.
- The transcript set itself is unchanged from 0.0.19: 51,491 RefSeq
  transcripts, genome hg18.
- `?transcripts` now describes the object's actual source (RefSeq via
  `TxDb.Hsapiens.UCSC.hg18.refGene`) instead of Ensembl build 75.
- New files: README, `LICENSE.note`, `inst/CITATION`, and tests for
  `transcripts`.

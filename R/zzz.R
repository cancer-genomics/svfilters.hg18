## Startup message recording the OncoKB snapshot, as committed to MSK when
## permission to redistribute the OncoKB-derived columns was requested. The text
## is CancerGenes::oncokb_snapshot_message("svfilters.hg18"), pasted in as a
## fixed string so the installed package does not depend on CancerGenes;
## tests/testthat/test_transcripts.R checks it against metadata(transcripts).
oncokb_message <- paste(
  "transcripts carries OncoKB-derived annotation (cancer gene list 2026-01-06;",
  "biomarker levels (Tx/Dx/Px/R) 2026-01-12), a frozen snapshot. OncoKB",
  "(https://www.oncokb.org) is maintained by MSK; see",
  "citation(\"svfilters.hg18\") and https://www.oncokb.org/terms. For current",
  "OncoKB content, obtain it from https://www.oncokb.org under your own",
  "registration.")

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(oncokb_message)
}

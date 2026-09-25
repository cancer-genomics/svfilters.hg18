context("transcripts")

test_that("transcripts carries both cancer-gene schemas", {
  data(transcripts, package = "svfilters.hg18", envir = environment())
  expect_identical(names(mcols(transcripts)),
                   c("tx_id", "tx_name", "gene_name", "cancer_connection",
                     "biol_sign", "clinically_significant", "cancer_gene"))
  expect_identical(unique(genome(transcripts)), "hg18")
  expect_identical(length(unique(transcripts$gene_name[transcripts$cancer_gene])),
                   1078L)
  clinical <- !is.na(transcripts$clinically_significant)
  expect_identical(length(unique(transcripts$gene_name[clinical])), 170L)
  expect_true(all(transcripts$cancer_gene[clinical]))
  expect_false(any(transcripts$clinically_significant[clinical] == ""))
  expect_identical(length(unique(transcripts$gene_name[transcripts$cancer_connection])),
                   195L)
})

test_that("OncoKB provenance travels with the object and the startup message", {
  data(transcripts, package = "svfilters.hg18", envir = environment())
  prov <- metadata(transcripts)$oncokb
  expect_identical(prov$terms, "https://www.oncokb.org/terms")
  dates <- format(prov$snapshots$date)
  expect_identical(dates, c("2026-01-06", "2026-01-12"))
  msg <- svfilters.hg18:::oncokb_message
  for (d in dates) expect_true(grepl(d, msg, fixed = TRUE))
  expect_true(grepl(prov$terms, msg, fixed = TRUE))
  expect_true(grepl("frozen snapshot", msg, fixed = TRUE))
})

library(ovariantn)
dp <- projectOvariantn()
data(germline_wbc)
rfiles <- file.path(dp["rearrangements/germline"], paste0(germline_wbc, ".rds"))
rear_list <- lapply(rfiles, readRDS)
lb <- lapply(rear_list, linkedBins)
lt <- lapply(lb, function(x) x$linked.to)
lb <- GRangesList(lb)
lt <- GRangesList(lt)
lb <- unlist(lb)
lt <- unlist(lt)
## takes the place lymphoblast_rear_hg18
normalblood_rear <- reduce(c(granges(lb), lt))
save(normalblood_rear, file="~/Software/svpackages/svfilters.hg18/data/normalblood_rear.rda")

library(svfilters.hg18)
##data(lymphoblast_rear) ## seqinfo is not hg18
##grl <- GRangesList(lymphoblast=lymphoblast_rear,
##                   normalblood=normalblood_rear)


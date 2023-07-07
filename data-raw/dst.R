# data notes
# -----------
# open metrics.R and run everything
# then, in the same R session, run this file

dst <- rbindlist(lapply(single_results, function(x){
  # extract results of the pairwise tests
  p <- x$pairwise_tests
  p$phenotype <- x$phenotype
  p$study <- x$study
  p$ancestry <- x$ancestry
  p$bbid <- x$bbid
  p$method_x <- adjust_methods_levels(p$method_x)
  p$method_y <- adjust_methods_levels(p$method_y)
  return(p)
}), fill=TRUE, use.names=TRUE)

dst <- rename_phenotypes(dst)

dst[,ci_low_diff:=beta_diff-1.96*se_diff]
dst[,ci_high_diff:=beta_diff+1.96*se_diff]

# sample overlap
dst <- dst[!(bbid == 'ukbb' & phenotype == 'Height' & ancestry == 'EUR')]
dst <- dst[!(bbid == 'ukbb' & phenotype == 'AD' & ancestry == 'EUR')]
dst <- dst[!(bbid == 'ebb' & phenotype == 'eGFR')]

# let R pick the best compression scheme
# tools::resaveRdaFiles("data/", compress="auto")
# tools::checkRdaFiles("data/") -> xz for this dataset
usethis::use_data(dst, overwrite=TRUE, compress="xz")

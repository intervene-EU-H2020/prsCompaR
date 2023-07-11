# pv_mrg
# like dst, but meta-analysed for effect size estimates
# run order: metrics.R, meta_res.R, this file (in same R session)

pv_mrg <-  rbind(pv_eur, pv_sas)
setnames(pv_mrg, c('rn', 'iy'), c('method_x', 'method_y'))
pv_mrg[,c('method_x','method_y'):=list(gsub('^method','',method_x),gsub('^method','',method_y))]
pv_mrg[,c('method_type_x','method_type_y'):=list(gsub('^.+_','',method_x),gsub('^.+_','',method_y))]
pv_mrg[,c('method_x','method_y'):=list(gsub('_.*$','',method_x),gsub('_.*$','',method_y))]
pv_mrg <- rename_phenotypes(pv_mrg)
pv_mrg$method_x <- adjust_methods_levels(pv_mrg$method_x)
pv_mrg[,ci_low_diff:=beta_diff - 1.96*se_diff]
pv_mrg[,ci_high_diff:=beta_diff + 1.96*se_diff]

# export processed data --------------------------------------------------------
# let R pick the best compression scheme
# tools::resaveRdaFiles("data/", compress="auto")
# tools::checkRdaFiles("data/") -> xz for this dataset
usethis::use_data(pv_mrg, overwrite = TRUE, version=2, compress="xz")


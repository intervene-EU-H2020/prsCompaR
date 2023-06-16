get_methods_levels <- function(exclude_multiprs=FALSE) {
  methods_levels <- c('pt.clump','dbslmm','ldpred2','lassosum','megaprs',
                      'prscs','sbayesr','UKBB.EnsPRS')
  if (exclude_multiprs) {
    return(methods_levels[-length(methods_levels)])
  } else {
    return(methods_levels)
  }
}

adjust_methods_levels <- function(x, exclude_multiprs=FALSE){
  x <- gsub('pt\\.clump\\.nested', 'pt.clump', x)
  x[grepl('UKBB-[A-Z]+\\.MultiPRS',x)] <- 'UKBB.EnsPRS'
  x <- factor(x, levels=get_methods_levels(exclude_multiPRS))
  if(any(is.na(x))){
    stop()
  }
  return(x)
}

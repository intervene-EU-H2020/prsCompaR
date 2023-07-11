# run after metrics.R (in the same session)

library(Matrix)
library(metafor)

# parameters -------------------------------------------------------------------

meta_analysis_use_test <- 't'
if (meta_analysis_use_test == 'z'){
  meta_analysis_dfs <- 'residual'
} else {
  stopifnot(meta_analysis_use_test == 't')
  meta_analysis_dfs <- 'contain'
}

methods_type_levels <- c('pt.clump_auto','pt.clump_CV','dbslmm_auto','ldpred2_auto','ldpred2_CV','lassosum_auto','lassosum_CV','megaprs_auto','megaprs_CV','prscs_auto','prscs_CV','sbayesr_auto','UKBB.EnsPRS_CV')

# ------------------------------------------------------------------------------

single_results <- lapply(1:nrow(results_metadata), loader_fun)

setnames(results_metadata, 'pheno', 'phenotype', skip_absent = TRUE)
results_metadata$i <- 1:nrow(results_metadata)

study_phenotype_combinations <- unique(results_metadata[,list(study,phenotype)])

# exclusions
metadata_ma <- data.table::copy(results_metadata)
# sample overlap
metadata_ma <- metadata_ma[!(bbid == 'ukbb' & phenotype == 'Height' & ancestry == 'EUR')]
metadata_ma <- metadata_ma[!(bbid == 'ukbb' & phenotype == 'Alzheimers_disease' & ancestry == 'EUR')]
metadata_ma <- metadata_ma[!(bbid == 'ebb' & phenotype == 'Creatinine_eGFR')]

# outliers / poor performance
metadata_ma <- metadata_ma[!(bbid == 'gnh' & phenotype %in% c('T1D','HbA1c'))]
metadata_ma <- metadata_ma[!(bbid == 'hunt' & phenotype %in% c('T1D'))]

get_block_correl_matrix <- function(r){
  return(as.matrix(bdiag(lapply(r, '[[','scor'))))
}

get_v <- function(r){

  V <- diag(r$metrics$SD) %*% r$scor %*% diag(r$metrics$SD)
  return(V)
}

results_list <- list()

rma_mv_dfs <- meta_analysis_dfs
rma_mv_test <- meta_analysis_use_test


pairwise_z_tests <- function(m_mv, calculate_t=FALSE){
  # for meta-analysis, takes an object of class rmv.mv, returns pairwise differences between effect size estimates
  # adjusts for the correlation between effect size estimates
  # performs z-tests

  c <- coef(m_mv)
  dst <- data.table(as.matrix(dist(c, method='manhattan')), keep.rownames = T)
  mvcov <- data.table(data.frame(vcov(m_mv, type = "fixed")), keep.rownames = T)
  d <- diag(vcov(m_mv))
  names(d) <- row.names(vcov(m_mv))
  dst <- melt(dst, id.vars = 'rn', variable.name = 'iy', value.name='beta_diff')
  mvcov <- melt(mvcov, id.vars = 'rn', variable.name = 'iy', value.name='cov_x_y')
  dst <- merge(dst, mvcov)

  dst[,var_x:=d[rn]]
  dst[,var_y:=d[iy]]
  dst[,beta_x:=coef(m_mv)[rn]]
  dst[,beta_y:=coef(m_mv)[iy]]
  dst$beta_diff <- with(dst, ifelse(beta_x > beta_y, beta_diff, -beta_diff))
  dst[,se_diff:=sqrt(var_x+var_y-2*cov_x_y)]
  dst[,z_diff:=beta_diff/se_diff]

  dst[,p_z_greater:=pnorm(-z_diff)]
  dst[,p_z_smaller:=pnorm(z_diff)]
  dst[,p_z_equal:=2*pnorm(abs(z_diff), lower.tail = F)]


  if (calculate_t){

    # hack the anova function to calculate the same thing we did above, but use t-tests
    # t-tests might be stringent, z-tests might be permissive...

    n_coef <- m_mv$p
    coef_names <- names(d)

    X <- do.call(rbind, lapply(1:n_coef, function(x){
      # defines the hypothesis we wish to test
      # returns all pairwise combinations of coeffiecients ( coef_a - coef_b )
      X <- matrix(data = rep(0,n_coef^2), nrow = n_coef, ncol = n_coef)
      diag(X) <- -1
      X[,x] <- 1
      return(X)
    }))

    tmp <- anova(m_mv, X = X)

    result <- data.table(beta_diff_t = tmp$Xb[,1], df_t=tmp$ddf, tval=tmp$zval, p_t_equal = tmp$pval, p_t_smaller = pt(tmp$zval, df=tmp$ddf[1], lower.tail = T), p_t_greater = pt(tmp$zval, df=tmp$ddf[1], lower.tail = F) )

    result$rn <- coef_names[rep(1:n_coef, each=n_coef)]
    result$iy <- coef_names[rep(1:n_coef, times=n_coef)]

    result <- result[iy != rn]

    dst <- merge(dst, result, by = c('rn','iy'), all.x = T)

  }


  return(dst)

}

for (i in 1:nrow(study_phenotype_combinations)){


  # select the results we want to meta-analyse
  p <- study_phenotype_combinations$phenotype[i]
  s <- study_phenotype_combinations$study[i]

  select_i <- metadata_ma[phenotype == p & study == s]$i
  res <- single_results[select_i]

  testdata <- rbindlist(lapply(res, function(x){
    # get the metrics
    m <- x$metrics
    m$study <- x$study
    m$phenotype <- x$phenotype
    m$ancestry <- x$ancestry
    m$bbid <- x$bbid
    return(m)
  }))

  testdata$method <- paste(testdata$method, testdata$method_type, sep='_')
  testdata$method <- adjust_methods_levels_2(testdata$method)

  # variance covariance matrix of the effect-size estimates
  V <- as.matrix(bdiag(lapply(res, get_v)))

  ancestries <- unique(testdata$ancestry)

  # estimate the effects of methods and biobank within ancestries
  ancestry_specific_results <- lapply(ancestries, function(x){

    # subset the data to the selected ancestry
    selection <- testdata$ancestry == x

    # subset the data and exclude pt+clump
    selection_no_ptclump <- (testdata$ancestry == x) & !grepl('pt.clump', testdata$method)
    testdata_no_ptclump <- testdata[selection_no_ptclump]
    V_no_ptclump <- V[selection_no_ptclump, selection_no_ptclump]

    # subset the data and exclude ensemble and pt+clump
    selection_no_ptclump_no_ensemble <- (testdata$ancestry == x) & !(grepl('pt.clump', testdata$method) | (grepl('MultiPRS|EnsPRS', testdata$method)))
    testdata_no_ptclump_no_ensemble <- testdata[selection_no_ptclump_no_ensemble]
    V_no_ptclump_no_ensemble <- V[selection_no_ptclump_no_ensemble, selection_no_ptclump_no_ensemble]

    if(any(grepl('pt.clump', testdata_no_ptclump_no_ensemble$method))){
      stop(print(testdata_no_ptclump_no_ensemble$method))
    }


    testdata <- testdata[selection]
    V <- V[selection, selection]

    # only perform the meta-analysis if there are results for more than 1 biobank
    bbids <- unique(testdata$bbid)

    if (length(bbids) > 1){

      # a model that jointly estimates the biobank and method effects (not used)
      # m_mv_all <- rma.mv(yi = BETA, V = V, mods = ~ method + bbid, dfs=rma_mv_dfs, data=testdata, test=rma_mv_test)
      # a model that only estimates the method effects, and treats the biobanks as a random effect, without the intercept (baseline)
      m_mv_method <-  tryCatch( {
        rma.mv(yi = BETA, V = V, mods = ~ method - 1, random = list(~ 1 | bbid), dfs=rma_mv_dfs, data=testdata, test=rma_mv_test, cvvc=F)
      },
      error = function(cond){
        rma.mv(yi = BETA, V = V, mods = ~ method - 1, random = list(~ 1 | bbid), dfs=rma_mv_dfs, data=testdata, test=rma_mv_test, cvvc=F, control=list(rel.tol=1e-8, iter.max=1000))
      })

      # pairwise z-tests for the model parameters
      pv_method <- pairwise_z_tests(m_mv_method, calculate_t = T)


      return(list(
        m_mv_method=m_mv_method,
        pv_method=pv_method
      ))

    } else {

      warning(paste0('Phenotype ', p, ' only available for one biobank (',bbids,') for ancestry ',x))

      m_mv_method <- rma.mv(yi = BETA, V = V, mods = ~ method - 1, dfs=rma_mv_dfs, data=testdata, test=rma_mv_test)
      pv_method <- pairwise_z_tests(m_mv_method, calculate_t = T)

      return(list(
        m_mv_method=m_mv_method,
        pv_method=pv_method
      ))

    }

  })

  names(ancestry_specific_results) <- ancestries

  results <- list()
  results$testdata <- testdata
  results$V <- V
  results$study <- s
  results$phenotype <- p
  results$ancestry_specific_results <- ancestry_specific_results

  results_list[[paste(s,p,sep=':')]] <- results
}

extract_methods_results <- function(x, ancestry, model='m_mv_method'){
  if (length(x$ancestry_specific_results[[ancestry]])>0){

    r <- x$ancestry_specific_results[[ancestry]][[model]]
    ci <- confint.default(r)

    biobanks <- paste0(r$`s.levels`[[1]], collapse=',')
    if (biobanks == ''){
      biobanks <- unique(r$data$bbid)
      tau_ci <- c(NA,NA,NA)
    } else {
      tau_ci <- confint(r)[[1]][2,]
    }

    return(
      data.frame(
        study=x$study,
        phenotype=x$phenotype,
        ancestry = ancestry,
        method = row.names(r$beta),
        beta = r$beta[,1],
        se = r$se,
        ci_low = ci[,1],
        ci_high = ci[,2],
        tau_biobank = tau_ci[1],
        tau_biobank_ci_low = tau_ci[2],
        tau_biobank_ci_high = tau_ci[3],
        tstat = r$zval,
        ddf = r$ddf,
        dfs = r$dfs,
        pval = r$pval,
        zstat = r$beta[,1] / r$se,
        pval_z = 2*pnorm(abs(r$beta[,1] / r$se), lower.tail = F),
        QE = r$QE,
        QEp = r$QEp,
        biobanks = biobanks
      )
    )
  }
}

extract_methods_results_pv <- function(x, ancestry, model='m_mv_method'){

  if (length(x$ancestry_specific_results[[ancestry]])>0){

    if (model == 'm_mv_method'){
      r <- x$ancestry_specific_results[[ancestry]]$pv_method
      m <- x$ancestry_specific_results[[ancestry]]$m_mv_method
      biobanks <- paste0(m$`s.levels`[[1]], collapse=',')
      if (biobanks == ''){
        biobanks <- unique(m$data$bbid)
      }

    } else {
      stop('invalid model specified.')
    }

    biobanks <- paste0(m$`s.levels`[[1]], collapse=',')
    if (biobanks == ''){
      biobanks <- unique(m$data$bbid)
    }

    r$study <- x$study
    r$phenotype <- x$phenotype
    r$ancestry <- ancestry
    r$biobanks <- biobanks

    return(r)
  }
}

results_eur <- rbindlist(lapply(results_list, extract_methods_results, ancestry='EUR', model='m_mv_method'))
results_sas <- rbindlist(lapply(results_list, extract_methods_results, ancestry='SAS', model='m_mv_method'))

pv_eur <- rbindlist(lapply(results_list, extract_methods_results_pv, ancestry='EUR', model='m_mv_method'))
pv_sas <- rbindlist(lapply(results_list, extract_methods_results_pv, ancestry='SAS', model='m_mv_method'))

meta_res <- rbindlist(list(results_eur, results_sas))
meta_res[,tuning_type:=str_split(method, '_', simplify = T)[,2]]
meta_res[,method:=gsub('^method','',str_split(method, '_', simplify = T)[,1])]

meta_res <- meta_res[order(study,phenotype,ancestry,-beta)]
meta_res <- rename_phenotypes(meta_res)

meta_res$method <- adjust_methods_levels(meta_res$method)

# export processed data --------------------------------------------------------
# let R pick the best compression scheme
# tools::resaveRdaFiles("data/", compress="auto")
# tools::checkRdaFiles("data/") -> gzip for this dataset
usethis::use_data(meta_res, overwrite = TRUE, version=2, compress="gzip")


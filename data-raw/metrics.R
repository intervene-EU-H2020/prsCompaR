library(stringr)
library(data.table)

source("data-raw/utils.R")

# functions --------------------------------------------------------------------

load_results_metadata <- function(basepath, biobanks=c('ukbb','ebb','finngen','GNH','hunt'), ancestries=c('AFR', 'AMR', 'EAS', 'EUR', 'SAS')){

  studies_path <- paste0(basepath,'config/studies_for_methods_comparison.tsv')
  studies <- fread(studies_path, sep='\t')

  phenos <- strsplit(studies$name,',')

  superpop <- ancestries

  study_ids <- rep(studies$study_id, times=sapply(phenos, length))
  phenos <- unlist(phenos)
  infiles <- list()

  # collecting all the files containing performance metrics and score-score correlations

  for (b in biobanks){
    for (s in superpop){
      for (i in seq_along(phenos)){

        # original workflow GenoPred metrics
        assoc_filepath <- paste0(basepath, 'results/',b,'/PRS_evaluation/',study_ids[i],'/',s,'/',study_ids[i],'.',phenos[i],'.',s,'.AllMethodComp.assoc.txt')
        bm_filepath <- paste0(basepath, 'results/',b,'/PRS_evaluation/',study_ids[i],'/',s,'/',study_ids[i],'.',phenos[i],'.',s,'.AllMethodComp.best_models.tsv')
        pe_filepath <- paste0(basepath, 'results/',b,'/PRS_evaluation/',study_ids[i],'/',s,'/',study_ids[i],'.',phenos[i],'.',s,'.AllMethodComp.pred_eval.txt')

        b_lower <- tolower(b)

        # dowstream workflow
        # re-calculated metrics
        mean_sd_full_filepath <- paste0(basepath,'01_prspipe_followup_results/',b_lower,'/downstream_workflows/results/metrics_and_scor_full/',study_ids[i],'/',study_ids[i],'.',phenos[i],'.',s,'.mean_sd.tsv')
        metrics_full_filepath <- paste0(basepath,'01_prspipe_followup_results/',b_lower,'/downstream_workflows/results/metrics_and_scor_full/',study_ids[i],'/',study_ids[i],'.',phenos[i],'.',s,'.metrics.tsv')
        mean_sd_traintest_filepath <- paste0(basepath,'01_prspipe_followup_results/',b_lower,'/downstream_workflows/results/metrics_and_scor_train_test/',study_ids[i],'/',study_ids[i],'.',phenos[i],'.',s,'.mean_sd.tsv')
        metrics_traintest_filepath <- paste0(basepath,'01_prspipe_followup_results/',b_lower,'/downstream_workflows/results/metrics_and_scor_train_test/',study_ids[i],'/',study_ids[i],'.',phenos[i],'.',s,'.metrics.tsv')
        score_score_correl_filepath <- paste0(basepath,'01_prspipe_followup_results/',b_lower,'/downstream_workflows/results/metrics_and_scor_train_test/',study_ids[i],'/',study_ids[i],'.',phenos[i],'.',s,'.scor.tsv')


        # nested pTclump, GenoPred metrics
        ptclump_nested_pe_filepath <- paste0(basepath,'01_prspipe_followup_results/',b_lower,'/downstream_workflows/results/PRS_evaluation/pt_clump_nested/',study_ids[i],'/',s,'/',study_ids[i],'.',phenos[i],'.',s,'.pred_eval.txt')
        ptclump_nested_assoc_filepath <- paste0(basepath,'01_prspipe_followup_results/',b_lower,'/downstream_workflows/results/PRS_evaluation/pt_clump_nested/',study_ids[i],'/',s,'/',study_ids[i],'.',phenos[i],'.',s,'.assoc.txt')

        if (file.exists(bm_filepath)){
          wc <- length(readLines(bm_filepath))
          if (wc > 1){
            infiles[[bm_filepath]] <- data.frame(bbid=tolower(b),
                                                 ancestry=s,
                                                 pheno=phenos[i],
                                                 study=study_ids[i],
                                                 assoc_path=assoc_filepath,
                                                 bm_path=bm_filepath,
                                                 pe_path=pe_filepath,
                                                 mean_sd_full_path=mean_sd_full_filepath,
                                                 metrics_full_path=metrics_full_filepath,
                                                 mean_sd_traintest_path=mean_sd_traintest_filepath,
                                                 metrics_traintest_path=metrics_traintest_filepath,
                                                 score_score_correl_path=score_score_correl_filepath,
                                                 ptclump_nested_pe_path=ptclump_nested_pe_filepath,
                                                 ptclump_nested_assoc_path=ptclump_nested_assoc_filepath
            )
          } else {
            print(paste(bm_filepath,'is empty. skipping.'))
          }
        }
      }
    }
  }

  file_metadata <- rbindlist(infiles)

  return(file_metadata)

}

get_tags <- function(s,p,a){
  return(tags_select[[paste(s,p,a,sep=':')]])
}

read_results <- function(r){

  # pred_eval.txt file loading --------------------------------------------
  # pred_eval files to not contain the UKBB MultiPRS

  pe <- fread(r$pe_path)
  pe$tag <- pe$Model
  for (m in methods){
    method_group <- paste0(m,'_group')
    pe[tag == method_group, tag:=paste0(r$study,'_',r$pheno,'_',m,'_',r$ancestry,'.MultiPRS')]
  }
  pe$tag <- gsub('_group$','',pe$tag)
  pe$tag <- gsub('.+\\.PredFile[0-9]+\\.','',pe$tag)
  pe$method <- str_extract(pe$Model, pattern = paste(methods, collapse='|'))
  pe[!grepl('MultiPRS',tag),tag:=paste(tag,method,sep='_')]
  pe$tag <- gsub('[[:punct:]]','.',pe$tag)
  pe$method <- NULL

  # pred_eval.txt file loading (nested ptclump) -----------------------------
  # these pred_eval files only contain pt+clump results
  # the nested ptclump MultiPRS is not evaluated in any of the other files

  pe_ptcn <- fread(r$ptclump_nested_pe_path)
  pe_ptcn$tag <- pe_ptcn$Model

  pe_ptcn[Model == 'All_group', tag:=paste0(r$study,'_',r$pheno,'_pt.clump.nested_',r$ancestry,'.MultiPRS')]
  pe_ptcn$tag <- gsub('_group$','',pe_ptcn$tag)
  pe_ptcn[!grepl('MultiPRS',tag),tag:=paste(tag,'pt.clump.nested',sep='_')]
  pe_ptcn$tag <- gsub('[[:punct:]]','.',pe_ptcn$tag)

  # assoc.txt file loading --------------------------------------------------
  # assoc files do not contain MultiPRS or UKBB MultiPRS

  assoc <- fread(r$assoc_path)
  assoc$tag <- assoc$Predictor
  assoc$method <- str_extract(assoc$tag, pattern = paste(methods, collapse='|'))
  assoc$tag <- gsub('^Group_.*PredFile[0-9]+\\.','',assoc$tag)
  assoc[,tag:=paste(tag,method,sep='_')]
  assoc$tag <- gsub('[[:punct:]]','.',assoc$tag)
  assoc$method <- NULL

  # assoc.txt file loading (nested ptclump) ----------------------------------
  # these pred_eval files only contain pt+clump results

  assoc_ptcn <- fread(r$ptclump_nested_assoc_path)
  assoc_ptcn$tag <- assoc_ptcn$Predictor
  assoc_ptcn$tag <- gsub('^Group_.*nested\\.','',assoc_ptcn$tag)
  assoc_ptcn[,tag:=paste(tag,'pt.clump.nested',sep='_')]
  assoc_ptcn$tag <- gsub('[[:punct:]]','.',assoc_ptcn$tag)

  # best_models.tsv file loading -----------------------------------------------
  # best_models files do not contain UKBB MultiPRS

  bm <- fread(r$bm_path)
  for (m in methods){
    method_group <- paste0(m,'_group')
    bm[tag == method_group, tag:=paste0(r$study,'_',r$pheno,'_',m,'_',r$ancestry,'.MultiPRS')]
  }
  bm$tag <- gsub('_group$','',bm$tag)
  bm$tag <- gsub('.+\\.PredFile[0-9]+\\.','',bm$tag)
  bm[!grepl('MultiPRS',tag),tag:=paste(tag,Method,sep='_')]
  bm$tag <- gsub('[[:punct:]]','.',bm$tag)

  # best models for nested ptclump --------------------------------------------

  bm_ptcn <- data.table::copy(pe_ptcn)
  bm_ptcn$group <- 'pt.clump.nested'
  bm_ptcn$method <- 'pt.clump.nested'
  bm_ptcn$is_multi <- ifelse(bm_ptcn$Model == 'All_group', 'yes', 'no')
  bm_ptcn$method_type <- ifelse(bm_ptcn$Model == 'All_group', 'MultiPRS', 'CV')
  bm_ptcn$label <- ifelse(bm_ptcn$Model == 'All_group', 'pT+clump.nested.MultiPRS', 'pT+clump.nested.CV')
  bm_ptcn$Model <- ifelse(bm_ptcn$Model == 'All_group', 'MultiPRS', 'CV')
  bm_ptcn$Method <- 'pt.clump.nested'
  bm_ptcn[,select:=CrossVal_R==max(CrossVal_R),by=method_type]
  bm_ptcn <- bm_ptcn[select == TRUE]
  bm_ptcn$select <- NULL


  # mean_sd.tsv file loading --------------------------------------------------
  # mean-sd (full)

  m_sd <- fread(r$mean_sd_full_path)
  m_sd$tag <- gsub('[[:punct:]]','.',m_sd$variable)
  setnames(m_sd, gsub('_train','_all',colnames(m_sd)))

  # metrics.tsv file loading --------------------------------------------------
  # re-calculated metrics (full)

  met <- fread(r$metrics_full_path)
  met$tag <- gsub('[[:punct:]]','.',met$predictor)

  # mean_sd.tsv file loading --------------------------------------------------
  # mean-sd (train-test)

  m_sd_tt <- fread(r$mean_sd_traintest_path)
  m_sd_tt$tag <- gsub('[[:punct:]]','.',m_sd_tt$variable)

  # metrics.tsv file loading --------------------------------------------------
  # re-calculated metrics (train-test)

  met_tt <- fread(r$metrics_traintest_path)
  met_tt$tag <- gsub('[[:punct:]]','.',met_tt$predictor)

  # scor.tsv file loading -----------------------------------------------------
  # score score correlations

  scor <- fread(r$score_score_correl_path)
  tag <- scor$V1
  scor$V1 <- NULL
  tag <- gsub('[[:punct:]]','.',tag)
  scor <- as.matrix(scor)
  colnames(scor) <- tag
  row.names(scor) <- tag
  rm(tag)


  ## merging -------------------------------------------------------------------

  assoc <- rbindlist(list(assoc, assoc_ptcn), use.names=T, fill=T)
  pe <- rbindlist(list(pe, pe_ptcn), use.names=T, fill=T)
  bm <- rbindlist(list(bm, bm_ptcn), use.names=T, fill=T)

  m_sd <- merge(m_sd, m_sd_tt, by = c('variable','tag'))
  met <- rbindlist(list(met, met_tt))



  return(list(pred_eval=pe, assoc=assoc, bm=bm, m_sd=m_sd, met=met, scor=scor, bbid=r$bbid, ancestry=r$ancestry, study=r$study, phenotype=r$pheno))

}


adjust_methods_levels <- function(x){
  x <- gsub('pt\\.clump\\.nested', 'pt.clump', x)
  x[grepl('UKBB-[A-Z]+\\.MultiPRS',x)] <- 'UKBB.EnsPRS'
  x <- factor(x, levels=methods_levels)
  if(any(is.na(x))){
    stop()
  }
  return(x)
}

adjust_methods_levels_2 <- function(x){
  x <- gsub('pt\\.clump\\.nested', 'pt.clump', x)
  x[grepl('UKBB-[A-Z]+\\.MultiPRS',x)] <- 'UKBB.EnsPRS_CV'
  x <- factor(x, levels=methods_type_levels)
  if(any(is.na(x))){
    stop('adjusting methods levels introduced NA')
  }
  return(x)
}

get_single_results <- function(i){

  r <- read_results(results_metadata[i])

  # get metrics calculated on the full sample, or just the test set in case of ukbb

  tags <- get_tags(r$study, r$pheno, 'EUR') # always get the EUR tags
  method_type <- tags$method_type
  tags <- tags$tag

  if (r$bbid == 'ukbb'){
    if (r$ancestry == 'EUR'){
      performance <- r$met[tag %in% tags & split == 'test']
    } else {
      performance <- r$met[tag %in% tags & split == 'all']
    }
  } else {
    performance <- r$met[tag %in% tags & split == 'all']
  }

  if (any(duplicated(performance$tag))){
    stop('duplicated model tags!')
  }

  performance <- performance[match(tags, tag)]
  performance$method <- sapply(strsplit(performance$predictor,'_'), function(x){x[length(x)]})
  performance$method_type <- method_type

  # get score-score correlations
  score_cor <- r$scor[tags, tags]

  if (r$phenotype == 'Chronic_kidney_disease'){
    # the PRS for these are flipped, because eGFR is negatively related to CKD
    # the ensemble PRS is not affected
    performance[!grepl('MultiPRS', method),c('BETA','BETA_CI_HIGH','BETA_CI_LOW'):=list(-BETA, -BETA_CI_LOW, -BETA_CI_HIGH)]
    performance[!grepl('MultiPRS', method),c('OR','OR_CI_HIGH','OR_CI_LOW'):=list(1/OR, 1/OR_CI_LOW, 1/OR_CI_HIGH)]
    performance[!grepl('MultiPRS', method),c('R2_OBS_CI_BOOT_HIGH','R2_OBS_CI_BOOT_LOW'):=list(R2_OBS_CI_BOOT_LOW, R2_OBS_CI_BOOT_HIGH)]
    if (!exclude_multiprs){
      # correlations also have to be flipped
      score_cor[grepl('MultiPRS', performance$method),] <- -1 * score_cor[grepl('MultiPRS', performance$method),]
      score_cor[,grepl('MultiPRS', performance$method)] <- -1 * score_cor[,grepl('MultiPRS', performance$method)]
    }
  }

  # calculate r2 on the liability scale
  if (r$ancestry == 'EUR'){
    k <- prevalence$V2[match(r$pheno, prevalence$V1)]
  } else if (r$ancestry == 'SAS') {
    k <- prevalence_sas$V2[match(r$pheno, prevalence_sas$V1)]
  } else {
    # use the sample prevalence
    k <- performance$N_CAS / performance$N
    stopifnot(all(k == k[1]))
    k <- k[1]
  }
  # this option is set in ../results_loading.R
  if (use_sample_prevalence){
    k <- performance$N_CAS / performance$N
    stopifnot(all(k == k[1]))
    k <- k[1]
  }

  p <- performance$N_CAS / performance$N
  stopifnot(all(p == p[1]))
  performance[,h2l:=h2l_R2(k, R2_OBS, p)]
  performance[,h2l_ci_low:=h2l_R2(k, R2_OBS_CI_BOOT_LOW, p)]
  performance[,h2l_ci_high:=h2l_R2(k, R2_OBS_CI_BOOT_HIGH, p)]
  performance[,h2l_k:=k]
  performance[,h2l_p:=p]

  # calculate p-values for the differences in BETA, adjusting for the correlations between scores
  dst <- data.table(as.matrix(dist(performance$BETA, method = 'manhattan')))
  dst <- melt(dst, variable.name = 'iy', measure.vars = colnames(dst), value.name='beta_diff')
  dst$iy <- as.numeric(as.character(dst$iy))
  dst$ix <- rep(1:nrow(performance), nrow(performance))
  dst$se_x <- performance$SD[dst$ix]
  dst$beta_x <- performance$BETA[dst$ix]
  dst$se_y <- performance$SD[dst$iy]
  dst$beta_y <- performance$BETA[dst$iy]
  dst$beta_diff <- with(dst, ifelse(beta_x > beta_y, beta_diff, -beta_diff))
  dst$cor_x_y <- melt(data.table(score_cor), measure.vars = 1:ncol(score_cor))$value

  dst[,se_diff:=sqrt(se_x^2 + se_y^2 - 2*se_x*se_y*cor_x_y)]
  dst[,z_diff:=beta_diff/se_diff]
  dst[,p_z_greater:=pnorm(-z_diff)]
  dst[,p_z_smaller:=pnorm(z_diff)]
  dst[,p_z_equal:=2*pnorm(abs(z_diff), lower.tail = F)]

  # carry over the sample size
  dst$n_x <- performance$N[dst$ix]
  dst$n_y <- performance$N[dst$iy]
  suppressWarnings({
    dst$ncas_x <- performance$N_CAS[dst$ix]
    dst$ncas_y <- performance$N_CAS[dst$iy]
    dst$ncon_x <- performance$N_CON[dst$ix]
    dst$ncon_y <- performance$N_CON[dst$iy]
  })

  # carry over the other metrics for relative comparisons
  suppressWarnings({
    dst$auc_x <- performance$AUC_MEDIAN[dst$ix]
    dst$auc_y <- performance$AUC_MEDIAN[dst$iy]
    dst$h2l_x <- performance$h2l[dst$ix]
    dst$h2l_y <- performance$h2l[dst$iy]
    dst$r2_x <- performance$R2_OBS[dst$ix]
    dst$r2_y <- performance$R2_OBS[dst$iy]
  })

  # add information about the method
  dst[,method_x:=performance$method[ix]]
  dst[,method_type_x:=performance$method_type[ix]]
  dst[,method_y:=performance$method[iy]]
  dst[,method_type_y:=performance$method_type[iy]]

  # carry over the model tags
  dst$tag_x <- performance$tag[dst$ix]
  dst$tag_y <- performance$tag[dst$iy]


  if (!(all(row.names(score_cor) == performance$tag))){
    stop('failed sanity check, correlation matrix mismatch')
  }
  if (!(all(colnames(score_cor) == performance$tag))){
    stop('failed sanity check, correlation matrix mismatch')
  }

  result <- list(metrics=performance, scor=score_cor, pairwise_tests=dst, study=r$study, phenotype=r$phenotype, ancestry=r$ancestry, bbid=r$bbid)

  return(result)

}

get_block_correl_matrix <- function(r){
  # used in the function below, returns a block-diagonal correlation matrix between scores to adjust for correlated effect-size estimates
  return(as.matrix(bdiag(lapply(r, '[[','scor'))))
}

# function that fetches the tags of the scores we want to compare
get_auto_cv_and_multiPRS_model_tags <- function(i){

  r <- read_results(results_metadata[i])

  automatic_model_tags <- r$bm[method_type %in% c('PseudoVal')]
  automatic_model_tags <- automatic_model_tags[!duplicated(label)]$tag
  automatic_model_tags <- c(automatic_model_tags, paste0(r$study,'.0.1e.08.pt.clump'))
  automatic_model_tags <- data.table(tag=automatic_model_tags, method_type = 'auto')

  tuned_model_tags <- r$bm[(method_type == 'CV') | ((method_type == 'PseudoVal') & (method %in% c('dbslmm','sbayesr')))]
  tuned_model_tags <- tuned_model_tags[!grepl('pt.clump$',tag)]
  tuned_model_tags <- tuned_model_tags[!duplicated(label),list(tag, method_type)]

  if (!exclude_multiprs){
    multiprs_tags <- r$met[split=='all'][grepl('All.UKBB.EUR',tag), list(tag)]
    multiprs_tags$method_type <- 'CV'
    tuned_model_tags <- rbind(tuned_model_tags, multiprs_tags)
  }

  tags <- rbindlist(list(automatic_model_tags, tuned_model_tags))
  tags <- tags[!(duplicated(tag))] # drop the tuned model if it is equal to the automatic model

  if(any(duplicated(tags$tag))){
    stop('duplicated model tags!')
  }

  if (r$phenotype %in% auto_only_phenotypes){
    tags <- tags[method_type == 'auto']
  }

  return(tags)

}

get_tags <- function(s,p,a){
  return(tags_select[[paste(s,p,a,sep=':')]])
}

h2l_R2 <- function(k, r2, p) {

  # K baseline disease risk
  # r2 from a linear regression model attributable to genomic profile risk score
  # P proportion of sample that are cases
  # calculates proportion of variance explained on the liability scale
  # from ABC at http://www.complextraitgenomics.com/software/
  # Lee SH, Goddard ME, Wray NR, Visscher PM. (2012) A better coefficient of determination for genetic profile analysis. Genet Epidemiol. 2012 Apr;36(3):214-24.

  if (length(k)==0){
    return(NA)
  }

  if (is.null(k)){
    return(NA)
  }

  if (is.na(k)){
    return(NA)
  }

  x= qnorm(1-k)
  z= dnorm(x)
  i=z/k
  C= k*(1-k)*k*(1-k)/(z^2*p*(1-p))
  theta= i*((p-k)/(1-k))*(i*((p-k)/(1-k))-x)
  h2l_R2 = C*r2 / (1 + C*theta*r2)
  return(h2l_R2)
}

# parameters
# ------------------------------------------------------------------------------
methods_levels <- c('pt.clump','dbslmm','ldpred2','lassosum','megaprs','prscs','sbayesr','UKBB.EnsPRS')
exclude_multiprs <- FALSE
use_sample_prevalence <- FALSE
basepath <- file.path(getwd(), "data-raw/")
methods <-
  c(
    'All',
    'lassosum',
    'ldpred2',
    'megaprs',
    'pt.clump.nested',
    'pt.clump',
    'dbslmm',
    'sbayesr',
    'prscs'
  )

# download data
if (!file.exists("data-raw/results.tar.gz")) {
  url <- "https://figshare.com/ndownloader/files/40418924?private_link=c9ea515995d8f61c74a0"
  download.file(url, "data-raw/results.tar.gz")
  untar("data-raw/results.tar.gz", exdir = "data-raw/")
}

# load data
# ------------------------------------------------------------------------------
prevalence <- fread(paste0(basepath, '/config/pop_prevalence.tsv'), sep='\t')
prevalence_sas <- fread(paste0(basepath, '/config/pop_prevalence_sas.tsv'), sep='\t')
results_metadata <- load_results_metadata(basepath, ancestries=c("EUR", "SAS"))
file_cols <- grep('path', colnames(results_metadata), value=TRUE)

if (exclude_multiprs){
  methods_levels <- methods_levels[-length(methods_levels)]
}

# general exclusions
results_metadata <- results_metadata[study != 'GCST003401']
results_metadata <- results_metadata[!(study == 'GCST007954' & pheno == 'T2D')]
# auto-only phenotypes (no CV)
auto_only_phenotypes <- c('Alzheimers_disease','AD','Height')


i_ukbb <- which(results_metadata$bbid == 'ukbb' & results_metadata$ancestry == 'EUR')
tags_select <- lapply(i_ukbb, get_auto_cv_and_multiPRS_model_tags)
names(tags_select)<-paste(results_metadata[i_ukbb]$study,results_metadata[i_ukbb]$pheno,results_metadata[i_ukbb]$ancestry,sep=':')

loader_fun <- get_single_results
single_results <- lapply(1:nrow(results_metadata), loader_fun)

metrics <- rbindlist(lapply(single_results, function(x){
  # extract metrics
  p <- x$metrics
  p$phenotype <- x$phenotype
  p$study <- x$study
  p$ancestry <- x$ancestry
  p$bbid <- x$bbid
  return(p)
}), fill=TRUE, use.names=TRUE)

# exclusions:
# sample overlap
metrics <- metrics[!(bbid == 'ukbb' & phenotype == 'Height' & ancestry == 'EUR')]
metrics <- metrics[!(bbid == 'ebb' & phenotype == 'Creatinine_eGFR')]
metrics <- metrics[!(bbid == 'ukbb' & phenotype == 'Alzheimers_disease' & ancestry == 'EUR')]

metrics <- rename_phenotypes(metrics)
metrics$method <- adjust_methods_levels(metrics$method)

# label continuous / binary traits
metrics$is_continuous <- ifelse(is.na(metrics$AUC_MEDIAN), TRUE, FALSE)

# export processed data --------------------------------------------------------
# let R pick the best compression scheme
# tools::resaveRdaFiles("data/", compress="auto")
# tools::checkRdaFiles("data/") -> gzip for this dataset
usethis::use_data(metrics, overwrite = TRUE, version=2, compress="gzip")


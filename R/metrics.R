#' Polygenic risk score performance metrics table for single biobanks
#'
#' A subset of data from... (science stuff here)
#'
#' @format
#' A data frame with 950 rows and 49 columns:
#' \describe{
#'   \item{predictor}{unique score identifier}
#'   \item{BETA}{linear/logistic regression coef}
#'   \item{SD}{standard deviation of BETA}
#'   \item{T_STAT}{for continuous traits, the t-statistic used to produce PVAL}
#'   \item{PVAL}{p-value for BETA (H_0: BETA = 0)}
#'   \item{BETA_CI_LOW}{95\% confidence intervals for BETA, lower bound}
#'   \item{BETA_CI_HIGH}{95\% confidence intervals for BETA, upper bound}
#'
#'   \item{R2_OBS}{variance explained on the observed scale (only really meaningful for continuous traits)}
#'   \item{R2_OBS_CI_LOW}{95\% CI lower bound for R2_OBS (linear regression based)}
#'   \item{R2_OBS_CI_HIGH}{95\% CI upper bound for R2_OBS (linear regression based)}
#'
#'   \item{N}{total sample size (number of cases + number of controls)}
#'   \item{split}{the data split used to calculate the metrics (train/test/all)}
#'   \item{tag}{similar to predictor, but uglier}
#'   \item{method}{PRS method}
#'   \item{method_type}{\strong{tuning type} (auto/CV)}
#'
#'   \item{h2l}{variance explained on the liability scale (this should be used instead of R2 for binary traits, based on R2_OBS)}
#'   \item{h2l_ci_low}{95\% CI lower bound for the variance explained on the liability scale (based on R2_OBS_BOOT_CI_LOW)}
#'   \item{h2l_ci_high}{95\% CI upper bound for the variance explained on the liability scale (based on R2_OBS_BOOT_CI_HIGH)}
#'   \item{h2l_k}{sample prevalence used to calculate h2l}
#'   \item{h2l_p}{population prevalence used to calculate h2l}
#'
#'   \item{phenotype}{the evaluated endpoint}
#'   \item{study}{the GWAS study identifier (GWAS catalog)}
#'
#'   \item{ancestry}{the target ancestry}
#'   \item{bbid}{the target biobank}
#'
#'   \item{Z_VAL}{for binary endpoints, the test statistic used to produce PVAL}
#'
#'   \item{BETA_OBS}{linear regression coefficient on the observed (probability) scale (calculated for binary endpoints only)}
#'   \item{SD_OBS}{standard deviation of BETA_OBS}
#'   \item{T_VAL_OBS}{test statistic used to calculate PVAL_OBS}
#'   \item{PVAL_OBS}{p-value for BETA_OBS (H_0: BETA_OBS = 0)}
#'   \item{BETA_OBS_CI_LOW}{95\% confidence intervals for BETA_OBS, lower bound}
#'   \item{BETA_OBS_CI_HIGH}{95\% confidence intervals for BETA_OBS, upper bound}
#'
#'   \item{R_OBS_BOOT_T0}{correlation coefficient of the PRS with the observed probability-scale outcome (binary endpoints only) for the original data}
#'   \item{R_OBS_BOOT_MEDIAN}{Median correlation coefficient of the PRS with the observed probability-scale outcome (binary endpoints only) for 1000 bootstrap samples}
#'   \item{R_OBS_BOOT_MEAN}{Mean correlation coefficient of the PRS with the observed probability-scale outcome (binary endpoints only) for 1000 bootstrap samples}
#'   \item{R_OBS_BOOT_SD}{Standard deviation of correlation coefficient of the PRS with the observed probability-scale outcome (binary endpoints only) for 1000 bootstrap samples}
#'   \item{R2_OBS_BOOT_MEDIAN}{Median variance explained on the probability-scale outcome (binary endpoints only) for 1000 bootstrap samples}
#'   \item{R2_OBS_BOOT_MEAN}{Mean variance explained on the probability-scale outcome (binary endpoints only) for 1000 bootstrap samples}
#'
#'   \item{R_OBS_CI_BOOT_LOW}{95\% CI lower bound for R_OBS based on 1000 bootstrap samples}
#'   \item{R_OBS_CI_BOOT_HIGH}{95\% CI upper bound for R_OBS based on 1000 bootstrap samples}
#'   \item{R2_OBS_CI_BOOT_LOW}{95\% CI lower bound for R2_OBS based on 1000 bootstrap samples}
#'   \item{R2_OBS_CI_BOOT_HIGH}{95\% CI upper bound for R2_OBS based on 1000 bootstrap samples}
#'
#'   \item{AUC_CI_LOW}{95\% CI lower bound for the AUROC}
#'   \item{AUC_MEDIAN}{AUROC for binary endpoints}
#'   \item{AUC_CI_HIGH}{95\% CI upper bound for the AUROC}
#'
#'   \item{OR}{Odds ratio per standard deviation (based on BETA) for binary endpoints}
#'   \item{OR_CI_LOW}{95\% CI lower bound for OR}
#'   \item{OR_CI_HIGH}{95\% CI upper bound for OR}
#'
#'   \item{N_CAS}{number of cases for binary endpoints}
#'   \item{N_CON}{number of controls for binary endpoints}
#' }
#' @source <https://link/to/paper>
"metrics"

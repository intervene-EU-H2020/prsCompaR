# extract r2 fields correctly for the PGS Catalog metrics update (PGP000517)
# 2024-11-11

# devtools::install_github("intervene-EU-H2020/pgsCompaR")

library(tidyverse)
library(pgsCompaR)

data(metrics, overwrite = TRUE)

# 945 in submission template - 5 rows are dropped from submission data?
dim(metrics)

# set R2 to be a consistent measurement
r2 <- metrics %>%
  mutate(consistent_R2 = ifelse(is_continuous, R2_OBS, h2l),
         consistent_R2_lower = ifelse(is_continuous, R2_OBS_CI_LOW, h2l_ci_low),
         consistent_R2_upper = ifelse(is_continuous, R2_OBS_CI_HIGH, h2l_ci_high))

# 1. recode metrics$method. data(metrics) -> submission template:
# pt.clump -> pt_clump
# pt.clump -> pt_clump_nested (need to extract)
# UKBB.EnsPRS -> UKBB-EUR.MultiPRS
# (other stay the same)

r2$methods_harmonised <- str_replace_all(sapply(metrics$predictor, function(x) tail(strsplit(x, "_")[[1]], 1)), fixed("."), "_")
r2$methods_harmonised <- fct_recode(r2$methods_harmonised, `UKBB-EUR.MultiPRS` = "UKBB-EUR_MultiPRS")
table(r2$methods_harmonised)

# 2. recode metrics$phenotype. data(metrics) -> submission template:
# Prostate cancer -> Prostate_cancer
# Breast cancer -> Breast_cancer
# (everything else is the same)
r2$phenotype <- str_replace(metrics$phenotype, fixed(" "), "_")
table(r2$phenotype)

# 3. create score name (on PGS Catalog submission template)
r2$score_name <-stringr::str_glue_data(r2, "INTERVENE_{methods_harmonised}.{method_type}_{phenotype}")

# 4. create sample set ID (on PGS Catalog submission template)
r2$bbid <- fct_recode(metrics$bbid, UKB = "ukbb", EB = "ebb", `G&H` = "gnh", HUNT = "hunt", FinnGen = "finngen")
r2$sample_set_id <- stringr::str_glue_data(r2, "{bbid}_{ancestry}_{phenotype}")

r2 %>%
  select(score_name, sample_set_id, consistent_R2, consistent_R2_lower, consistent_R2_upper) %>%
  write_tsv("r2.tsv")

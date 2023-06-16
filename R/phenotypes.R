get_auto_only_phenotypes <- function() {
  return(c('Alzheimers_disease','AD','Height'))
}

rename_phenotypes <- function(x){
  plotdata <- data.table::copy(x)
  if ('phenotype' %in% colnames(plotdata)){
    plotdata[phenotype=='Chronic_kidney_disease', phenotype:='CKD']
    plotdata[phenotype=='Alzheimers_disease', phenotype:='AD']
    plotdata[phenotype=='Urate_Gout', phenotype:='Gout']
    plotdata[phenotype=='Stroke_excluding_SAH', phenotype:='Stroke']
    plotdata[phenotype=='Inflammatory_bowel_disease', phenotype:='IBD']
    plotdata[phenotype=='Seropositive_rheumatoid_arthritis', phenotype:='Arthritis']
    plotdata[phenotype=='Breast_cancer', phenotype:='Breast cancer']
    plotdata[phenotype=='Prostate_cancer', phenotype:='Prostate cancer']
    plotdata[phenotype=='Seropositive_rheumatoid_arthritis', phenotype:='Arthritis']
    plotdata[phenotype=='HDL_cholesterol', phenotype:='HDL']
    plotdata[phenotype=='Creatinine_eGFR', phenotype:='eGFR']
  }
  if ('pheno' %in% colnames(plotdata)){
    plotdata[,phenotype:=pheno]
    plotdata[pheno=='Chronic_kidney_disease', phenotype:='CKD']
    plotdata[pheno=='Alzheimers_disease', phenotype:='AD']
    plotdata[pheno=='Urate_Gout', phenotype:='Gout']
    plotdata[pheno=='Stroke_excluding_SAH', phenotype:='Stroke']
    plotdata[pheno=='Inflammatory_bowel_disease', phenotype:='IBD']
    plotdata[pheno=='Seropositive_rheumatoid_arthritis', phenotype:='Arthritis']
    plotdata[pheno=='Breast_cancer', phenotype:='Breast cancer']
    plotdata[pheno=='Prostate_cancer', phenotype:='Prostate cancer']
    plotdata[pheno=='Seropositive_rheumatoid_arthritis', phenotype:='Arthritis']
    plotdata[pheno=='HDL_cholesterol', phenotype:='HDL']
    plotdata[pheno=='Creatinine_eGFR', phenotype:='eGFR']
  }
  return(plotdata)
}

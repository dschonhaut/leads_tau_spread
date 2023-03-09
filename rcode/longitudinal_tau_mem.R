library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(tidyverse)

# --------------------------
# Import and format the data
data_dir = "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
filepath = file.path(data_dir, "tau_rois.csv")
tau_rois = read.csv(filepath)

# Replace missing values with NA
cols_to_replace = c("minority", "apoe_genotype", "apoe4_alleles", "cdrsob_baseline", "mmse_baseline")
for (col in cols_to_replace) {
  tau_rois[,col] = replace(tau_rois[,col], tau_rois[,col] == -999, NA)
  missing_values = sum(is.na(tau_rois[,col]))
  cat(col, ": ", missing_values, " missing values\n")
}
for (col in cols_to_replace) {
  na_count = sum(is.na(tau_rois[,col]))
  non_na_count = sum(!is.na(tau_rois[,col]))
  cat(col, ": ", na_count, " NA values and ", non_na_count, " non-NA values\n")
}

# Subset the dataframe
tau_rois$subjhemroi = paste(tau_rois$subj, tau_rois$hemroi, sep="_")
tau_rois$cdrsob_gt4 = tau_rois$cdrsob_baseline > 4
tau_rois_bline = tau_rois[(tau_rois$visit==1),]
tau_rois_bline_eoad_con = tau_rois_bline[(tau_rois_bline$dx %in% c("EOAD", "CON")),]
tau_rois_bline_eoad_eonad = tau_rois_bline[(tau_rois_bline$dx %in% c("EOAD", "EOnonAD")),]
tau_rois_bline_eonad_con = tau_rois_bline[(tau_rois_bline$dx %in% c("EOnonAD", "CON")),]
tau_rois_bline_eoad = tau_rois_bline[(tau_rois_bline$dx=="EOAD"),]
tau_rois_eoad = tau_rois[(tau_rois$dx=="EOAD") & (tau_rois$pet_visits>1),]
tau_rois_eoad_all = tau_rois[(tau_rois$dx=="EOAD"),]

# Convert strings to factors
tau_rois_bline = as.data.frame(unclass(tau_rois_bline), stringsAsFactors=TRUE)
tau_rois_bline_eoad_con = as.data.frame(unclass(tau_rois_bline_eoad_con), stringsAsFactors=TRUE)
tau_rois_bline_eoad_eonad = as.data.frame(unclass(tau_rois_bline_eoad_eonad), stringsAsFactors=TRUE)
tau_rois_bline_eoad_eonad$dx = ordered(tau_rois_bline_eoad_eonad$dx, c("EOnonAD", "EOAD"))
tau_rois_bline_eonad_con = as.data.frame(unclass(tau_rois_bline_eonad_con), stringsAsFactors=TRUE)
tau_rois_bline_eoad = as.data.frame(unclass(tau_rois_bline_eoad), stringsAsFactors=TRUE)
tau_rois_eoad = as.data.frame(unclass(tau_rois_eoad), stringsAsFactors=TRUE)
tau_rois_eoad_all = as.data.frame(unclass(tau_rois_eoad_all), stringsAsFactors=TRUE)

eoad_subjs = levels(as.factor(tau_rois_eoad$subj))
hemrois = levels(as.factor(tau_rois_eoad$hemroi))
subjhemrois = levels(as.factor(tau_rois_eoad$subjhemroi))

# -------------------------------------------------------
# LMs for baseline comparisons between diagnostic groups

## ------------
## EOAD vs. CON
mod_icpt   = lm(value ~ 1,
                data=tau_rois_bline_eoad_con)
mod_dx     = lm(value ~ 1 + dx,
                data=tau_rois_bline_eoad_con)
mod_dx_roi = lm(value ~ 1 + dx + hemroi,
                  data=tau_rois_bline_eoad_con)
mod_dx.roi = lm(value ~ 1 + dx + hemroi + dx:hemroi,
                  data=tau_rois_bline_eoad_con)

# Test: Tau PET SUVR values vary between diagnostic groups
anova(mod_dx, mod_icpt)
summary(mod_dx)[["coefficients"]]

# Test: Tau PET values vary between diagnostic groups, differentially across regions.
anova(mod_dx.roi, mod_dx_roi)

# Run a separate model in each region and extract the coefficients.
coef_df = data.frame()
for (hemroi in hemrois) {
  dat = tau_rois_bline_eoad_con[(tau_rois_bline_eoad_con$hemroi==hemroi),]
  mod = lm(value ~ 1 + dx, data=dat)
  ci = confint(mod, method="boot") # 95% confidence interval by parametric bootstrapping
  colnames(ci) = c("ci2.5", "ci97.5")
  coef_vals = summary(mod)[["coefficients"]]
  coef_names = rownames(coef_vals)
  for (coef_name in coef_names[2:length(coef_names)]) {
    coef_row = data.frame(hemroi = hemroi, nobs = nrow(dat), coef = coef_name, t(coef_vals[coef_name,]), t(ci[coef_name,]))
    coef_df = rbind(coef_df, coef_row) # append to the dataframe
  }
}
names(coef_df)[names(coef_df)=="Estimate"] = "est"
names(coef_df)[names(coef_df)=="Std..Error"] = "sem"
names(coef_df)[names(coef_df)=="t.value"] = "tval"
names(coef_df)[names(coef_df)=="Pr...t.."] = "pval"
coef_df = relocate(coef_df, any_of(c("ci2.5", "ci97.5")), .after = 5)

# FDR-correct P-values and determine significance.
alpha = 0.05
coef_df["pval_fdr"] = p.adjust(coef_df[,"pval"], method = "BH")
coef_df["sig_fdr"] = coef_df["pval_fdr"] < alpha

coef_eoad_con_bline = coef_df
View(coef_eoad_con_bline)
if (TRUE) {
  filepath = file.path(data_dir, "coef_eoad_con_bline.csv")
  write.csv(coef_eoad_con_bline, filepath, row.names=FALSE)
  cat("Saved", filepath)
}

# Checks
mean(tau_rois[(tau_rois$dx=="EOAD") & (tau_rois$visit==1) & (tau_rois$hemroi=="L_entorhinal"), "value"])
mean(tau_rois[(tau_rois$dx=="EOAD") & (tau_rois$visit==1) & (tau_rois$hemroi=="R_amygdala"), "value"])
mean(tau_rois[(tau_rois$dx=="EOAD") & (tau_rois$visit==1) & (tau_rois$hemroi=="R_hippocampus"), "value"])

## ----------------
## EOAD vs. EOnonAD
mod_icpt   = lm(value ~ 1,
                data=tau_rois_bline_eoad_eonad)
mod_dx     = lm(value ~ 1 + dx,
                data=tau_rois_bline_eoad_eonad)
mod_dx_roi = lm(value ~ 1 + dx + hemroi,
                data=tau_rois_bline_eoad_eonad)
mod_dx.roi = lm(value ~ 1 + dx + hemroi + dx:hemroi,
                data=tau_rois_bline_eoad_eonad)

# Test: Tau PET SUVR values vary between diagnostic groups
anova(mod_dx, mod_icpt)
summary(mod_dx)[["coefficients"]]

# Test: Tau PET values vary between diagnostic groups, differentially across regions.
anova(mod_dx.roi, mod_dx_roi)

# Run a separate model in each region and extract the coefficients.
coef_df = data.frame()
for (hemroi in hemrois) {
  dat = tau_rois_bline_eoad_eonad[(tau_rois_bline_eoad_eonad$hemroi==hemroi),]
  mod = lm(value ~ 1 + dx, data=dat)
  ci = confint(mod, method="boot") # 95% confidence interval by parametric bootstrapping
  colnames(ci) = c("ci2.5", "ci97.5")
  coef_vals = summary(mod)[["coefficients"]]
  coef_names = rownames(coef_vals)
  for (coef_name in coef_names[2:length(coef_names)]) {
    coef_row = data.frame(hemroi = hemroi, nobs = nrow(dat), coef = coef_name, t(coef_vals[coef_name,]), t(ci[coef_name,]))
    coef_df = rbind(coef_df, coef_row) # append to the dataframe
  }
}
names(coef_df)[names(coef_df)=="Estimate"] = "est"
names(coef_df)[names(coef_df)=="Std..Error"] = "sem"
names(coef_df)[names(coef_df)=="t.value"] = "tval"
names(coef_df)[names(coef_df)=="Pr...t.."] = "pval"
coef_df = relocate(coef_df, any_of(c("ci2.5", "ci97.5")), .after = 5)

# FDR-correct P-values and determine significance.
alpha = 0.05
coef_df["pval_fdr"] = p.adjust(coef_df[,"pval"], method = "BH")
coef_df["sig_fdr"] = coef_df["pval_fdr"] < alpha

coef_eoad_eonad_bline = coef_df
View(coef_eoad_eonad_bline)

## ---------------
## EOnonAD vs. CON
mod_icpt   = lm(value ~ 1,
                data=tau_rois_bline_eonad_con)
mod_dx     = lm(value ~ 1 + dx,
                data=tau_rois_bline_eonad_con)
mod_dx_roi = lm(value ~ 1 + dx + hemroi,
                data=tau_rois_bline_eonad_con)
mod_dx.roi = lm(value ~ 1 + dx + hemroi + dx:hemroi,
                data=tau_rois_bline_eonad_con)

# Test: Tau PET SUVR values vary between diagnostic groups
anova(mod_dx, mod_icpt)
summary(mod_dx)[["coefficients"]]

# Test: Tau PET values vary between diagnostic groups, differentially across regions.
anova(mod_dx.roi, mod_dx_roi)

# Run a separate model in each region and extract the coefficients.
coef_df = data.frame()
for (hemroi in hemrois) {
  dat = tau_rois_bline_eonad_con[(tau_rois_bline_eonad_con$hemroi==hemroi),]
  mod = lm(value ~ 1 + dx, data=dat)
  ci = confint(mod, method="boot") # 95% confidence interval by parametric bootstrapping
  colnames(ci) = c("ci2.5", "ci97.5")
  coef_vals = summary(mod)[["coefficients"]]
  coef_names = rownames(coef_vals)
  for (coef_name in coef_names[2:length(coef_names)]) {
    coef_row = data.frame(hemroi = hemroi, nobs = nrow(dat), coef = coef_name, t(coef_vals[coef_name,]), t(ci[coef_name,]))
    coef_df = rbind(coef_df, coef_row) # append to the dataframe
  }
}
names(coef_df)[names(coef_df)=="Estimate"] = "est"
names(coef_df)[names(coef_df)=="Std..Error"] = "sem"
names(coef_df)[names(coef_df)=="t.value"] = "tval"
names(coef_df)[names(coef_df)=="Pr...t.."] = "pval"
coef_df = relocate(coef_df, any_of(c("ci2.5", "ci97.5")), .after = 5)

# FDR-correct P-values and determine significance.
alpha = 0.05
coef_df["pval_fdr"] = p.adjust(coef_df[,"pval"], method = "BH")
coef_df["sig_fdr"] = coef_df["pval_fdr"] < alpha

coef_eonad_con_bline = coef_df
View(coef_eonad_con_bline)

# --------------------------------------
# MEMs for longitudinal analysis in EOAD
# mod_roi      = lmer(value ~ -1 + hemroi                                                            + (1|subjhemroi),
#                     data=tau_rois_eoad, REML=FALSE)
# mod_roi_time = lmer(value ~ -1 + hemroi + years_from_baseline_pet                                  + (1|subjhemroi),
#                     data=tau_rois_eoad, REML=FALSE)
# mod_roi.time = lmer(value ~ -1 + hemroi + years_from_baseline_pet + hemroi:years_from_baseline_pet + (1|subjhemroi),
#                     data=tau_rois_eoad, REML=FALSE)
mod_roi      = lmer(value ~ -1 + hemroi                                                            + (1|subjhemroi),
                    data=tau_rois_eoad_all, REML=FALSE)
mod_roi_time = lmer(value ~ -1 + hemroi + years_from_baseline_pet                                  + (1|subjhemroi),
                    data=tau_rois_eoad_all, REML=FALSE)
mod_roi.time = lmer(value ~ -1 + hemroi + years_from_baseline_pet + hemroi:years_from_baseline_pet + (1|subjhemroi),
                    data=tau_rois_eoad_all, REML=FALSE)

# Test: Tau PET values change linearly with time from baseline, consistently across regions.
anova(mod_roi_time, mod_roi)
summary(mod_roi_time)[["coefficients"]]["years_from_baseline_pet",]
ci = confint(mod_roi_time, method="boot") # 95% confidence interval by parametric bootstrapping
ci

# Test: Tau PET values change linearly with time from baseline, at variable rates across regions.
anova(mod_roi.time, mod_roi_time)

# Run a separate model in each region and extract the coefficients.
coef_df = data.frame()
for (hemroi in hemrois) {
  # dat = tau_rois_eoad[(tau_rois_eoad$hemroi==hemroi),]
  dat = tau_rois_eoad_all[(tau_rois_eoad_all$hemroi==hemroi),]
  mod = lmer(value ~ 1 + years_from_baseline_pet + (1|subj), data=dat, REML=TRUE)
  ci = confint(mod, method="boot") # 95% confidence interval by parametric bootstrapping
  colnames(ci) = c("ci2.5", "ci97.5")
  coef_vals = summary(mod)[["coefficients"]]
  coef_names = rownames(coef_vals)
  for (coef_name in coef_names[1:length(coef_names)]) {
    coef_row = data.frame(hemroi = hemroi, nrow = nrow(dat), coef = coef_name, t(coef_vals[coef_name,]), t(ci[coef_name,]))
    coef_df = rbind(coef_df, coef_row) # append to the dataframe
  }
}
names(coef_df)[names(coef_df)=="Estimate"] = "est"
names(coef_df)[names(coef_df)=="Std..Error"] = "sem"
names(coef_df)[names(coef_df)=="t.value"] = "tval"
names(coef_df)[names(coef_df)=="Pr...t.."] = "pval"
coef_df = relocate(coef_df, any_of(c("ci2.5", "ci97.5")), .after = 5)

# FDR-correct P-values and determine significance.
alpha = 0.05
coef_df["pval_fdr"] = p.adjust(coef_df[,"pval"], method = "BH")
coef_df["sig_fdr"] = coef_df["pval_fdr"] < alpha

# coef_eoad_time = coef_df
# View(coef_eoad_time)
coef_eoad_time_all = coef_df
View(coef_eoad_time_all)
if (TRUE) {
  filepath = file.path(data_dir, "coef_eoad_time_all.csv")
  write.csv(coef_eoad_time_all, filepath, row.names=FALSE)
  cat("Saved", filepath)
}

# Merge the dfs.
# merged_df = merge(coef_eoad_con_bline, coef_eoad_time[(coef_eoad_time$coef=="years_from_baseline_pet"),], by="hemroi")
# merged_df = merge(distinct(tau_rois[,c("hemroi", "lobe")]), merged_df, by="hemroi")
# View(merged_df)
merged_df_all = merge(coef_eoad_con_bline, coef_eoad_time_all[(coef_eoad_time_all$coef=="years_from_baseline_pet"),], by="hemroi")
merged_df_all = merge(distinct(tau_rois[,c("hemroi", "lobe")]), merged_df_all, by="hemroi")
View(merged_df_all)
if (TRUE) {
  filepath = file.path(data_dir, "merged_df_all.csv")
  write.csv(merged_df_all, filepath, row.names=FALSE)
  cat("Saved", filepath)
}

# ***********************************************************
# -----------------------------------------------------------
# MEMs for cog/demo intearction with longitudinal tau in EOAD
# cdrsob_baseline, age_at_pet, sex_m, apoe4_alleles, fbb_centiloids_baseline
tau_rois_eoad_ = tau_rois[(tau_rois$dx=="EOAD") & (!is.na(tau_rois$cdrsob_gt4)),]
tau_rois_eoad_ = as.data.frame(unclass(tau_rois_eoad_), stringsAsFactors=TRUE)
mod_roi.time         = lmer(value ~ -1 + hemroi + years_from_baseline_pet                   + hemroi:years_from_baseline_pet                                                                                                                     + (1|subjhemroi),
                            data=tau_rois_eoad_, REML=FALSE)
mod_roi.time_cdrgt4  = lmer(value ~ -1 + hemroi + years_from_baseline_pet + cdrsob_gt4 + hemroi:years_from_baseline_pet                                                                                                                     + (1|subjhemroi),
                            data=tau_rois_eoad_, REML=FALSE)
mod_roi.time._cdrgt4 = lmer(value ~ -1 + hemroi + years_from_baseline_pet + cdrsob_gt4 + hemroi:years_from_baseline_pet + hemroi:cdrsob_gt4                                                                                            + (1|subjhemroi),
                            data=tau_rois_eoad_, REML=FALSE)
mod_roi.time_.cdrgt4 = lmer(value ~ -1 + hemroi + years_from_baseline_pet + cdrsob_gt4 + hemroi:years_from_baseline_pet                          + years_from_baseline_pet:cdrsob_gt4                                                  + (1|subjhemroi),
                            data=tau_rois_eoad_, REML=FALSE)
mod_roi.time.cdrgt4  = lmer(value ~ -1 + hemroi + years_from_baseline_pet + cdrsob_gt4 + hemroi:years_from_baseline_pet + hemroi:cdrsob_gt4 + years_from_baseline_pet:cdrsob_gt4 + hemroi:years_from_baseline_pet:cdrsob_gt4 + (1|subjhemroi),
                            data=tau_rois_eoad_, REML=FALSE)

# Test: Tau PET values covary linearly with CDR-SOB, consistently across regions.
anova(mod_roi.time_cdrgt4, mod_roi.time)
summary(mod_roi.time_cdrgt4)[["coefficients"]]

# Test: Tau PET values covary linearly with CDR-SOB, variably across regions.
anova(mod_roi.time._cdrgt4, mod_roi.time_cdrgt4)

# Test: Rate of change in tau PET covaries linearly with CDR-SOB, consistently across regions.
anova(mod_roi.time_.cdrgt4, mod_roi.time_cdrgt4)
View(summary(mod_roi.time_.cdrgt4)[["coefficients"]])

# Test: Rate of change in tau PET covaries linearly with CDR-SOB, variably across regions.
anova(mod_roi.time.cdrgt4, mod_roi.time._cdrgt4)
anova(mod_roi.time.cdrgt4, mod_roi.time_.cdrgt4)

# Run a separate model in each region and extract the coefficients.
coef_df = data.frame()
for (hemroi in hemrois) {
  dat = tau_rois_eoad_[(tau_rois_eoad_$hemroi==hemroi),]
  mod = lmer(value ~ 1 + years_from_baseline_pet * cdrsob_gt4 + (1|subj), data=dat, REML=TRUE)
  ci = confint(mod, method="boot") # 95% confidence interval by parametric bootstrapping
  colnames(ci) = c("ci2.5", "ci97.5")
  coef_vals = summary(mod)[["coefficients"]]
  coef_names = rownames(coef_vals)
  for (coef_name in coef_names[1:length(coef_names)]) {
    coef_row = data.frame(hemroi = hemroi, nrow = nrow(dat), coef = coef_name, t(coef_vals[coef_name,]), t(ci[coef_name,]))
    coef_df = rbind(coef_df, coef_row) # append to the dataframe
  }
}
names(coef_df)[names(coef_df)=="Estimate"] = "est"
names(coef_df)[names(coef_df)=="Std..Error"] = "sem"
names(coef_df)[names(coef_df)=="t.value"] = "tval"
names(coef_df)[names(coef_df)=="Pr...t.."] = "pval"
coef_df = relocate(coef_df, any_of(c("ci2.5", "ci97.5")), .after = 5)

# FDR-correct P-values and determine significance.
alpha = 0.05
coef_df["pval_fdr"] = p.adjust(coef_df[,"pval"], method = "BH")
coef_df["sig_fdr"] = coef_df["pval_fdr"] < alpha

coef_eoad_time_cdrgt4 = coef_df
View(coef_eoad_time_cdrgt4)
View(coef_eoad_time_cdrgt4[(coef_eoad_time_cdrgt4$coef=="years_from_baseline_pet:cdrsob_gt4"),])
View(coef_eoad_time_cdrgt4[(coef_eoad_time_cdrgt4$coef=="years_from_baseline_pet:cdrsob_gt4") & (coef_eoad_time_cdrgt4$pval<0.05),])
View(coef_eoad_time_cdrgt4[(coef_eoad_time_cdrgt4$coef=="cdrsob_gt4") & (coef_eoad_time_cdrgt4$pval<0.05),])
if (TRUE) {
  filepath = file.path(data_dir, "coef_eoad_time_cdrgt4.csv")
  write.csv(coef_eoad_time_cdrgt4, filepath, row.names=FALSE)
  cat("Saved", filepath)
}
# ***********************************************************

# -----------------------------------------------------------
# MEMs for cog/demo intearction with longitudinal tau in EOAD
# cdrsob_baseline, age_at_pet, sex_m, apoe4_alleles, fbb_centiloids_baseline
tau_rois_eoad_ = tau_rois[(tau_rois$dx=="EOAD") & (tau_rois$pet_visits>1) & (!is.na(tau_rois$cdrsob_baseline)),]
tau_rois_eoad_ = as.data.frame(unclass(tau_rois_eoad_), stringsAsFactors=TRUE)
mod_roi.time         = lmer(value ~ -1 + hemroi + years_from_baseline_pet                   + hemroi:years_from_baseline_pet                                                                                                                     + (1|subjhemroi),
                            data=tau_rois_eoad_, REML=FALSE)
mod_roi.time_cdrsob  = lmer(value ~ -1 + hemroi + years_from_baseline_pet + cdrsob_baseline + hemroi:years_from_baseline_pet                                                                                                                     + (1|subjhemroi),
                            data=tau_rois_eoad_, REML=FALSE)
mod_roi.time._cdrsob = lmer(value ~ -1 + hemroi + years_from_baseline_pet + cdrsob_baseline + hemroi:years_from_baseline_pet + hemroi:cdrsob_baseline                                                                                            + (1|subjhemroi),
                            data=tau_rois_eoad_, REML=FALSE)
mod_roi.time_.cdrsob = lmer(value ~ -1 + hemroi + years_from_baseline_pet + cdrsob_baseline + hemroi:years_from_baseline_pet                          + years_from_baseline_pet:cdrsob_baseline                                                  + (1|subjhemroi),
                            data=tau_rois_eoad_, REML=FALSE)
mod_roi.time.cdrsob  = lmer(value ~ -1 + hemroi + years_from_baseline_pet + cdrsob_baseline + hemroi:years_from_baseline_pet + hemroi:cdrsob_baseline + years_from_baseline_pet:cdrsob_baseline + hemroi:years_from_baseline_pet:cdrsob_baseline + (1|subjhemroi),
                            data=tau_rois_eoad_, REML=FALSE)
mod_roi.time.cdrsob2 = lmer(value ~ -1 + hemroi * years_from_baseline_pet * cdrsob_baseline + (1|subjhemroi),
                            data=tau_rois_eoad_, REML=FALSE)

# Test: Tau PET values covary linearly with CDR-SOB, consistently across regions.
anova(mod_roi.time_cdrsob, mod_roi.time)
summary(mod_roi.time_cdrsob)[["coefficients"]]["cdrsob_baseline",]

# Test: Tau PET values covary linearly with CDR-SOB, variably across regions.
anova(mod_roi.time._cdrsob, mod_roi.time_cdrsob)

# Test: Rate of change in tau PET covaries linearly with CDR-SOB, consistently across regions.
anova(mod_roi.time_.cdrsob, mod_roi.time_cdrsob)
View(summary(mod_roi.time_.cdrsob)[["coefficients"]])

# Test: Rate of change in tau PET covaries linearly with CDR-SOB, variably across regions.
anova(mod_roi.time.cdrsob, mod_roi.time._cdrsob)
anova(mod_roi.time.cdrsob, mod_roi.time_.cdrsob)

# Run a separate model in each region and extract the coefficients.
coef_df = data.frame()
for (hemroi in hemrois) {
  dat = tau_rois_eoad_[(tau_rois_eoad_$hemroi==hemroi),]
  mod = lmer(value ~ 1 + years_from_baseline_pet * cdrsob_baseline + (1|subj), data=dat, REML=TRUE)
  ci = confint(mod, method="boot") # 95% confidence interval by parametric bootstrapping
  colnames(ci) = c("ci2.5", "ci97.5")
  coef_vals = summary(mod)[["coefficients"]]
  coef_names = rownames(coef_vals)
  for (coef_name in coef_names[1:length(coef_names)]) {
    coef_row = data.frame(hemroi = hemroi, nrow = nrow(dat), coef = coef_name, t(coef_vals[coef_name,]), t(ci[coef_name,]))
    coef_df = rbind(coef_df, coef_row) # append to the dataframe
  }
}
names(coef_df)[names(coef_df)=="Estimate"] = "est"
names(coef_df)[names(coef_df)=="Std..Error"] = "sem"
names(coef_df)[names(coef_df)=="t.value"] = "tval"
names(coef_df)[names(coef_df)=="Pr...t.."] = "pval"
coef_df = relocate(coef_df, any_of(c("ci2.5", "ci97.5")), .after = 5)

# FDR-correct P-values and determine significance.
alpha = 0.05
coef_df["pval_fdr"] = p.adjust(coef_df[,"pval"], method = "BH")
coef_df["sig_fdr"] = coef_df["pval_fdr"] < alpha

coef_eoad_time_cdrsob = coef_df
View(coef_eoad_time_cdrsob)
View(coef_eoad_time_cdrsob[(coef_eoad_time_cdrsob$coef=="years_from_baseline_pet:cdrsob_baseline"),])
View(coef_eoad_time_cdrsob[(coef_eoad_time_cdrsob$coef=="years_from_baseline_pet:cdrsob_baseline") & (coef_eoad_time_cdrsob$pval<0.05),])
View(coef_eoad_time_cdrsob[(coef_eoad_time_cdrsob$coef=="cdrsob_baseline") & (coef_eoad_time_cdrsob$pval<0.05),])












# --------
# Plotting
ggplot(coef_df, aes(x="est", y="hemroi")) +
  geom_point() +
  geom_errorbarh(aes(xmin="ci2.5", xmax="ci97.5"), height=0.2)

cor.test(merged_df$est.x, merged_df$est.y)
(ggplot(merged_df, aes(est.x, est.y))
 + geom_point(size=3) + geom_smooth(method=lm)
)

ggplot(merged_df, aes(est.x, est.y)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  xlab("Baseline tau PET") +
  ylab("Change in tau PET/year") +
  scale_x_continuous(limits = c(2, 16.5), breaks = c(2, 6, 10, 14)) +
  scale_y_continuous(limits = c(-0.1, 1.5), breaks = c(0, 0.5, 1, 1.5)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

# sex, age, apoe4, cdrsob


# Post-hoc: In which regions do tau PET values increase significantly over time?
summary(mod_icpt)
summary(mod_roi)
summary(mod_roi_time)
summary(mod_roi.time)

mod_base = lmer(value ~ -1 + hemroi + hemroi:years_from_baseline_pet + (1|subjhemroi),
                data=tau_rois_eoad, REML=FALSE)
mod_base_ = lmer(value ~ -1 + hemroi + hemroi:years_from_baseline_pet + (1|subj:hemroi),
                 data=tau_rois_eoad, REML=FALSE)
print(summary(mod_base))
print(summary(mod_base_))
print(anova(mod_base, mod_base_))

print(c(length(resid(mod_base)[tau_rois_eoad$hemroi=="L_entorhinal"]^2),
        mean(resid(mod_base)[tau_rois_eoad$hemroi=="L_entorhinal"]^2)))
View(coef(mod_base)$subjhemroi[grepl("L_entorhinal", rownames(coef(mod_base)$subjhemroi)),
                               c("(Intercept)", "hemroiL_entorhinal", "hemroiL_entorhinal:years_from_baseline_pet")]
     )


mod_base = lmer(value ~ -1 + hemroi + years_from_baseline_pet +
                        hemroi:years_from_baseline_pet +
                        (1|subjhemroi),
                data=tau_rois_eoad, REML=FALSE)
mod_age  = lmer(value ~ -1 + hemroi + years_from_baseline_pet + age_at_pet +
                        hemroi:years_from_baseline_pet + age_at_pet:years_from_baseline_pet +
                        (1|subjhemroi),
                data=tau_rois_eoad, REML=FALSE)
print(anova(mod_base, mod_age))
View(summary(mod_base))
View(summary(mod_base)[["coefficients"]])
View(summary(mod_age))
View(summary(mod_age)[["coefficients"]])

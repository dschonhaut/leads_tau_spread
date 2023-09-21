library(tidyverse)
library(mgcv)
library(generalr)

# Setup

## Define parameters
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
data_file <- file.path(data_dir, "tau-eoad-re_2023-05-16.csv")
save_dir <- "/Users/dschonhaut/box/projects/leads_tau_spread/analysis/r_models/gam/bward_sel/dv_chg_suvr_re"
# today <- format(Sys.Date(), "%Y-%m-%d")
today <- "2023-06-24"
overwrite = FALSE
verbose = TRUE

## Load and format the data
check_cols <- c("subj", "roi", "suvr_bl_re", "chg_yr_re",
                "age_at_ftp_bl", "sex", "fbb_cl_bl",
                "apoe4_alleles", "cdr_sb_bl")
df <- read.csv(data_file)
df <- as.data.frame(unclass(df), stringsAsFactors = TRUE)
df <- df[complete.cases(df[,check_cols]), check_cols]
df <- data.frame(lapply(df, function(x) if (is.factor(x)) droplevels(x) else x))
df$age_at_ftp_bl_mc <- df$age_at_ftp_bl - mean(df$age_at_ftp_bl)
df$apoe4_alleles <- factor(df$apoe4_alleles, ordered = TRUE)
contrasts(df$roi) <- 'contr.sum'
contrasts(df$sex) <- 'contr.sum'

# Model fitting: backward selection

# Round 0

# Full model:
#   DV: Change in SUVR
#   IV: BL SUVR, BL Age, BL FBB-CL, BL CDR-SB, APOE4, Sex, ROI,
#       ROI x IV interactions
#   RE: Subject
mod_name <- "mod_long_full_r0_ml"
mod_long_full_r0_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r0_ml, filepath, overwrite, verbose)

# Round 1

# Subtract baseline age x ROI
mod_name <- "mod_long_full_r1_sub_age_x_roi_ml"
mod_long_full_r1_sub_age_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r1_sub_age_x_roi_ml, filepath, overwrite, verbose)

# Subtract baseline FBB-CL x ROI
mod_name <- "mod_long_full_r1_sub_fbb_x_roi_ml"
mod_long_full_r1_sub_fbb_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r1_sub_fbb_x_roi_ml, filepath, overwrite, verbose)

# Subtract baseline CDR-SB x ROI
mod_name <- "mod_long_full_r1_sub_cdr_x_roi_ml"
mod_long_full_r1_sub_cdr_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r1_sub_cdr_x_roi_ml, filepath, overwrite, verbose)

# Subtract APOE4 x ROI
mod_name <- "mod_long_full_r1_sub_apoe_x_roi_ml"
mod_long_full_r1_sub_apoe_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    # apoe4_alleles:roi +
    sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r1_sub_apoe_x_roi_ml, filepath, overwrite, verbose)

# Subtract Sex x ROI
mod_name <- "mod_long_full_r1_sub_sex_x_roi_ml"
mod_long_full_r1_sub_sex_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r1_sub_sex_x_roi_ml, filepath, overwrite, verbose)

# Compare models
mod_long_full_r0_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r0_ml", "_", today, ".rds"))
  )
mod_long_full_r1_sub_age_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r1_sub_age_x_roi_ml", "_", today, ".rds"))
  )
mod_long_full_r1_sub_fbb_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r1_sub_fbb_x_roi_ml", "_", today, ".rds"))
  )
mod_long_full_r1_sub_cdr_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r1_sub_cdr_x_roi_ml", "_", today, ".rds"))
  )
mod_long_full_r1_sub_apoe_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r1_sub_apoe_x_roi_ml", "_", today, ".rds"))
  )
mod_long_full_r1_sub_sex_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r1_sub_sex_x_roi_ml", "_", today, ".rds"))
  )

par(mfrow = c(2, 2))
gam.check(mod_long_full_r0_ml)
gam.check(mod_long_full_r1_sub_age_x_roi_ml)
gam.check(mod_long_full_r1_sub_fbb_x_roi_ml)
gam.check(mod_long_full_r1_sub_cdr_x_roi_ml)
gam.check(mod_long_full_r1_sub_apoe_x_roi_ml)
gam.check(mod_long_full_r1_sub_sex_x_roi_ml)

AIC(
  mod_long_full_r0_ml,
  mod_long_full_r1_sub_age_x_roi_ml,
  mod_long_full_r1_sub_fbb_x_roi_ml,
  mod_long_full_r1_sub_cdr_x_roi_ml,
  mod_long_full_r1_sub_apoe_x_roi_ml,
  mod_long_full_r1_sub_sex_x_roi_ml
  )

anova_test <- "F"

r1_lrt_sub_age_x_roi <- anova.gam(mod_long_full_r0_ml, mod_long_full_r1_sub_age_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r1_lrt_sub_age_x_roi", "_", today, ".rds"))
save_obj(r1_lrt_sub_age_x_roi, filepath, overwrite, verbose)

r1_lrt_sub_fbb_x_roi <- anova.gam(mod_long_full_r0_ml, mod_long_full_r1_sub_fbb_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r1_lrt_sub_fbb_x_roi", "_", today, ".rds"))
save_obj(r1_lrt_sub_fbb_x_roi, filepath, overwrite, verbose)

r1_lrt_sub_cdr_x_roi <- anova.gam(mod_long_full_r0_ml, mod_long_full_r1_sub_cdr_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r1_lrt_sub_cdr_x_roi", "_", today, ".rds"))
save_obj(r1_lrt_sub_cdr_x_roi, filepath, overwrite, verbose)

r1_lrt_sub_apoe_x_roi <- anova.gam(mod_long_full_r0_ml, mod_long_full_r1_sub_apoe_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r1_lrt_sub_apoe_x_roi", "_", today, ".rds"))
save_obj(r1_lrt_sub_apoe_x_roi, filepath, overwrite, verbose)

r1_lrt_sub_sex_x_roi <- anova.gam(mod_long_full_r0_ml, mod_long_full_r1_sub_sex_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r1_lrt_sub_sex_x_roi", "_", today, ".rds"))
save_obj(r1_lrt_sub_sex_x_roi, filepath, overwrite, verbose)

# Round 2

# Subtract sex x ROI + baseline age x ROI
mod_name <- "mod_long_full_r2_sub_age_x_roi_ml"
mod_long_full_r2_sub_age_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r2_sub_age_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + baseline FBB-CL x ROI
mod_name <- "mod_long_full_r2_sub_fbb_x_roi_ml"
mod_long_full_r2_sub_fbb_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r2_sub_fbb_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + baseline CDR-SB x ROI
mod_name <- "mod_long_full_r2_sub_cdr_x_roi_ml"
mod_long_full_r2_sub_cdr_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r2_sub_cdr_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + APOE4 x ROI
mod_name <- "mod_long_full_r2_sub_apoe_x_roi_ml"
mod_long_full_r2_sub_apoe_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    # apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r2_sub_apoe_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex
mod_name <- "mod_long_full_r2_sub_sex_ml"
mod_long_full_r2_sub_sex_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r2_sub_sex_ml, filepath, overwrite, verbose)

# Compare models
mod_long_full_r2_sub_age_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r2_sub_age_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r2_sub_fbb_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r2_sub_fbb_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r2_sub_cdr_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r2_sub_cdr_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r2_sub_apoe_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r2_sub_apoe_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r2_sub_sex_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r2_sub_sex_ml", "_", today, ".rds"))
)

par(mfrow = c(2, 2))
gam.check(mod_long_full_r2_sub_age_x_roi_ml)
gam.check(mod_long_full_r2_sub_fbb_x_roi_ml)
gam.check(mod_long_full_r2_sub_cdr_x_roi_ml)
gam.check(mod_long_full_r2_sub_apoe_x_roi_ml)
gam.check(mod_long_full_r2_sub_sex_ml)

AIC(
  mod_long_full_r1_sub_sex_x_roi_ml,
  mod_long_full_r2_sub_age_x_roi_ml,
  mod_long_full_r2_sub_fbb_x_roi_ml,
  mod_long_full_r2_sub_cdr_x_roi_ml,
  mod_long_full_r2_sub_apoe_x_roi_ml,
  mod_long_full_r2_sub_sex_ml
)

anova_test <- "F"

r2_lrt_sub_age_x_roi <- anova.gam(mod_long_full_r1_sub_sex_x_roi_ml, mod_long_full_r2_sub_age_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r2_lrt_sub_age_x_roi", "_", today, ".rds"))
save_obj(r2_lrt_sub_age_x_roi, filepath, overwrite, verbose)

r2_lrt_sub_fbb_x_roi <- anova.gam(mod_long_full_r1_sub_sex_x_roi_ml, mod_long_full_r2_sub_fbb_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r2_lrt_sub_fbb_x_roi", "_", today, ".rds"))
save_obj(r2_lrt_sub_fbb_x_roi, filepath, overwrite, verbose)

r2_lrt_sub_cdr_x_roi <- anova.gam(mod_long_full_r1_sub_sex_x_roi_ml, mod_long_full_r2_sub_cdr_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r2_lrt_sub_cdr_x_roi", "_", today, ".rds"))
save_obj(r2_lrt_sub_cdr_x_roi, filepath, overwrite, verbose)

r2_lrt_sub_apoe_x_roi <- anova.gam(mod_long_full_r1_sub_sex_x_roi_ml, mod_long_full_r2_sub_apoe_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r2_lrt_sub_apoe_x_roi", "_", today, ".rds"))
save_obj(r2_lrt_sub_apoe_x_roi, filepath, overwrite, verbose)

r2_lrt_sub_sex <- anova.gam(mod_long_full_r1_sub_sex_x_roi_ml, mod_long_full_r2_sub_sex_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r2_lrt_sub_sex", "_", today, ".rds"))
save_obj(r2_lrt_sub_sex, filepath, overwrite, verbose)

# Round 3

# Subtract sex x ROI + sex + baseline age x ROI
mod_name <- "mod_long_full_r3_sub_age_x_roi_ml"
mod_long_full_r3_sub_age_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r3_sub_age_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex + baseline FBB-CL x ROI
mod_name <- "mod_long_full_r3_sub_fbb_x_roi_ml"
mod_long_full_r3_sub_fbb_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r3_sub_fbb_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex + baseline CDR-SB x ROI
mod_name <- "mod_long_full_r3_sub_cdr_x_roi_ml"
mod_long_full_r3_sub_cdr_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r3_sub_cdr_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex + APOE4 x ROI
mod_name <- "mod_long_full_r3_sub_apoe_x_roi_ml"
mod_long_full_r3_sub_apoe_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    # apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r3_sub_apoe_x_roi_ml, filepath, overwrite, verbose)

# Compare models
mod_long_full_r3_sub_age_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r3_sub_age_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r3_sub_fbb_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r3_sub_fbb_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r3_sub_cdr_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r3_sub_cdr_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r3_sub_apoe_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r3_sub_apoe_x_roi_ml", "_", today, ".rds"))
)

par(mfrow = c(2, 2))
gam.check(mod_long_full_r3_sub_age_x_roi_ml)
gam.check(mod_long_full_r3_sub_fbb_x_roi_ml)
gam.check(mod_long_full_r3_sub_cdr_x_roi_ml)
gam.check(mod_long_full_r3_sub_apoe_x_roi_ml)

AIC(
  mod_long_full_r2_sub_sex_ml,
  mod_long_full_r3_sub_age_x_roi_ml,
  mod_long_full_r3_sub_fbb_x_roi_ml,
  mod_long_full_r3_sub_cdr_x_roi_ml,
  mod_long_full_r3_sub_apoe_x_roi_ml
)

anova_test <- "F"

r3_lrt_sub_age_x_roi <- anova.gam(mod_long_full_r2_sub_sex_ml, mod_long_full_r3_sub_age_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r3_lrt_sub_age_x_roi", "_", today, ".rds"))
save_obj(r3_lrt_sub_age_x_roi, filepath, overwrite, verbose)

r3_lrt_sub_fbb_x_roi <- anova.gam(mod_long_full_r2_sub_sex_ml, mod_long_full_r3_sub_fbb_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r3_lrt_sub_fbb_x_roi", "_", today, ".rds"))
save_obj(r3_lrt_sub_fbb_x_roi, filepath, overwrite, verbose)

r3_lrt_sub_cdr_x_roi <- anova.gam(mod_long_full_r2_sub_sex_ml, mod_long_full_r3_sub_cdr_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r3_lrt_sub_cdr_x_roi", "_", today, ".rds"))
save_obj(r3_lrt_sub_cdr_x_roi, filepath, overwrite, verbose)

r3_lrt_sub_apoe_x_roi <- anova.gam(mod_long_full_r2_sub_sex_ml, mod_long_full_r3_sub_apoe_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r3_lrt_sub_apoe_x_roi", "_", today, ".rds"))
save_obj(r3_lrt_sub_apoe_x_roi, filepath, overwrite, verbose)

# Round 4

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline age x ROI
mod_name <- "mod_long_full_r4_sub_age_x_roi_ml"
mod_long_full_r4_sub_age_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r4_sub_age_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline FBB-CL
mod_name <- "mod_long_full_r4_sub_fbb_ml"
mod_long_full_r4_sub_fbb_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    # s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r4_sub_fbb_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline CDR-SB x ROI
mod_name <- "mod_long_full_r4_sub_cdr_x_roi_ml"
mod_long_full_r4_sub_cdr_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r4_sub_cdr_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + APOE4 x ROI
mod_name <- "mod_long_full_r4_sub_apoe_x_roi_ml"
mod_long_full_r4_sub_apoe_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    # apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r4_sub_apoe_x_roi_ml, filepath, overwrite, verbose)

# Compare models
mod_long_full_r4_sub_age_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r4_sub_age_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r4_sub_fbb_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r4_sub_fbb_ml", "_", today, ".rds"))
)
mod_long_full_r4_sub_cdr_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r4_sub_cdr_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r4_sub_apoe_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r4_sub_apoe_x_roi_ml", "_", today, ".rds"))
)

par(mfrow = c(2, 2))
gam.check(mod_long_full_r4_sub_age_x_roi_ml)
gam.check(mod_long_full_r4_sub_fbb_ml)
gam.check(mod_long_full_r4_sub_cdr_x_roi_ml)
gam.check(mod_long_full_r4_sub_apoe_x_roi_ml)

AIC(
  mod_long_full_r3_sub_fbb_x_roi_ml,
  mod_long_full_r4_sub_age_x_roi_ml,
  mod_long_full_r4_sub_fbb_ml,
  mod_long_full_r4_sub_cdr_x_roi_ml,
  mod_long_full_r4_sub_apoe_x_roi_ml
)

anova_test <- "F"

r4_lrt_sub_age_x_roi <- anova.gam(mod_long_full_r3_sub_fbb_x_roi_ml, mod_long_full_r4_sub_age_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r4_lrt_sub_age_x_roi", "_", today, ".rds"))
save_obj(r4_lrt_sub_age_x_roi, filepath, overwrite, verbose)

r4_lrt_sub_fbb <- anova.gam(mod_long_full_r3_sub_fbb_x_roi_ml, mod_long_full_r4_sub_fbb_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r4_lrt_sub_fbb", "_", today, ".rds"))
save_obj(r4_lrt_sub_fbb, filepath, overwrite, verbose)

r4_lrt_sub_cdr_x_roi <- anova.gam(mod_long_full_r3_sub_fbb_x_roi_ml, mod_long_full_r4_sub_cdr_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r4_lrt_sub_cdr_x_roi", "_", today, ".rds"))
save_obj(r4_lrt_sub_cdr_x_roi, filepath, overwrite, verbose)

r4_lrt_sub_apoe_x_roi <- anova.gam(mod_long_full_r3_sub_fbb_x_roi_ml, mod_long_full_r4_sub_apoe_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r4_lrt_sub_apoe_x_roi", "_", today, ".rds"))
save_obj(r4_lrt_sub_apoe_x_roi, filepath, overwrite, verbose)

# Round 5

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline FBB-CL +
# baseline age x ROI
mod_name <- "mod_long_full_r5_sub_age_x_roi_ml"
mod_long_full_r5_sub_age_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    # s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r5_sub_age_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline FBB-CL +
# baseline CDR-SB x ROI
mod_name <- "mod_long_full_r5_sub_cdr_x_roi_ml"
mod_long_full_r5_sub_cdr_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    # s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r5_sub_cdr_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline FBB-CL +
# APOE4 x ROI
mod_name <- "mod_long_full_r5_sub_apoe_x_roi_ml"
mod_long_full_r5_sub_apoe_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    # s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    # apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r5_sub_apoe_x_roi_ml, filepath, overwrite, verbose)

# Compare models
mod_long_full_r5_sub_age_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r5_sub_age_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r5_sub_cdr_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r5_sub_cdr_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r5_sub_apoe_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r5_sub_apoe_x_roi_ml", "_", today, ".rds"))
)

par(mfrow = c(2, 2))
gam.check(mod_long_full_r5_sub_age_x_roi_ml)
gam.check(mod_long_full_r5_sub_cdr_x_roi_ml)
gam.check(mod_long_full_r5_sub_apoe_x_roi_ml)

AIC(
  mod_long_full_r4_sub_fbb_ml,
  mod_long_full_r5_sub_age_x_roi_ml,
  mod_long_full_r5_sub_cdr_x_roi_ml,
  mod_long_full_r5_sub_apoe_x_roi_ml
)

anova_test <- "F"

r5_lrt_sub_age_x_roi <- anova.gam(mod_long_full_r4_sub_fbb_ml, mod_long_full_r5_sub_age_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r5_lrt_sub_age_x_roi", "_", today, ".rds"))
save_obj(r5_lrt_sub_age_x_roi, filepath, overwrite, verbose)

r5_lrt_sub_cdr_x_roi <- anova.gam(mod_long_full_r4_sub_fbb_ml, mod_long_full_r5_sub_cdr_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r5_lrt_sub_cdr_x_roi", "_", today, ".rds"))
save_obj(r5_lrt_sub_cdr_x_roi, filepath, overwrite, verbose)

r5_lrt_sub_apoe_x_roi <- anova.gam(mod_long_full_r4_sub_fbb_ml, mod_long_full_r5_sub_apoe_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r5_lrt_sub_apoe_x_roi", "_", today, ".rds"))
save_obj(r5_lrt_sub_apoe_x_roi, filepath, overwrite, verbose)

# Round 6

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline FBB-CL +
# baseline CDR-SB x ROI + baseline age x ROI
mod_name <- "mod_long_full_r6_sub_age_x_roi_ml"
mod_long_full_r6_sub_age_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    # s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r6_sub_age_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline FBB-CL +
# baseline CDR-SB x ROI + baseline CDR-SB
mod_name <- "mod_long_full_r6_sub_cdr_ml"
mod_long_full_r6_sub_cdr_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    # s(fbb_cl_bl) +
    # s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r6_sub_cdr_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline FBB-CL +
# baseline CDR-SB x ROI + APOE4 x ROI
mod_name <- "mod_long_full_r6_sub_apoe_x_roi_ml"
mod_long_full_r6_sub_apoe_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    # s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    # apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r6_sub_apoe_x_roi_ml, filepath, overwrite, verbose)

# Compare models
mod_long_full_r6_sub_age_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r6_sub_age_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r6_sub_cdr_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r6_sub_cdr_ml", "_", today, ".rds"))
)
mod_long_full_r6_sub_apoe_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r6_sub_apoe_x_roi_ml", "_", today, ".rds"))
)

par(mfrow = c(2, 2))
gam.check(mod_long_full_r6_sub_age_x_roi_ml)
gam.check(mod_long_full_r6_sub_cdr_ml)
gam.check(mod_long_full_r6_sub_apoe_x_roi_ml)

AIC(
  mod_long_full_r5_sub_cdr_x_roi_ml,
  mod_long_full_r6_sub_age_x_roi_ml,
  mod_long_full_r6_sub_cdr_ml,
  mod_long_full_r6_sub_apoe_x_roi_ml
)

anova_test <- "F"

r6_lrt_sub_age_x_roi <- anova.gam(mod_long_full_r5_sub_cdr_x_roi_ml, mod_long_full_r6_sub_age_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r6_lrt_sub_age_x_roi", "_", today, ".rds"))
save_obj(r6_lrt_sub_age_x_roi, filepath, overwrite, verbose)

r6_lrt_sub_cdr <- anova.gam(mod_long_full_r5_sub_cdr_x_roi_ml, mod_long_full_r6_sub_cdr_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r6_lrt_sub_cdr", "_", today, ".rds"))
save_obj(r6_lrt_sub_cdr, filepath, overwrite, verbose)

r6_lrt_sub_apoe_x_roi <- anova.gam(mod_long_full_r5_sub_cdr_x_roi_ml, mod_long_full_r6_sub_apoe_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r6_lrt_sub_apoe_x_roi", "_", today, ".rds"))
save_obj(r6_lrt_sub_apoe_x_roi, filepath, overwrite, verbose)

# Round 7

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline FBB-CL +
# baseline CDR-SB x ROI + baseline CDR-SB + baseline age x ROI
mod_name <- "mod_long_full_r7_sub_age_x_roi_ml"
mod_long_full_r7_sub_age_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    # s(fbb_cl_bl) +
    # s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r7_sub_age_x_roi_ml, filepath, overwrite, verbose)

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline FBB-CL +
# baseline CDR-SB x ROI + baseline CDR-SB + APOE4 x ROI
mod_name <- "mod_long_full_r7_sub_apoe_x_roi_ml"
mod_long_full_r7_sub_apoe_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    # s(fbb_cl_bl) +
    # s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    # apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r7_sub_apoe_x_roi_ml, filepath, overwrite, verbose)

# Compare models
mod_long_full_r7_sub_age_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r7_sub_age_x_roi_ml", "_", today, ".rds"))
)
mod_long_full_r7_sub_apoe_x_roi_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r7_sub_apoe_x_roi_ml", "_", today, ".rds"))
)

par(mfrow = c(2, 2))
gam.check(mod_long_full_r7_sub_age_x_roi_ml)
gam.check(mod_long_full_r7_sub_apoe_x_roi_ml)

AIC(
  mod_long_full_r6_sub_cdr_ml,
  mod_long_full_r7_sub_age_x_roi_ml,
  mod_long_full_r7_sub_apoe_x_roi_ml
)

anova_test <- "F"

r7_lrt_sub_age_x_roi <- anova.gam(mod_long_full_r6_sub_cdr_ml, mod_long_full_r7_sub_age_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r7_lrt_sub_age_x_roi", "_", today, ".rds"))
save_obj(r7_lrt_sub_age_x_roi, filepath, overwrite, verbose)

r7_lrt_sub_apoe_x_roi <- anova.gam(mod_long_full_r6_sub_cdr_ml, mod_long_full_r7_sub_apoe_x_roi_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r7_lrt_sub_apoe_x_roi", "_", today, ".rds"))
save_obj(r7_lrt_sub_apoe_x_roi, filepath, overwrite, verbose)

# Round 8

# Test for linear associations

# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline FBB-CL +
# baseline CDR-SB x ROI + baseline CDR-SB
# Linearize age and age x ROI
mod_name <- "mod_long_full_r8_lin_age_ml"
mod_long_full_r8_lin_age_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    age_at_ftp_bl +
    # s(fbb_cl_bl) +
    # s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    age_at_ftp_bl:roi +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r8_lin_age_ml, filepath, overwrite, verbose)

# Compare models
mod_long_full_r8_lin_age_ml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r8_lin_age_ml", "_", today, ".rds"))
)

par(mfrow = c(2, 2))
gam.check(mod_long_full_r8_lin_age_ml)

AIC(
  mod_long_full_r6_sub_cdr_ml,
  mod_long_full_r8_lin_age_ml
)

anova_test <- "F"

r8_lrt_lin_age <- anova.gam(mod_long_full_r6_sub_cdr_ml, mod_long_full_r8_lin_age_ml, test = anova_test)
filepath <- file.path(save_dir, paste0("r8_lrt_lin_age", "_", today, ".rds"))
save_obj(r8_lrt_lin_age, filepath, overwrite, verbose)

# Round 9

# Refit the winning model using gam and REML
# Subtract sex x ROI + sex + baseline FBB-CL x ROI + baseline FBB-CL +
# baseline CDR-SB x ROI + baseline CDR-SB
# Linearize age and age x ROI
mod_name <- "mod_long_full_r9_lin_age_reml"
mod_long_full_r9_lin_age_reml <- gam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    age_at_ftp_bl_mc +
    # s(fbb_cl_bl) +
    # s(cdr_sb_bl) +
    apoe4_alleles +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    age_at_ftp_bl_mc:roi +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "REML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r9_lin_age_reml, filepath, overwrite, verbose)
sum_mod_long_full_r9_lin_age_reml <- summary(mod_long_full_r9_lin_age_reml)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_mod_long_full_r9_lin_age_reml, filepath, overwrite, verbose)

df$apoe4_alleles_int <- as.numeric(as.character(df$apoe4_alleles))
mod_name <- "mod_long_full_r9_lin_age_reml"
mod_long_full_r9_lin_age_reml <- gam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    age_at_ftp_bl_mc +
    # s(fbb_cl_bl) +
    # s(cdr_sb_bl) +
    apoe4_alleles_int +
    # sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    age_at_ftp_bl_mc:roi +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles_int:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "REML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r9_lin_age_reml, filepath, overwrite, verbose)
sum_mod_long_full_r9_lin_age_reml <- summary(mod_long_full_r9_lin_age_reml)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_mod_long_full_r9_lin_age_reml, filepath, overwrite, verbose)

# Check the model
mod_long_full_r9_lin_age_reml <- readRDS(
  file.path(save_dir, paste0("mod_long_full_r9_lin_age_reml", "_", today, ".rds"))
)
sum_mod_long_full_r9_lin_age_reml <- readRDS(
  file.path(save_dir, paste0("sum_mod_long_full_r9_lin_age_reml", "_", today, ".rds"))
)

par(mfrow = c(2, 2))
gam.check(mod_long_full_r9_lin_age_reml)
AIC(mod_long_full_r9_lin_age_reml)
sum_mod_long_full_r9_lin_age_reml

# ----------------
mod_name <- "mod_long_full_r9_lin_age_roifs_reml"
mod_long_full_r9_lin_age_roifs_reml <- gam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    age_at_ftp_bl_mc +
    # s(fbb_cl_bl) +
    # s(cdr_sb_bl) +
    apoe4_alleles_int +
    # sex +
    # roi +
    # ROI interactions
    s(suvr_bl_re, roi, bs = "fs") +
    age_at_ftp_bl_mc:roi +
    # s(fbb_cl_bl, by = roi, m = 1, id = 3) +
    # s(cdr_sb_bl, by = roi, m = 1, id = 4) +
    apoe4_alleles_int:roi +
    # sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "REML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r9_lin_age_roifs_reml, filepath, overwrite, verbose)
sum_mod_long_full_r9_lin_age_roifs_reml <- summary(mod_long_full_r9_lin_age_roifs_reml)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_mod_long_full_r9_lin_age_roifs_reml, filepath, overwrite, verbose)

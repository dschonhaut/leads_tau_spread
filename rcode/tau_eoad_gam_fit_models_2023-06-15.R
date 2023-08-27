library(tidyverse)
library(mgcv)
library(lme4)
library(lmerTest)
library(emmeans)
library(effects)
library(generalr)

# Setup

## Define parameters
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
data_file <- file.path(data_dir, "tau-eoad-re_2023-05-16.csv")
save_dir <- "/Users/dschonhaut/box/projects/leads_tau_spread/analysis/r_models/gam"
today <- format(Sys.Date(), "%Y-%m-%d")
overwrite = TRUE
verbose = TRUE

## Load and format the data
check_cols <- c("subj", "roi", "suvr_bl_re", "chg_yr_re",
                "age_at_ftp_bl", "sex", "fbb_cl_bl",
                "apoe4_alleles", "cdr_sb_bl")
df <- read.csv(data_file)
df <- as.data.frame(unclass(df), stringsAsFactors = TRUE)
df <- df[complete.cases(df[,check_cols]), check_cols]
df <- data.frame(lapply(df, function(x) if (is.factor(x)) droplevels(x) else x))
df$apoe4_alleles <- factor(df$apoe4_alleles, ordered = TRUE)
contrasts(df$roi) <- 'contr.sum'
contrasts(df$sex) <- 'contr.sum'

# Fit models

## DV: Baseline SUVR
### IVs: ROI
# Get the data
dat <- data.frame(df)

# Fit the models
mod_bl_roi <- lmer(suvr_bl_re ~ roi + (1|subj), data = dat, REML = FALSE)
mod_bl_icpt <- lmer(suvr_bl_re ~ 1 + (1|subj), data = dat, REML = FALSE)

# Compare models
AIC(mod_bl_roi, mod_bl_icpt)
anova(mod_bl_roi, mod_bl_icpt)

# Get partial residuals of the fixed effects
eff <- data.frame(effect("roi", mod_bl_roi))
dat <- merge(dat, eff[c("roi", "fit")], by = "roi") %>% rename(fixest = fit)
dat$resid <- resid(mod_bl_roi)
dat$partial <- dat$resid + dat$fixest
save_csv(
  dat,
  file.path(data_dir, paste0("dat_bl_roi", "_", today, ".csv")),
  overwrite = overwrite
  )

# Estimate marginal means
emm_bl_roi <- emmeans(mod_bl_roi, ~ roi, infer = TRUE, adjust = "none")

# Save the EMM dataframe
data.frame(emm_bl_roi) %>%
  arrange(emmean) %>%
  select(-SE, -df) %>%
  rename(
    suvr_bl_re = emmean,
    lower = asymp.LCL,
    upper = asymp.UCL,
    z = z.ratio,
    p_fdr = p.value
  ) %>%
  save_csv(
    file.path(data_dir, paste0("emm_bl_roi", "_", today, ".csv")),
    overwrite = overwrite
    )

data.frame(emm_bl_roi) %>%
  arrange(emmean) %>%
  select(-SE, -df) %>%
  rename(
    suvr_bl_re = emmean,
    lower = asymp.LCL,
    upper = asymp.UCL,
    z = z.ratio,
    p_fdr = p.value
    ) %>%
  plot(comparisons = TRUE)

# Calculate pairwise comparisons
emm_pairs_bl_roi <- data.frame(pairs(emm_bl_roi, adjust="BH"))
sum(emm_pairs_bl_roi["p.value"] < 0.05)/nrow(emm_pairs_bl_roi)
# pwpm(emm_bl_roi, adjust = "BH")
pwpp(emm_bl_roi, adjust = "BH")

## DV: Change in SUVR
### IVs: ROI
# Get the data
dat <- data.frame(df)

# Fit the models
mod_long_roi <- lmer(chg_yr_re ~ roi + (1|subj), data = dat, REML = FALSE)
mod_long_icpt <- lmer(chg_yr_re ~ 1 + (1|subj), data = dat, REML = FALSE)

# Compare models
AIC(mod_long_roi, mod_long_icpt)
anova(mod_long_roi, mod_long_icpt)

# Get partial residuals of the fixed effects
eff <- data.frame(effect("roi", mod_long_roi))
dat <- merge(dat, eff[c("roi", "fit")], by = "roi") %>% rename(fixest = fit)
dat$resid <- resid(mod_long_roi)
dat$partial <- dat$resid + dat$fixest
save_csv(
  dat,
  file.path(data_dir, paste0("dat_long_roi", "_", today, ".csv")),
  overwrite = overwrite
  )
# ci_long_roi <- confint(mod_long_roi, method = "boot")

# Estimate marginal means
emm_long_roi <- emmeans(mod_long_roi, ~ roi, infer = TRUE, adjust = "none")

# Save the EMM dataframe.
data.frame(emm_long_roi) %>%
  arrange(emmean) %>%
  select(-SE, -df) %>%
  rename(
    chg_yr_re = emmean,
    lower = asymp.LCL,
    upper = asymp.UCL,
    z = z.ratio,
    p_fdr = p.value
    ) %>%
  save_csv(
    file.path(data_dir, paste0("emm_long_roi", "_", today, ".csv")),
    overwrite = TRUE
    )

# Calculate pairwise comparisons
emm_pairs_long_roi <- data.frame(pairs(emm_long_roi, adjust="BH"))
sum(emm_pairs_long_roi["p.value"] < 0.05)/nrow(emm_pairs_long_roi)
# pwpm(emm_long_roi, adjust = "BH")
pwpp(emm_long_roi, adjust = "BH")

## DV: Change in SUVR
### IVs: Baseline SUVR, ROI
mod_name <- "mod_long_lbl_ml"
mod_long_lbl_ml <- bam(
  chg_yr_re ~
    # Main effects
    suvr_bl_re +
    roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_lbl_ml, filepath, overwrite, verbose)
sum_mod_long_lbl_ml <- summary(mod_long_lbl_ml)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_mod_long_lbl_ml, filepath, overwrite, verbose)

mod_name <- "mod_long_sbl_ml"
mod_long_sbl_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_sbl_ml, filepath, overwrite, verbose)
sum_mod_long_sbl_ml <- summary(mod_long_sbl_ml)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_mod_long_sbl_ml, filepath, overwrite, verbose)

mod_name <- "mod_long_lbl_x_roi_ml"
mod_long_lbl_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    suvr_bl_re +
    roi +
    # ROI interactions
    suvr_bl_re:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_lbl_x_roi_ml, filepath, overwrite, verbose)
sum_mod_long_lbl_x_roi_ml <- summary(mod_long_lbl_x_roi_ml)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_mod_long_lbl_x_roi_ml, filepath, overwrite, verbose)

mod_name <- "mod_long_sbl_x_roi_ml"
mod_long_sbl_x_roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_sbl_x_roi_ml, filepath, overwrite, verbose)
sum_mod_long_sbl_x_roi_ml <- summary(mod_long_sbl_x_roi_ml)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_mod_long_sbl_x_roi_ml, filepath, overwrite, verbose)

### Tests
par(mfrow = c(2, 2))
gam.check(mod_long_sbl_ml)
gam.check(mod_long_sbl_x_roi_ml)
gam.check(mod_long_lbl_x_roi_ml)

#### Is there a baseline SUVR x ROI interaction?
anova.gam(mod_long_sbl_ml, mod_long_sbl_x_roi_ml, test = "LRT")

#### Is change in SUVR ~ baseline SUVR a linear association?
anova.gam(mod_long_lbl_x_roi_ml, mod_long_sbl_x_roi_ml, test = "LRT")

### IVs: Baseline age, Baseline SUVR, ROI
mod_name <- ""
mod_long_sage_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

mod_name <- ""
mod_long_sage.x.roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

mod_name <- ""
mod_long_age_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    age_at_ftp_bl +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

mod_name <- ""
mod_long_age.x.roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    age_at_ftp_bl +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    age_at_ftp_bl:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

### IVs: Baseline FBB-CL, Baseline SUVR, ROI
mod_long_sfbb_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(fbb_cl_bl) +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

mod_long_sfbb.x.roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(fbb_cl_bl) +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(fbb_cl_bl, by = roi, m = 1, id = 2) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

mod_long_fbb_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    fbb_cl_bl +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

mod_long_fbb.x.roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    fbb_cl_bl +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    fbb_cl_bl:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

### IVs: Baseline CDR-SB, Baseline SUVR, ROI
mod_long_scdr_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

mod_long_scdr.x.roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    s(age_at_ftp_bl) +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    s(age_at_ftp_bl, by = roi, m = 1, id = 2) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

mod_long_age_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    age_at_ftp_bl +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

mod_long_age.x.roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    age_at_ftp_bl +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    age_at_ftp_bl:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

### IVs: APOE4 allele count, Baseline SUVR, ROI
mod_long_apoe_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    apoe4_alleles +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

mod_long_apoe.x.roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    apoe4_alleles +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    apoe4_alleles:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

### IVs: Sex, Baseline SUVR, ROI
mod_long_sex_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

mod_long_sex.x.roi_ml <- bam(
  chg_yr_re ~
    # Main effects
    s(suvr_bl_re) +
    sex +
    roi +
    # ROI interactions
    s(suvr_bl_re, by = roi, m = 1, id = 1) +
    sex:roi +
    # Random effect covariates
    s(subj, bs = "re"),
  data = df,
  method = "ML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(xx, filepath, overwrite, verbose)
sum_xx <- summary(xx)
filepath <- file.path(save_dir, paste0("sum_", mod_name, "_", today, ".rds"))
save_obj(sum_xx, filepath, overwrite, verbose)

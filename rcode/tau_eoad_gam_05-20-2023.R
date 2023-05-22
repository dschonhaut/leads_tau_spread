library(tidyverse)
library(lme4)
library(lmerTest)
library(mgcv)
library(gratia)
library(viridis)
library(broom)

# ------------------------------------------------------------------------------
# Import and format the data
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
tau <- read.csv(file.path(data_dir, "tau-eoad-re_2023-05-16.csv"))
tau <- as.data.frame(unclass(tau), stringsAsFactors = TRUE)
tau$apoe4_alleles <- factor(tau$apoe4_alleles, ordered = TRUE)
subjs <- levels(as.factor(tau$subj))
rois <- levels(as.factor(tau$roi))
blue2 <- "#2E45B8"
red2 <- "#E3120B"

# Model the data.
check_cols <- c("subj", "roi", "suvr_bl_re", "chg_yr_re",
                "age_at_ftp_bl", "sex", "apoe4_alleles",
                "fbb_cl_bl", "cdr_sb_bl")
df <- tau
df <- df[complete.cases(df[,check_cols]),]

# -----------------------------------------------
mod_bltau_roi <- gam(
  chg_yr_re ~
    s(suvr_bl_re, bs="tp") +
    s(roi, bs="re") +
    s(suvr_bl_re, by=roi, bs="tp", m=1) +
    s(subj, bs="re"),
  data=df,
  method="ML",
  family=gaussian,
  select=TRUE
)
layout(matrix(1:4,ncol=2,byrow=TRUE)); gam.check(mod_bltau_roi) + abline(0, 1, col=red2); layout(1)
summary(mod_bltau_roi)
plot.gam(
  mod_bltau_roi,
  unconditional=T,
  scale=0,
  residuals=T,
  seWithMean=T,
  # shift=fixef(mod1)["(Intercept)"],
  shade=T,
  scheme=2,
  pages=4
)
plot.gam(
  mod_bltau_roi,
  unconditional=T,
  scale=0,
  residuals=T,
  seWithMean=T,
  shift=fixef(mod_bltau_roi)["(Intercept)"],
  shade=T,
  scheme=2,
  select=1
)

# -----------------------------------------------
mod_bltau_roi_age <- gam(
  chg_yr_re ~
    s(suvr_bl_re, bs="tp") +
    s(age_at_ftp_bl, bs="tp") +
    s(roi, bs="re") +
    s(suvr_bl_re, by=roi, bs="tp", m=1) +
    s(age_at_ftp_bl, by=roi, bs="tp", m=1) +
    s(subj, bs="re"),
  data=df,
  method="ML",
  family=gaussian,
  select=TRUE
)
layout(matrix(1:4,ncol=2,byrow=T)); gam.check(mod_bltau_roi_age) + abline(0, 1, col=red2); layout(1)
View(concurvity(mod_bltau_roi_age, full=F)[["estimate"]])
View(concurvity(mod_bltau_roi_age, full=T))
summary(mod_bltau_roi_age)
appraise(mod_bltau_roi_age, method="simulate")
View(smooth_estimates(mod_bltau_roi_age, "s(roi)") %>% add_confint() %>% arrange(-est))
View(smooth_estimates(mod_bltau_roi_age, "s(suvr_bl_re)") %>% add_confint() %>% arrange(-est))
View(smooth_estimates(mod_bltau_roi_age, "s(age_at_ftp_bl)") %>% add_confint() %>% arrange(-est))
plot.gam(
  mod_bltau_roi_age,
  unconditional=T,
  scale=0,
  residuals=T,
  seWithMean=T,
  # shift=fixef(mod1)["(Intercept)"],
  shade=T,
  scheme=2,
  pages=6
)
layout(1)
plot.gam(
  mod_bltau_roi_age,
  unconditional=T,
  scale=0,
  residuals=T,
  seWithMean=T,
  shift=fixef(mod_bltau_roi_age)["(Intercept)"],
  shade=T,
  scheme=2,
  select=2,
  # xlab="Baseline FTP SUVR",
  ylab="FTP SUVR change/year"
) + abline(h=0, col=blue2)
draw(mod_bltau_roi_age,
     scales="fixed",
     residuals=T,
     # shift=fixef(mod_bltau_roi_age)["(Intercept)"],
     # size=1,
     # coverage=1,
     overall_uncertainty=T,
     select=3) #+ ylim(-0.6, 0.6)

# -----------------------------------------------
mod_full_lin <- gam(
  chg_yr_re ~
    # Fixed effects
    roi +
    apoe4_alleles +
    sex +
    suvr_bl_re +
    age_at_ftp_bl +
    fbb_cl_bl +
    cdr_sb_bl +
    # Interactions
    apoe4_alleles:roi +
    sex:roi +
    suvr_bl_re:roi +
    age_at_ftp_bl:roi +
    fbb_cl_bl:roi +
    cdr_sb_bl:roi +
    # Random effects
    s(subj, bs="re"),
  data=df,
  method="ML",
  family=gaussian,
  select=TRUE
)
mod_full_gam <- gam(
  chg_yr_re ~
    # Fixed effects: discrete factors
    roi +
    apoe4_alleles +
    sex +
    # Fixed effects: smooth
    s(suvr_bl_re, bs="tp") +
    s(age_at_ftp_bl, bs="tp") +
    s(fbb_cl_bl, bs="tp") +
    s(cdr_sb_bl, bs="tp") +
    # Factor/smooth interactions
    s(apoe4_alleles, by=roi, m=1, bs="tp") +
    s(sex, by=roi, m=1, bs="tp") +
    s(suvr_bl_re, by=roi, m=1, bs="tp") +
    s(age_at_ftp_bl, by=roi, m=1, bs="tp") +
    s(fbb_cl_bl, by=roi, m=1, bs="tp") +
    s(cdr_sb_bl, by=roi, m=1, bs="tp") +
    # Random effects
    s(subj, bs="re"),
  data=df,
  method="ML",
  family=gaussian,
  select=TRUE
)
layout(matrix(1:4,ncol=2,byrow=TRUE)); gam.check(mod_full_gam) + abline(0, 1, col="red"); layout(1)
summary(mod_full_gam)
summary(mod_full_lin)
anova(mod_full_gam, mod_full_lin)
plot.gam(
  mod_bltau_roi_age,
  residuals=F,
  seWithMean=T,
  # shift=fixef(mod1)["(Intercept)"],
  shade=T,
  scheme=2,
  pages=4
)

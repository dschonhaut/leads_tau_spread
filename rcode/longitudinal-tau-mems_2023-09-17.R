library(tidyverse)
library(gridExtra)
library(patchwork)
library(broom)
library(viridis)
library(performance)
library(sjPlot)
library(sjstats)
library(gtsummary)
library(emmeans)
library(dlookr)
library(boot)
library(lme4)
library(lmerTest)
library(mgcv)
library(gratia)

# ------------------------------------------------------------------------------
# Import and format the data
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
save_dir <- "/Users/dschonhaut/box/projects/leads_tau_spread/analysis/r_models/gam/idv_rois"
today <- format(Sys.Date(), "%Y-%m-%d")
overwrite = FALSE
verbose = TRUE
tau_eoad_long <- read_csv(file.path(data_dir, "tau-rois-agg_eoad-long_formatted_2023-10-22.csv"))
tau_eoad_long <- tau_eoad_long %>% filter(parc %in% c("metarois", "fsroi_bilat"))
# mean_age = (tau_eoad_long %>% distinct(subj, .keep_all=T) %>% summarize(age_at_ftp_bl_mc=mean(age_at_ftp_bl)))[[1]]
# tau_eoad_long$age_at_ftp_bl_mc <- tau_eoad_long$age_at_ftp_bl - mean_age

# tau_eoad_long <- tau_eoad_long[(tau_eoad_long["parc"] == "roi24"), ]

# Convert strings to factors
tau_eoad_long <- as.data.frame(unclass(tau_eoad_long), stringsAsFactors = TRUE)
tau_eoad_long$apoe4_alleles <- factor(tau_eoad_long$apoe4_alleles, ordered = TRUE)
contrasts(tau_eoad_long$sex) <- 'contr.sum'

eoad_subjs <- levels(as.factor(tau_eoad_long$subj))
rois <- levels(as.factor(tau_eoad_long$roi))

# Bootstrap first attempt
# ------------------------------------------------------------------------------


# Fit mixed effects models
# ------------------------------------------------------------------------------
# Store the input dataframe as df.
df <- tibble(tau_eoad_long)


# check_cols <- c("subj", "roi", "suvr",
#                 "age_at_ftp_bl", "age_at_ftp_bl_mc", "sex", "fbb_cl_bl",
#                 "apoe4_alleles", "cdr_sb_bl", "mmse_bl", "suvr_bl", "ftp_yrs_from_bl")
# df <- df[complete.cases(df[,check_cols]), ]
# df <- df %>%
#   filter(
#     dx == "EOAD",
#     ftp_visits >= 2
#   )

# Initialize an empty dataframe to store the results..
mm_fits <- data.frame(
  subj = character(),
  roi = character(),
  fixed_icpt = numeric(),
  fixed_time = numeric(),
  rdm_icpt = numeric(),
  rdm_time = numeric(),
  stringsAsFactors = FALSE
)

# For each roi in df$roi, fit a linear mixed-effects model like this:
# suvr ~ ftp_yrs_from_bl + (1|subj) + (0 + ftp_yrs_from_bl|subj)
# then use the fitted model to append a row to the mm_fits dataframe
# for each subject and roi.
# rois <- c("temporal", "frontal")
for (roi in rois) {
  # Select data for the current ROI
  df_roi <- df[df$roi == roi, ]

  # Fit the linear mixed-effects model
  # mm <- lmer(
  #   suvr ~ ftp_yrs_from_bl + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  #   data = df_roi
  # )
  mm <- lmer(
    suvr_roi_mc ~ ftp_yrs_from_bl + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
    data = df_roi
  )

  # Extract fixed and random effects
  mm_fixefs <- fixef(mm)
  mm_ranefs <- ranef(mm)$subj

  # Create a dataframe for the current roi's results
  mm_fits_roi <- data.frame(
    subj = rownames(mm_ranefs),
    roi = roi,
    fixed_icpt = mm_fixefs[["(Intercept)"]],
    fixed_time = mm_fixefs[["ftp_yrs_from_bl"]],
    rdm_icpt = mm_ranefs[, "(Intercept)"],
    rdm_time = mm_ranefs[, "ftp_yrs_from_bl"]
  )

  # Append the results to the 'mm_fits' dataframe
  mm_fits <- bind_rows(mm_fits, mm_fits_roi)
}

mm_fits <- mm_fits %>%
  mutate(
    icpt = fixed_icpt + rdm_icpt,
    time = fixed_time + rdm_time
  )
mm_fits <- mm_fits %>%
  inner_join(
    df %>%
      distinct(subj, roi, .keep_all=TRUE),
    by = c("subj", "roi")
    ) %>%
  select(
    subj, parc, roi, dx, ftp_visits,
    age_at_ftp_bl, age_at_ftp_bl_mc, sex,
    apoe4_alleles, apoe4_alleles_mc,
    fbb_cl_bl, fbb_cl_bl_mc,
    mmse_bl, mmse_bl_mc,
    cdr_sb_bl, cdr_sb_bl_mc,
    fixed_icpt, fixed_time,
    rdm_icpt, rdm_time,
    icpt, time
  )

# Preview summary stats
mm_fits %>%
  group_by(roi) %>%
  summarize(
    n_subj = n_distinct(subj),
    min_icpt = min(icpt),
    max_icpt = max(icpt),
    mean_icpt = mean(icpt),
    sd_icpt = sd(icpt),
    min_time = min(time),
    max_time = max(time),
    mean_time = mean(time),
    sd_time = sd(time)
  ) %>%
  arrange(-mean_time) %>%
  view()

# Save the dataframe with subject- and ROI-specific SUVR intercepts and slopes
mm_outputf <- file.path(data_dir, "tau-rois-mms_all-eoad-gt1ftp_2023-10-19.csv")
write_csv(
  mm_fits,
  mm_outputf,
)

# mean suvr chg/yr from v1 to v2
mean((df_roi[df_roi$visit==2,"suvr"] - df_roi[df_roi$visit==2,"suvr_bl"]) / df_roi[df_roi$visit==2,"ftp_yrs_from_bl"])

# Single additional predictors
roi <- "caudalmiddlefrontal"

df <- data.frame(tau_eoad_long[complete.cases(tau_eoad_long[,c("suvr")]) & (tau_eoad_long$roi==roi), ])
mm <- lmer(
  suvr ~ ftp_yrs_from_bl + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  data = df
)
summary(mm)

df <- data.frame(tau_eoad_long[complete.cases(tau_eoad_long[,c("age_at_ftp_bl_mc")]) & (tau_eoad_long$roi==roi), ])
mm_age <- lmer(
  suvr ~ ftp_yrs_from_bl * age_at_ftp_bl_mc + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  data = df
)
summary(mm_age)

df <- data.frame(tau_eoad_long[complete.cases(tau_eoad_long[,c("sex")]) & (tau_eoad_long$roi==roi), ])
mm_sex <- lmer(
  suvr ~ ftp_yrs_from_bl * sex + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  data = df
)
summary(mm_sex)

df <- data.frame(tau_eoad_long[complete.cases(tau_eoad_long[,c("apoe4_alleles_mc")]) & (tau_eoad_long$roi==roi), ])
mm_apoe <- lmer(
  suvr ~ ftp_yrs_from_bl * apoe4_alleles_mc + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  data = df
)
summary(mm_apoe)

df <- data.frame(tau_eoad_long[complete.cases(tau_eoad_long[,c("fbb_cl_bl_mc")]) & (tau_eoad_long$roi==roi), ])
mm_cl <- lmer(
  suvr ~ ftp_yrs_from_bl * fbb_cl_bl_mc + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  data = df
)
summary(mm_cl)

df <- data.frame(tau_eoad_long[complete.cases(tau_eoad_long[,c("mmse_bl_mc")]) & (tau_eoad_long$roi==roi), ])
mm_mmse <- lmer(
  suvr ~ ftp_yrs_from_bl * mmse_bl_mc + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  data = df
)
summary(mm_mmse)

df <- data.frame(tau_eoad_long[complete.cases(tau_eoad_long[,c("cdr_sb_bl_mc")]) & (tau_eoad_long$roi==roi), ])
mm_cdr <- lmer(
  suvr ~ ftp_yrs_from_bl * cdr_sb_bl_mc + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  data = df
)
summary(mm_cdr)

df <- data.frame(tau_eoad_long[complete.cases(tau_eoad_long[,c("age_at_ftp_bl_mc", "sex", "apoe4_alleles_mc", "fbb_cl_bl_mc", "mmse_bl_mc")]) & (tau_eoad_long$roi==roi), ])
mm_all <- lmer(
  suvr ~
    ftp_yrs_from_bl +
    age_at_ftp_bl_mc +
    sex +
    apoe4_alleles_mc +
    fbb_cl_bl_mc +
    mmse_bl_mc +
    ftp_yrs_from_bl:age_at_ftp_bl_mc +
    ftp_yrs_from_bl:sex +
    ftp_yrs_from_bl:apoe4_alleles_mc +
    ftp_yrs_from_bl:fbb_cl_bl_mc +
    ftp_yrs_from_bl:mmse_bl_mc +
    (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  data = df
)
summary(mm_all)

df <- data.frame(tau_eoad_long[complete.cases(tau_eoad_long[,c("age_at_ftp_bl_mc", "sex", "apoe4_alleles_mc", "fbb_cl_bl_mc", "mmse_bl_mc")]) & (tau_eoad_long$roi==roi), ])
mm_all2 <- lmer(
  suvr ~
    ftp_yrs_from_bl +
    age_at_ftp_bl_mc +
    sex +
    apoe4_alleles_mc +
    fbb_cl_bl_mc +
    mmse_bl_mc +
    ftp_yrs_from_bl:age_at_ftp_bl_mc +
    ftp_yrs_from_bl:sex +
    ftp_yrs_from_bl:apoe4_alleles_mc +
    ftp_yrs_from_bl:fbb_cl_bl_mc +
    ftp_yrs_from_bl:mmse_bl_mc +
    age_at_ftp_bl_mc:sex +
    age_at_ftp_bl_mc:apoe4_alleles_mc +
    age_at_ftp_bl_mc:fbb_cl_bl_mc +
    age_at_ftp_bl_mc:mmse_bl_mc +
    sex:apoe4_alleles_mc +
    sex:fbb_cl_bl_mc +
    sex:mmse_bl_mc +
    apoe4_alleles_mc:fbb_cl_bl_mc +
    apoe4_alleles_mc:mmse_bl_mc +
    fbb_cl_bl_mc:mmse_bl_mc +
    (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  data = df
)
summary(mm_all2)

df <- data.frame(tau_eoad_long[complete.cases(tau_eoad_long[,c("age_at_ftp_bl_mc", "sex", "apoe4_alleles_mc", "fbb_cl_bl_mc", "mmse_bl_mc")]) & (tau_eoad_long$roi==roi), ])
mm_all3 <- lmer(
  suvr ~
    ftp_yrs_from_bl +
    age_at_ftp_bl_mc +
    sex +
    apoe4_alleles_mc +
    fbb_cl_bl_mc +
    mmse_bl_mc +
    ftp_yrs_from_bl:age_at_ftp_bl_mc +
    ftp_yrs_from_bl:sex +
    ftp_yrs_from_bl:apoe4_alleles_mc +
    ftp_yrs_from_bl:fbb_cl_bl_mc +
    ftp_yrs_from_bl:mmse_bl_mc +
    age_at_ftp_bl_mc:sex +
    age_at_ftp_bl_mc:apoe4_alleles_mc +
    age_at_ftp_bl_mc:fbb_cl_bl_mc +
    age_at_ftp_bl_mc:mmse_bl_mc +
    sex:apoe4_alleles_mc +
    sex:fbb_cl_bl_mc +
    sex:mmse_bl_mc +
    apoe4_alleles_mc:fbb_cl_bl_mc +
    apoe4_alleles_mc:mmse_bl_mc +
    fbb_cl_bl_mc:mmse_bl_mc +
    ftp_yrs_from_bl:age_at_ftp_bl_mc:sex +
    ftp_yrs_from_bl:age_at_ftp_bl_mc:apoe4_alleles_mc +
    ftp_yrs_from_bl:age_at_ftp_bl_mc:fbb_cl_bl_mc +
    ftp_yrs_from_bl:age_at_ftp_bl_mc:mmse_bl_mc +
    ftp_yrs_from_bl:sex:apoe4_alleles_mc +
    ftp_yrs_from_bl:sex:fbb_cl_bl_mc +
    ftp_yrs_from_bl:sex:mmse_bl_mc +
    ftp_yrs_from_bl:apoe4_alleles_mc:fbb_cl_bl_mc +
    ftp_yrs_from_bl:apoe4_alleles_mc:mmse_bl_mc +
    ftp_yrs_from_bl:fbb_cl_bl_mc:mmse_bl_mc +
    (1 | subj) + (0 + ftp_yrs_from_bl | subj),
  data = df
)
summary(mm_all3)

tab_model(mm, mm_all, mm_all2, mm_all3, mm_age, mm_sex, mm_apoe, mm_cl, mm_mmse, mm_cdr,
          title=roi, digits=4, CSS = css_theme("cells"))

check_model(mm)

# ------------------------------------------------------------------------------
# GAM models
# Full model:
#   DV: Change in SUVR
#   IV: BL SUVR, BL Age, BL FBB-CL, BL CDR-SB, APOE4, Sex, ROI,
#       ROI x IV interactions
#   RE: Subject
## Load and format the data
mm_outputf <- file.path(data_dir, "tau-roi24-mms_all-eoad-gt1ftp_2023-10-18.csv")
check_cols <- c("age_at_ftp_bl", "sex", "fbb_cl_bl", "apoe4_alleles", "cdr_sb_bl", "mmse_bl")
df <- read.csv(mm_outputf)
df <- as.data.frame(unclass(df), stringsAsFactors = TRUE)
df <- df[complete.cases(df[,check_cols]), ]
df <- data.frame(lapply(df, function(x) if (is.factor(x)) droplevels(x) else x))
df$apoe4_alleles <- factor(df$apoe4_alleles, ordered = TRUE)
# contrasts(df$roi) <- 'contr.sum'
contrasts(df$sex) <- 'contr.sum'

roi <- "temporal"
df_roi <- as.data.frame(df[df$roi==roi,])
mod_name <- paste("mod_long_full_reml", roi, sep="_")
mod_long_full_temporal <- gam(
  time ~
    # Main effects
    s(icpt) +
    s(age_at_ftp_bl_mc) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex,
  data=df_roi,
  select=FALSE,
  method="REML"
)
mod_long_subicpt_temporal <- gam(
  time ~
    # Main effects
    # s(icpt) +
    s(age_at_ftp_bl_mc) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex,
  data=df_roi,
  select=FALSE,
  method="REML"
)
plot.gam(mod_long_full_temporal, residuals=T, shift=0.0836445, seWithMean=T, shade=T)

roi <- "frontal"
df_roi <- as.data.frame(df[df$roi==roi,])
mod_name <- paste("mod_long_full_reml", roi, sep="_")
mod_long_full_frontal <- gam(
  time ~
    # Main effects
    s(icpt) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex,
  data=df_roi,
  select=FALSE,
  method="REML"
)
mod_long_subicpt_frontal <- gam(
  time ~
    # Main effects
    # s(icpt) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex,
  data=df_roi,
  select=FALSE,
  method="REML"
)

mod_long_lin_reml <- bam(
  time ~
    # Main effects
    s(icpt) +
    age_at_ftp_bl +
    fbb_cl_bl +
    cdr_sb_bl +
    apoe4_alleles +
    sex,
  data = df,
  method = "fREML"
)
filepath <- file.path(save_dir, paste0(mod_name, "_", today, ".rds"))
save_obj(mod_long_full_r0_ml, filepath, overwrite, verbose)

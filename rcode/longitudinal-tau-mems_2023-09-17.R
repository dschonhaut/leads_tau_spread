library(tidyverse)
library(lme4)
library(lmerTest)
library(viridis)
library(broom)
library(gtsummary)
library(emmeans)
library(mgcv)
library(gratia)

# ------------------------------------------------------------------------------
# Import and format the data
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
save_dir <- "/Users/dschonhaut/box/projects/leads_tau_spread/analysis/r_models/gam/idv_rois"
today <- format(Sys.Date(), "%Y-%m-%d")
overwrite = FALSE
verbose = TRUE
tau_eoad_long <- read.csv(file.path(data_dir, "tau-rois-agg_eoad-long_2023-09-17.csv"))
tau_eoad_long <- tau_eoad_long[(tau_eoad_long["parc"] == "roi24"), ]

# Convert strings to factors
tau_eoad_long <- as.data.frame(unclass(tau_eoad_long), stringsAsFactors = TRUE)
tau_eoad_long$apoe4_alleles <- factor(tau_eoad_long$apoe4_alleles, ordered = TRUE)

eoad_subjs <- levels(as.factor(tau_eoad_long$subj))
rois <- levels(as.factor(tau_eoad_long$roi))

# Fit mixed effects models
# ------------------------------------------------------------------------------
# Store the input dataframe as df.
df <- data.frame(tau_eoad_long)
df <- df %>%
  filter(
    dx == "EOAD",
    ftp_visits > 1
  )

# Initialize an empty dataframe to store the results..
results <- data.frame(
  subj = character(),
  roi = character(),
  fixed_icpt = numeric(),
  fixed_slope_time = numeric(),
  rndm_icpt = numeric(),
  rndm_slope_time = numeric(),
  stringsAsFactors = FALSE
)

# For each roi in df$roi, fit a linear mixed-effects model like this:
# suvr ~ ftp_yrs_from_bl + (1|subj) + (0 + ftp_yrs_from_bl|subj)
# then use the fitted model to append a row to the results dataframe
# for each subject and roi.
# Loop through each roi and fit mixed effects models
for (roi in rois) {
  # Select data for the current ROI
  df_roi <- df[df$roi == roi, ]

  # Fit the linear mixed-effects model
  mem <- lmer(
    suvr ~ ftp_yrs_from_bl + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
    data = df_roi
  )

  # Extract fixed and random effects
  mem_fefs <- fixef(mem)
  mem_refs <- ranef(mem)$subj

  # Create a dataframe for the current roi's results
  results_roi <- data.frame(
    subj = rownames(mem_refs),
    roi = roi,
    fixed_icpt = mem_fefs["(Intercept)"],
    fixed_slope_time = mem_fefs["ftp_yrs_from_bl"],
    rndm_icpt = mem_refs[, "(Intercept)"],
    rndm_slope_time = mem_refs[, "ftp_yrs_from_bl"]
  )

  # Append the results to the 'results' dataframe
  results <- bind_rows(results, results_roi)
}

results <- results %>%
  mutate(
    icpt = fixed_icpt + rndm_icpt,
    slope = fixed_slope_time + rndm_slope_time
  )
results <- results %>%
  inner_join(
    df %>%
      group_by(subj, roi) %>%
      summarize(
        parc = "roi24",
        dx = unique(dx),
        age_at_ftp_bl = unique(age_at_ftp_bl),
        sex = unique(sex),
        ftp_visits = unique(ftp_visits),
        ftp_yrs_from_bl = max(ftp_yrs_from_bl),
        apoe4_alleles = unique(apoe4_alleles),
        fbb_cl_bl = unique(fbb_cl_bl),
        mmse_bl = unique(mmse_bl),
        cdr_sb_bl = unique(cdr_sb_bl),
      ),
    by = c("subj", "roi")
    ) %>%
  select(
    subj, roi, parc, dx, age_at_ftp_bl, sex, ftp_visits, ftp_yrs_from_bl,
    apoe4_alleles, fbb_cl_bl, mmse_bl, cdr_sb_bl, icpt, slope, fixed_icpt,
    fixed_slope_time, rndm_icpt, rndm_slope_time
  )

# Save the dataframe with subject- and ROI-specific SUVR intercepts and slopes
write.csv(
  results,
  file.path(data_dir, "tau-roi24-mems_all-eoad-gt1ftp_2023-09-18.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# GAM models
# Full model:
#   DV: Change in SUVR
#   IV: BL SUVR, BL Age, BL FBB-CL, BL CDR-SB, APOE4, Sex, ROI,
#       ROI x IV interactions
#   RE: Subject
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

roi <- "prcu"
df <- as.data.frame(results[results$roi==roi,])
mod_name <- paste("mod_long_full_reml", roi, sep="_")
mod_long_full_reml <- bam(
  slope ~
    # Main effects
    s(icpt) +
    s(age_at_ftp_bl) +
    s(fbb_cl_bl) +
    s(cdr_sb_bl) +
    apoe4_alleles +
    sex,
  data = df,
  select=TRUE,
  method = "fREML"
)
mod_long_lin_reml <- bam(
  slope ~
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

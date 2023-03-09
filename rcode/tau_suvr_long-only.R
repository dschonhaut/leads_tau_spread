# `radian` in Terminal to get R shell
library(lme4)
library(lmerTest)
library(ggplot2)
library(tidyverse)

# --------------------------
# Import and format the data
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
filepath <- file.path(data_dir, "tau_roi_suvrs_eoad_long.csv")
tau_suvrs <- read.csv(filepath)
# Replace missing values with NA
cols_to_replace <- c("minority", "apoe_genotype", "apoe4_alleles", "cdrsb_bline", "mmse_baseline")
tau_suvrs[cols_to_replace] <- lapply(tau_suvrs[cols_to_replace], function(x) {
  x[x == ""] <- NA
  x
})

# Convert strings to factors
tau_suvrs <- as.data.frame(unclass(tau_suvrs), stringsAsFactors = TRUE)

eoad_subjs <- levels(as.factor(tau_suvrs$subj))
hemrois <- levels(as.factor(tau_suvrs$hemroi))
subjhemrois <- levels(as.factor(tau_suvrs$subjhemroi))

# --------------------------------------
# MEMs for longitudinal analysis in EOAD
hemroi <- "L_inferiorparietal"
hemroi <- "L_frontalpole"
tau_slopes <- tau_suvrs[, -which(names(tau_suvrs) == "delta_suvr")]
ranef_df <- data.frame()
coef_df <- data.frame()
for (hemroi in hemrois) {
  # Get the data for this hemroi.
  dat <- tau_suvrs[(tau_suvrs$hemroi == hemroi), ]

  # Fit the model.
  mod <- lmer(delta_suvr ~ 1 + (1 | subj), data = dat, REML = TRUE)

  # Append the fixed effects.
  ci <- confint(mod, method = "boot") # 95% confidence interval by parametric bootstrapping
  colnames(ci) <- c("ci2.5", "ci97.5")
  coef_vals <- summary(mod)[["coefficients"]]
  coef_names <- rownames(coef_vals)
  for (coef_name in coef_names[1:length(coef_names)]) {
    coef_row <- data.frame(hemroi = hemroi, nobs = nrow(dat), coef = coef_name, t(coef_vals[coef_name, ]), t(ci[coef_name, ]))
    coef_df <- rbind(coef_df, coef_row) # append to the dataframe
  }

  # Append the subject random intercepts.
  ranef_df <- rbind(ranef_df, data.frame(hemroi = hemroi, t(coef(mod)$subj)))
}

# Write coef_df to file if it doesn't already exist
overwrite <- TRUE
outfile <- file.path(data_dir, "tau_roi_suvrs_eoad_long_slopes.csv")
if (overwrite || !file.exists(outfile)) {
  write.csv(coef_df, outfile, row.names = FALSE)
  print(paste("Wrote", outfile))
}

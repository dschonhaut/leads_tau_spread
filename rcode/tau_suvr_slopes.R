# `radian` in Terminal to get R shell
library(lme4)
library(lmerTest)
library(ggplot2)
library(tidyverse)

# --------------------------
# Import and format the data
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
filepath <- file.path(data_dir, "tau_rois.csv")
tau_rois <- read.csv(filepath)

# Replace missing values with NA
cols_to_replace <- c("minority", "apoe_genotype", "apoe4_alleles", "cdrsob_baseline", "mmse_baseline")
tau_rois[cols_to_replace] <- lapply(tau_rois[cols_to_replace], function(x) {
  x[x == ""] <- NA
  x
})

# Subset the dataframe
tau_rois$subjhemroi <- paste(tau_rois$subj, tau_rois$hemroi, sep = "_")
tau_rois$cdrsob_gt4 <- tau_rois$cdrsob_baseline > 4
tau_rois_bline <- tau_rois[(tau_rois$visit == 1), ]
tau_rois_bline_eoad_con <- tau_rois_bline[(tau_rois_bline$dx %in% c("EOAD", "CON")), ]
tau_rois_bline_eoad_eonad <- tau_rois_bline[(tau_rois_bline$dx %in% c("EOAD", "EOnonAD")), ]
tau_rois_bline_eonad_con <- tau_rois_bline[(tau_rois_bline$dx %in% c("EOnonAD", "CON")), ]
tau_rois_bline_eoad <- tau_rois_bline[(tau_rois_bline$dx == "EOAD"), ]
tau_rois_eoad <- tau_rois[(tau_rois$dx == "EOAD") & (tau_rois$pet_visits > 1), ]
tau_rois_eoad_all <- tau_rois[(tau_rois$dx == "EOAD"), ]

# Convert strings to factors
tau_rois_bline <- as.data.frame(unclass(tau_rois_bline), stringsAsFactors = TRUE)
tau_rois_bline_eoad_con <- as.data.frame(unclass(tau_rois_bline_eoad_con), stringsAsFactors = TRUE)
tau_rois_bline_eoad_eonad <- as.data.frame(unclass(tau_rois_bline_eoad_eonad), stringsAsFactors = TRUE)
tau_rois_bline_eoad_eonad$dx <- ordered(tau_rois_bline_eoad_eonad$dx, c("EOnonAD", "EOAD"))
tau_rois_bline_eonad_con <- as.data.frame(unclass(tau_rois_bline_eonad_con), stringsAsFactors = TRUE)
tau_rois_bline_eoad <- as.data.frame(unclass(tau_rois_bline_eoad), stringsAsFactors = TRUE)
tau_rois_eoad <- as.data.frame(unclass(tau_rois_eoad), stringsAsFactors = TRUE)
tau_rois_eoad_all <- as.data.frame(unclass(tau_rois_eoad_all), stringsAsFactors = TRUE)

eoad_subjs <- levels(as.factor(tau_rois_eoad$subj))
hemrois <- levels(as.factor(tau_rois_eoad$hemroi))
subjhemrois <- levels(as.factor(tau_rois_eoad$subjhemroi))

## ------------
## EOAD vs. CON
mod_icpt <- lm(
  value ~ 1,
  data <- tau_rois_bline_eoad_con
)
mod_dx <- lm(
  value ~ 1 + dx,
  data <- tau_rois_bline_eoad_con
)
mod_dx_roi <- lm(
  value ~ 1 + dx + hemroi,
  data <- tau_rois_bline_eoad_con
)
mod_dx.roi <- lm(
  value ~ 1 + dx + hemroi + dx:hemroi,
  data <- tau_rois_bline_eoad_con
)

# Test: Tau PET SUVR values vary between diagnostic groups
anova(mod_dx, mod_icpt)
summary(mod_dx)[["coefficients"]]

# Test: Tau PET values vary between diagnostic groups, differentially across regions.
anova(mod_dx.roi, mod_dx_roi)

# Run a separate model in each region and extract the coefficients.
coef_df <- data.frame()
for (hemroi in hemrois) {
  dat <- tau_rois_bline_eoad_con[(tau_rois_bline_eoad_con$hemroi == hemroi), ]
  mod <- lm(value ~ 1 + dx, data = dat)
  ci <- confint(mod, method = "boot") # 95% confidence interval by parametric bootstrapping
  colnames(ci) <- c("ci2.5", "ci97.5")
  coef_vals <- summary(mod)[["coefficients"]]
  coef_names <- rownames(coef_vals)
  for (coef_name in coef_names[2:length(coef_names)]) {
    coef_row <- data.frame(hemroi = hemroi, nobs = nrow(dat), coef = coef_name, t(coef_vals[coef_name, ]), t(ci[coef_name, ]))
    coef_df <- rbind(coef_df, coef_row) # append to the dataframe
  }
}
names(coef_df)[names(coef_df) == "Estimate"] <- "est"
names(coef_df)[names(coef_df) == "Std..Error"] <- "sem"
names(coef_df)[names(coef_df) == "t.value"] <- "tval"
names(coef_df)[names(coef_df) == "Pr...t.."] <- "pval"
coef_df <- relocate(coef_df, any_of(c("ci2.5", "ci97.5")), .after = 5)

# FDR-correct the p-values
coef_df$pval_fdr <- p.adjust(coef_df$pval, method = "fdr")

# Sort coef_df by p-value descending
coef_df <- coef_df[order(coef_df$pval, decreasing = TRUE), ]

# --------------------------------------
# MEMs for longitudinal analysis in EOAD
coef_df <- data.frame()
for (hemroi in hemrois) {
  dat <- tau_rois_eoad[(tau_rois_eoad$hemroi == hemroi), ]
  mod <- lmer(value ~ 1 + years_from_baseline_pet + (1 + years_from_baseline_pet | subj), data = dat, REML = TRUE)
  subj_coefs <- coef(mod)$subj
  fefs <- fixef(mod)

  # Change subj_coefs column names
  colnames(subj_coefs) <- c("icpt", "years")

  # Make rownames into the first column, called subj
  subj_coefs$subj <- rownames(subj_coefs)

  # Add hemroi as the second column
  subj_coefs$hemroi <- hemroi

  # Reorder columns like (subj, hemroi, icpt, years)
  subj_coefs <- subj_coefs[, c(3, 4, 1, 2)]

  # Reset rownames to NULL
  rownames(subj_coefs) <- NULL

  subj_coefs$icpt <- subj_coefs$icpt + fefs[1]
  subj_coefs$years <- subj_coefs$years + fefs[2]

  # Append subj_coefs to coef_df
  coef_df <- rbind(coef_df, subj_coefs)
}

# Print coef_df dimensions
dim(coef_df)

# Separate hemroi column into hem and roi columns
coef_df$hem <- substr(coef_df$hemroi, 1, 1)
coef_df$roi <- substr(coef_df$hemroi, 3, nchar(coef_df$hemroi))

# Save these columns as factors
coef_df$subj <- as.factor(coef_df$subj)
coef_df$hem <- as.factor(coef_df$hem)
coef_df$roi <- as.factor(coef_df$roi)
coef_df$hemroi <- as.factor(coef_df$hemroi)

# Reorder columns like (subj, hemroi, hem, roi, icpt, years)
coef_df <- coef_df[, c(1, 2, 5, 6, 3, 4)]

# Write coef_df to file if it doesn't already exist
overwrite = TRUE
outfile <- file.path(data_dir, "tau_rois_eoad_longitudinal.csv")
if (overwrite || !file.exists(outfile)) {
  write.csv(coef_df, outfile, row.names = FALSE)
  print(paste("Wrote", outfile))
}

# Create a facet graph with one scatterplot for each of the 36 ROIs in coef_df$roi
# For each scatterplot, x=icpt, y=years, and color=hem.
dat = coef_df

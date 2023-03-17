library(tidyverse)
library(mgcv)
library(gamm4)
library(gratia)
library(viridis)

# ------------------------------------------------------------------------------
# Import and format the data
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
filepath <- file.path(data_dir, "tau_rois_spec.csv")
tau_rois <- read.csv(filepath)
tau_rois <- subset(tau_rois, select = -X)
filepath <- file.path(data_dir, "tau_rois_spec_eoad_long2.csv")
tau_rois_long <- read.csv(filepath)

# Replace missing values with NA
cols_to_replace <- c("minority", "apoe_genotype", "apoe4_alleles", "apoe2_alleles", "cdrsob_baseline", "mmse_baseline")
for (col in cols_to_replace) {
  tau_rois[tau_rois[[col]] == "", col] <- NA
  tau_rois[tau_rois[[col]] == -999, col] <- NA
  tau_rois_long[tau_rois_long[[col]] == "", col] <- NA
  tau_rois_long[tau_rois_long[[col]] == -999, col] <- NA
}

# Subset the dataframe
tau_rois_eoad <- tau_rois[(tau_rois$dx == "EOAD") & (tau_rois$pet_visits > 1), ]

# Convert strings to factors
tau_rois_eoad <- as.data.frame(unclass(tau_rois_eoad), stringsAsFactors = TRUE)
tau_rois_long <- as.data.frame(unclass(tau_rois_long), stringsAsFactors = TRUE)
tau_rois_long$apoe4_allelesf <- factor(tau_rois_long$apoe4_alleles,
                                       # levels = c(0, 1, 2),
                                       # labels = c("0", "1", "2"),
                                       ordered = TRUE)
tau_rois_long$cdrsob_baselinef <- factor(tau_rois_long$cdrsob_baseline,
                                         ordered = TRUE)
tau_rois_long <- tau_rois_long %>%
  rename(sex = sex_m) %>%
  mutate(sex = ifelse(sex == 0, "f", "m"),
         sex = factor(sex, levels = c("f", "m")))
eoad_subjs <- levels(as.factor(tau_rois_eoad$subj))
rois <- levels(as.factor(tau_rois_eoad$roi))

# ------------------------------------------------------------------------------
mod_roi <- gam(
  delta_suvr_peryear ~
    s(last_suvr, bs = "cs", k = 6) +
    s(roi_full, bs = "re") +
    s(last_suvr, by = roi_full, bs = "cs", k = 6) +
    s(subj, bs = "re"),
  data = tau_rois_long,
  REML = TRUE,
  family = gaussian
)
summary(mod_roi)
gam.check(mod_roi)
plotdat <- plot.gam(mod_roi,
                    select=1,
                    shade=TRUE,
                    shift=fixef(mod_roi)["(Intercept)"],
                    residuals=FALSE,
                    xlim=c(0.95, 4),
                    ylim=c(-0.4, 0.4),
                    # xlab="Tau-SUVR at time(t)",
                    # ylab="ΔTau-SUVR/year: time(t) to time(t+1)",
                    rug=FALSE) +
  abline(h = 0, lty = 1)
plot.gam(mod_roi, pages=4)

# ------------------------------------------------------------------------------
mod_roi <- gam(delta_suvr_peryear ~
                 -1 +
                 roi +
                 s(last_suvr, by=roi_full, bs="cs") +
                 s(subj, bs="re"),
               data=tau_rois_long,
               na.action="na.omit",
               REML=TRUE,
               family=gaussian)
summary(mod_roi)
gam.check(mod_roi)

plot.gam(mod_roi, select=4)

fv <- fitted_values(mod_roi, data=tau_rois_long)
fv |>
  ggplot(aes(x = last_suvr, y = fitted, colour = roi_full)) +
  geom_hline(yintercept = 0) +
  geom_point(data = tau_rois_long, mapping = aes(y = delta_suvr_peryear), size = 0.2) +
  geom_ribbon(aes(x = last_suvr, ymin = lower, ymax = upper, fill = roi_full,
                  colour = NULL),
              alpha = 0.2) +
  geom_line() +
  ylim(-0.5, 0.5) +
  xlab("Baseline tau-PET SUVR") +
  ylab("ΔSUVR/year") +
  guides(colour = FALSE, fill = FALSE) +
  theme_minimal() +
  facet_wrap(~ roi_full)

# ------------------------------------------------------------------------------
keep_cols <- c("delta_suvr_peryear", "roi_full", "last_suvr",
               "cdrsob_baseline", "cdrsob_baselinef", "subj")
mod_cdr <- gam(
  delta_suvr_peryear ~
    s(cdrsob_baseline, bs="cs", k=6) +
    s(roi_full, bs="re") +
    s(cdrsob_baseline, by=roi_full, bs="cs", k=6) +
    s(subj, bs="re"),
  data=tau_rois_long[complete.cases(tau_rois_long[,keep_cols]),],
  REML=TRUE,
  family=gaussian
)
summary(mod_cdr)
gam.check(mod_cdr)
plot.gam(mod_cdr, pages=4)
plotdata <- plot.gam(mod_cdr,
                     select=1,
                     shade=TRUE,
                     xlim=c(0.5, 8),
                     # ylim=c(-0.5, 0.5),
                     xlab="Baseline CDR-SB",
                     ylab="ΔSUVR/year",
                     shift=fixef(mod_roi_full)["(Intercept)"],
                     residuals=TRUE,
                     seWithMean=FALSE)


mod_cdr_base <- gam(
  delta_suvr_peryear ~
    s(cdrsob_baseline, bs="cs", k=6) +
    s(last_suvr, bs="cs", k=6) +
    s(roi_full, bs="re") +
    s(cdrsob_baseline, by=roi_full, bs="cs", k=6) +
    s(last_suvr, by=roi_full, bs="cs", k=6) +
    s(subj, bs="re"),
  data=tau_rois_long[complete.cases(tau_rois_long[,keep_cols]),],
  REML=TRUE,
  family=gaussian
)
summary(mod_cdr_base)
gam.check(mod_cdr_base)
plot.gam(mod_cdr_base, pages=4)
plotdata <- plot.gam(mod_cdr_base,
                     select=1,
                     shade=TRUE,
                     xlim=c(0.5, 8),
                     # ylim=c(-0.5, 0.5),
                     xlab="Baseline CDR-SB",
                     ylab="ΔSUVR/year",
                     shift=fixef(mod_roi_full)["(Intercept)"],
                     residuals=TRUE,
                     seWithMean=FALSE)

# ------------------------------------------------------------------------------
keep_cols <- c("delta_suvr_peryear", "roi_full", "last_suvr",
               "apoe4_alleles", "subj")
mod_apoe <- gam(
  delta_suvr_peryear ~
    s(last_suvr, bs="cs", k=6) +
    s(apoe4_allelesf, bs="re") +
    s(roi_full, bs="re") +
    s(last_suvr, by=apoe4_allelesf, bs="cs", k=6) +
    s(last_suvr, by=roi_full, bs="cs", k=6) +
    s(subj, bs="re"),
  data=tau_rois_long[complete.cases(tau_rois_long[,keep_cols]),],
  REML=TRUE,
  family=gaussian
)
summary(mod_apoe)
gam.check(mod_apoe)
plot.gam(mod_apoe, pages=4)
plotdata <- plot.gam(mod_apoe,
                     select=1,
                     shade=TRUE,
                     xlim=c(0.5, 8),
                     # ylim=c(-0.5, 0.5),
                     xlab="Baseline CDR-SB",
                     ylab="ΔSUVR/year",
                     shift=fixef(mod_roi_full)["(Intercept)"],
                     residuals=TRUE,
                     seWithMean=FALSE)

# ------------------------------------------------------------------------------
keep_cols <- c("delta_suvr_peryear", "roi_full", "last_suvr",
               "age_at_pet", "subj", "visit")
mod_ageatpet <- gam(
  delta_suvr_peryear ~
    s(last_suvr, bs="cs", k=6) +
    s(age_at_pet, bs="cs", k=6) +
    s(roi_full, bs="re") +
    ti(last_suvr, age_at_pet, bs="cs") +
    s(last_suvr, by=roi_full, bs="cs", k=6) +
    s(age_at_pet, by=roi_full, bs="cs", k=6) +
    s(subj, bs="re"),
  data=tau_rois_long[complete.cases(tau_rois_long[,keep_cols]),],
  REML=TRUE,
  family=gaussian
)
summary(mod_ageatpet)
gam.check(mod_ageatpet)
plot.gam(mod_ageatpet, pages=4)
plotdat <- plot.gam(mod_ageatpet,
                    select=31, #1, 2, 23, 31
                    shade=TRUE,
                    shift=fixef(mod_ageatpet)["(Intercept)"],
                    residuals=TRUE,
                    # xlim=c(0.95, 4),
                    ylim=c(-0.5, 0.5),
                    # xlab="Tau-SUVR at time(t)",
                    # ylab="ΔTau-SUVR/year: time(t) to time(t+1)",
                    rug=FALSE) +
  abline(h = 0, lty = 1)

# ------------------------------------------------------------------------------
mod0 <- gamm4(value ~ roi + s(years_from_baseline_pet, by=roi, bs="cr"),
              random=~(1|subj) + (1|subj:roi),
              data=tau_rois_eoad, REML=TRUE, family=gaussian)
mod0 <- gamm4(value ~ roi +
                      s(cdrsob_baseline, by=roi) +
                      s(years_from_baseline_pet, by=roi) +
                      t2(cdrsob_baseline, years_from_baseline_pet, by=roi),
              random=~(1|subj) + (1|subj:roi),
              data=tau_rois_eoad[((tau_rois_eoad$cdrsob_baseline >= 0)),], REML=TRUE, family=gaussian)
summary(mod0$gam)
plot.gam(mod0$gam)
mod1 <- gam(value ~ roi +
                    s(years_from_baseline_pet, by=roi, bs="gp") +
                    s(subj, bs="re"),
            data=tau_rois_eoad, REML=TRUE, family=gaussian,
            drop.intercept=TRUE, select=FALSE)
summary(mod1)
plot(mod1, pages=1)

# ------------------------------------------------------------------------------
keep_cols <- c("delta_suvr_peryear", "roi", "last_suvr", "apoe4_alleles", "apoe2_alleles", "age_at_pet", "cdrsob_baseline", "fbb_centiloids_baseline", "subj")
mod_all <- gam(delta_suvr_peryear ~
                 -1 +
                 roi +
                 apoe2_alleles +
                 apoe4_alleles +
                 s(last_suvr) +
                 s(age_at_pet) +
                 s(cdrsob_baseline) +
                 s(fbb_centiloids_baseline) +
                 s(subj, bs="re"),
               data=tau_rois_long[complete.cases(tau_rois_long[,keep_cols]),],
               na.action="na.omit",
               REML=TRUE,
               family=gaussian)
summary(mod_all)
gam.check(mod_all)
plot.gam(mod_all, pages=1)

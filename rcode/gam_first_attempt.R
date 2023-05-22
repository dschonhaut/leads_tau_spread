library(tidyverse)
library(mgcv)
library(gamm4)
library(gratia)
library(viridis)
library(broom)

# ------------------------------------------------------------------------------
# Import and format the data
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
filepath <- file.path(data_dir, "tau_rois_spec.csv")
tau_rois <- read.csv(filepath)
tau_rois <- subset(tau_rois, select = -X)
filepath <- file.path(data_dir, "tau_rois_spec_eoad_long2.csv")
tau_rois_long <- read.csv(filepath)

# Replace missing values with NA
cols_to_replace <- c("minority", "apoe_genotype", "apoe4_alleles",
                     "apoe2_alleles", "cdrsob_baseline", "mmse_baseline")
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
tau_rois_long$minority <- factor(tau_rois_long$minority, ordered=FALSE)
tau_rois_long$apoe4_allelesf <- factor(tau_rois_long$apoe4_alleles,
                                       # levels = c(0, 1, 2),
                                       # labels = c("0", "1", "2"),
                                       ordered = TRUE)
tau_rois_long$apoe2_allelesf <- factor(tau_rois_long$apoe2_alleles,
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
# MODELS

# Delta tau ~ baseline tau * region
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
gam.check(mod_roi) + abline(0, 1, col="red")
appraise(mod_roi, method="simulate")
plotdat <- plot.gam(mod_roi,
                    select=1,
                    shade=TRUE,
                    shift=fixef(mod_roi)["(Intercept)"],
                    residuals=FALSE,
                    seWithMean=TRUE,
                    xlim=c(0.95, 4),
                    ylim=c(-0.4, 0.4),
                    # xlab="Tau-SUVR at time(t)",
                    # ylab="ΔTau-SUVR/year: time(t) to time(t+1)",
                    rug=FALSE) +
  abline(h = 0, lty = 1)

plot.gam(mod_roi, seWithMean=TRUE, pages=4)

# ------------------------------------------------------------------------------
# Delta tau ~ baseline tau * region
# mod_roi_alt <- gam(
#   delta_suvr_from_baseline ~
#     # s(years_from_baseline_pet, bs = "cs") +
#     # s(baseline_suvr, bs = "cs") +
#     # ti(years_from_baseline_pet, baseline_suvr, bs = "cs") +
#     te(years_from_baseline_pet, baseline_suvr, bs = "cs") +
#     # s(years_from_baseline_pet, by = roi_full, bs = "cs", k = 6) +
#     # s(baseline_suvr, by = roi_full, bs = "cs", k = 6) +
#     # ti(years_from_baseline_pet, baseline_suvr, by = roi_full, bs = "cs", k=12) +
#     # te(years_from_baseline_pet, baseline_suvr, by = roi_full, bs = "cs") +
#     # te(years_from_baseline_pet, baseline_suvr, roi_full, bs = "fs") +
#     s(roi_full, bs = "re") +
#     s(subj, bs = "re"),
#   data = tau_rois_long,
#   REML = TRUE,
#   family = gaussian
# )
mod_roi_alt <- gam(
  delta_suvr_from_baseline ~
    # PARAMETRIC EFFECTS
    # ------------------
    # roi_full +

    # SMOOTHS OF INTEREST
    # -------------------
    s(baseline_suvr, bs = "cs") +
    s(baseline_suvr, by=roi_full, bs = "cs") +
    # ti(years_from_baseline_pet, baseline_suvr, bs = "cs") +
    # te(years_from_baseline_pet, baseline_suvr, bs = "cs") +
    # s(years_from_baseline_pet, by = roi_full, bs = "cs", k = 6) +
    # s(baseline_suvr, by = roi_full, bs = "cs", k = 6) +
    # ti(years_from_baseline_pet, baseline_suvr, by = roi_full, bs = "cs", k=12) +
    # te(years_from_baseline_pet, baseline_suvr, by = roi_full, bs = "cs") +

    # NUISANCE VARIABLES
    # ------------------
    s(years_from_baseline_pet, bs = "cs") +
    s(years_from_baseline_pet, by=roi_full, bs = "cs") +
    s(roi_full, bs = "re") +
    s(subj, bs = "re"),
  data = tau_rois_long,
  REML = TRUE,
  family = gaussian,
  select = TRUE
)
summary(mod_roi_alt)
gam.check(mod_roi_alt) + abline(0, 1, col="red")
appraise(mod_roi_alt, method="simulate")

smooths(mod_roi_alt)
smooth_coefs(mod_roi_alt, term="s(baseline_suvr)")
est_grand_icpt <- fixef(mod_roi_alt)
est_roi_icpt <- smooth_estimates(mod_roi_alt, "s(roi_full)") %>% add_confint() %>% arrange(-est)
est_baseline_suvr <- smooth_estimates(mod_roi_alt, "s(baseline_suvr)") %>% add_confint()
est_years_from_baseline <- smooth_estimates(mod_roi_alt, "s(years_from_baseline_pet)") %>% add_confint()
est_roi_x_baseline_suvr <- smooth_estimates(mod_roi_alt, "s(years_from_baseline_pet):roi_fullPrecuneus") %>% add_confint() %>% arrange(-est)

plot.gam(mod_roi_alt,
         residuals=TRUE,
         seWithMean=TRUE,
         shade=TRUE,
         shift=fixef(mod_roi_alt)["(Intercept)"],
         ylim=c(-0.5, 0.4),
         select = 1)+
  abline(h = 0, lty = 1)
draw(mod_roi_alt,
     scales="fixed",
     residuals=TRUE,
     shift=fixef(mod_roi_alt)["(Intercept)"],
     # size=1,
     # coverage=1,
     overall_uncertainty=TRUE,
     select=1) +
  ylim(-0.6, 0.6)
draw()
plotdat <- plot.gam(mod_roi_alt,
                    select=2,
                    shade=TRUE,
                    #shift=fixef(mod_roi_alt)["(Intercept)"],
                    residuals=TRUE,
                    seWithMean=TRUE,
                    # xlim=c(0.95, 4),
                    ylim=c(-0.4, 0.4),
                    # xlab="Tau-SUVR at time(t)",
                    # ylab="ΔTau-SUVR/year: time(t) to time(t+1)",
                    rug=FALSE) +
  abline(h = 0, lty = 1)

# ------------------------------------------------------------------------------
mod_roi2 <- gam(delta_suvr_peryear ~
                roi_full +
                s(last_suvr, bs="cs") +
                s(last_suvr, by=roi_full, bs="cs"),
                  # s(subj, bs="re"),
               data=tau_rois_long,
               na.action="na.omit",
               REML=TRUE,
               family=gaussian)
summary(mod_roi2)
gam.check(mod_roi2)
tmod_roi2 <- tidy(mod_roi2, parametric=FALSE, conf.int=TRUE)
sm <- mod_roi2$smooth[[1]]
plotdat <- plot.gam(mod_roi2,
                    select=1,
                    shade=TRUE,
                    shift=fixef(mod_roi2)["(Intercept)"],
                    residuals=FALSE,
                    seWithMean=TRUE,
                    # xlim=c(0.95, 4),
                    # ylim=c(-0.4, 0.4),
                    # xlab="Tau-SUVR at time(t)",
                    # ylab="ΔTau-SUVR/year: time(t) to time(t+1)",
                    rug=FALSE) +
  abline(h = 0, lty = 1)

plot.gam(mod_roi2, seWithMean=TRUE)

fv <- fitted_values(mod_roi2, data=tau_rois_long)
fv |>
  ggplot(aes(x = last_suvr, y = fitted, colour = roi_full)) +
  geom_hline(yintercept = 0) +
  geom_point(data = tau_rois_long, mapping = aes(y = delta_suvr_peryear),
             size = 0.2) +
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
# Delta tau ~ baseline tau * CDR
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
plot.gam(mod_cdr, seWithMean=TRUE, pages=4)
plotdata <- plot.gam(mod_cdr,
                     select=1,
                     shade=TRUE,
                     xlim=c(0.5, 8),
                     # ylim=c(-0.5, 0.5),
                     xlab="Baseline CDR-SB",
                     ylab="ΔSUVR/year",
                     shift=fixef(mod_roi_full)["(Intercept)"],
                     residuals=TRUE,
                     seWithMean=TRUE)


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
plot.gam(mod_cdr_base, seWithMean=TRUE, pages=4)
plotdata <- plot.gam(mod_cdr_base,
                     select=1,
                     shade=TRUE,
                     xlim=c(0.5, 8),
                     # ylim=c(-0.5, 0.5),
                     xlab="Baseline CDR-SB",
                     ylab="ΔSUVR/year",
                     shift=fixef(mod_roi_full)["(Intercept)"],
                     residuals=TRUE,
                     seWithMean=TRUE)

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
plot.gam(mod_apoe, seWithMean=TRUE, pages=4)
plotdata <- plot.gam(mod_apoe,
                     select=1,
                     shade=TRUE,
                     xlim=c(0.5, 8),
                     # ylim=c(-0.5, 0.5),
                     xlab="Baseline CDR-SB",
                     ylab="ΔSUVR/year",
                     shift=fixef(mod_roi_full)["(Intercept)"],
                     residuals=TRUE,
                     seWithMean=TRUE)

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
plot.gam(mod_ageatpet, seWithMean=TRUE, pages=4)
plotdat <- plot.gam(mod_ageatpet,
                    select=31, #1, 2, 23, 31
                    shade=TRUE,
                    shift=fixef(mod_ageatpet)["(Intercept)"],
                    residuals=TRUE,
                    seWithMean=TRUE,
                    # xlim=c(0.95, 4),
                    ylim=c(-0.5, 0.5),
                    # xlab="Tau-SUVR at time(t)",
                    # ylab="ΔTau-SUVR/year: time(t) to time(t+1)",
                    rug=FALSE) +
  abline(h = 0, lty = 1)

# ------------------------------------------------------------------------------
# GAMM4
mod0 <- gamm4(value ~ roi + s(years_from_baseline_pet, by=roi, bs="cr"),
              random=~(1|subj) + (1|subj:roi),
              data=tau_rois_eoad, REML=TRUE, family=gaussian)
mod0 <- gamm4(value ~ roi +
                      s(cdrsob_baseline, by=roi) +
                      s(years_from_baseline_pet, by=roi) +
                      t2(cdrsob_baseline, years_from_baseline_pet, by=roi),
              random=~(1|subj) + (1|subj:roi),
              data=tau_rois_eoad[((tau_rois_eoad$cdrsob_baseline >= 0)),],
              REML=TRUE, family=gaussian)
summary(mod0$gam)
plot.gam(mod0$gam, seWithMean=TRUE)
mod1 <- gam(value ~ roi +
                    s(years_from_baseline_pet, by=roi, bs="gp") +
                    s(subj, bs="re"),
            data=tau_rois_eoad, REML=TRUE, family=gaussian,
            drop.intercept=TRUE, select=FALSE)
summary(mod1)
plot(mod1, pages=1)

# ------------------------------------------------------------------------------
# Full model
keep_cols <- c("delta_suvr_peryear", "roi", "last_suvr", "apoe4_alleles",
               "apoe2_alleles", "age_at_pet", "cdrsob_baseline",
               "fbb_centiloids_baseline", "subj")
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
plot.gam(mod_all, seWithMean=TRUE, pages=1)

mod_full <- gam(
  delta_suvr_from_baseline ~
    # PARAMETRIC COEFS
    # ----------------

  # SMOOTHS OF INTEREST
  # -------------------
  s(baseline_suvr, bs = "cs") +
    s(baseline_suvr, roi_full, bs = "fs") +
    s(age_at_baseline, bs = "cs") +
    ti(baseline_suvr, age_at_baseline, bs = "cs") +
    s(baseline_suvr, sex, bs = "fs") +
    s(cdrsob_baseline, bs = "cs") +
    ti(baseline_suvr, cdrsob_baseline, bs = "cs") +
    s(baseline_suvr, apoe4_allelesf, bs = "fs") +
    s(fbb_centiloids_baseline, bs = "cs") +
    ti(baseline_suvr, fbb_centiloids_baseline, bs = "cs") +

    # NUISANCE VARIABLES
    # ------------------
  # s(roi_full, bs = "re") +
  s(years_from_baseline_pet, bs = "cs") +
    s(years_from_baseline_pet, roi_full, bs = "fs") +
    s(subj, bs = "re"),

  data = tau_rois_long,
  REML = TRUE,
  family = gaussian,
  select = TRUE
)
summary(mod_full)
plot.gam(mod_full, seWithMean=TRUE, residuals=TRUE)
draw(mod_full, scales = "fixed", overall_uncertainty = TRUE)

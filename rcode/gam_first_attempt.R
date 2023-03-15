library(tidyverse)
library(mgcv)
library(gamm4)
library(gratia)

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
  tau_rois[tau_rois[[col]] == "" | tau_rois[[col]] == -999, col] <- NA
  tau_rois_long[tau_rois_long[[col]] == "" | tau_rois_long[[col]] == -999, col] <- NA
}

# Subset the dataframe
tau_rois_eoad <- tau_rois[(tau_rois$dx == "EOAD") & (tau_rois$pet_visits > 1), ]

# Convert strings to factors
tau_rois_eoad <- as.data.frame(unclass(tau_rois_eoad), stringsAsFactors = TRUE)
tau_rois_long <- as.data.frame(unclass(tau_rois_long), stringsAsFactors = TRUE)
eoad_subjs <- levels(as.factor(tau_rois_eoad$subj))
rois <- levels(as.factor(tau_rois_eoad$roi))

# ------------------------------------------------------------------------------

mod_roi <- gam(delta_suvr_peryear ~
                 -1 +
                 roi +
                 s(last_suvr, by=roi) +
                 s(subj, bs="re"),
               data=tau_rois_long,
               na.action="na.omit",
               REML=TRUE,
               family=gaussian)
summary(mod_roi)
gam.check(mod_roi)
plot.gam(mod_roi, select=4)
sm <- smooth_estimates(mod_roi, data=tau_rois_long) %>% add_confint()
ggplot(sm, aes(x=last_suvr, y=est, colour=roi)) +
  geom_ribbon(aes(ymin=lower_ci, ymax=upper_ci, colour=NULL, fill=roi), alpha=0.2) +
  geom_line() +
  facet_wrap(~ roi)

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

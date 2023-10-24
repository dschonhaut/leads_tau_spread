library(tictoc)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(patchwork)
library(broom)
library(performance)
library(sjstats)
library(gtsummary)
library(emmeans)
library(dlookr)
library(boot)
library(lme4)
library(lmerTest)
library(mgcv)
library(report)
library(gratia)
library(ggplot2)
library(sjPlot)
library(see)
library(viridis)
library(esquisse)
library(patchwork)

# Define default parameters.
data_dir <- "/Users/dschonhaut/Box/projects/leads_tau_spread/data/ssheets"
save_dir <- "/Users/dschonhaut/box/projects/leads_tau_spread/analysis/r_models/gam/idv_rois"
today <- format(Sys.Date(), "%Y-%m-%d")
overwrite = FALSE
verbose = TRUE

# Load dataframes.
tau_suvrsf <- file.path(data_dir, "tau-rois-agg_eoad-long_formatted_2023-10-21.csv")
tau_suvrs <- read_csv(tau_suvrsf)
tau_suvrs_widef <- file.path(data_dir, "tau-rois-agg_eoad-wide_formatted_2023-10-23.csv")
tau_suvrs_wide <- read_csv(tau_suvrs_widef)
tau_suvrs_wide$parc <- "metarois"
tau_mmsf <- file.path(data_dir, "tau-rois-mms_all-eoad-gt1ftp_2023-10-21.csv")
tau_mms <- read_csv(tau_mmsf)

# Clean dataframes.
transform_df <- function(df) {
  df %>%
    filter(parc %in% c("metarois")) %>%
    mutate(across(where(is.character), as.factor)) %>%
    mutate(apoe4_alleles = factor(apoe4_alleles, ordered = TRUE)) %>%
    { contrasts(.$sex) <- 'contr.sum'; . }
}
tau_suvrs <- transform_df(tau_suvrs)
tau_mms <- transform_df(tau_mms)
tau_suvrs_wide <- transform_df(tau_suvrs_wide)

# Merge the tau_suvrs (1 row per subj/visit/roi) dataframe with the mixed model fits.
right_cols <- c("subj", "parc", "roi", "fixed_icpt", "fixed_time",
                "rdm_icpt", "rdm_time", "icpt", "time")
if (!all(c("icpt", "time") %in% colnames(tau_suvrs))) {
  tau_suvrs <- tau_suvrs %>%
    left_join(tau_mms %>% select(all_of(right_cols)),
              by = c("subj", "parc", "roi"))
}
if (!all(c("icpt", "time") %in% colnames(tau_suvrs))) {
  tau_suvrs_wide <- tau_suvrs_wide %>%
    left_join(tau_suvrs_wide %>% select(all_of(right_cols)),
              by = c("subj", "parc", "roi"))
}

tau_suvrs2 <- tau_suvrs %>%
  filter(visit > 1) %>%
  mutate(suvr_last = suvr - suvr_chg_from_last) %>%
  group_by(subj, roi) %>%
  summarize(
    ftp_visits = unique(ftp_visits),
    suvr_bl = unique(suvr_bl),
    suvr_last = mean(suvr_last),
    suvr_annchg_from_last = mean(suvr_annchg_from_last),
    icpt = unique(icpt),
    time = unique(time),
    age_at_ftp_bl_mc = unique(age_at_ftp_bl_mc),
    fbb_cl_bl_mc = unique(fbb_cl_bl_mc),
    cdr_sb_bl_mc = unique(cdr_sb_bl_mc),
    apoe4_alleles_mc = unique(apoe4_alleles_mc),
    sex = unique(sex)
  ) %>%
  ungroup()

df <- as.tibble(tau_suvrs_wide)
mod_name <- paste("mod_long_full_reml", roi, sep="_")
gmod <- gam(
  suvr_annchg_from_last_temporal ~
    # Main effects
    s(suvr_bl_temporal),
    # s(age_at_ftp_bl_mc) +
    # s(fbb_cl_bl_mc) +
    # s(cdr_sb_bl_mc) +
    # apoe4_alleles +
    # sex,
  data=df,
  method="REML"
)
bs = "cs"
gmod_t_by_t <- gam(suvr_annchg_from_last_temporal ~ s(suvr_bl_temporal, bs = bs), data=df, method="REML")
gmod_p_by_p <- gam(suvr_annchg_from_last_parietal ~ s(suvr_bl_parietal, bs = bs), data=df, method="REML")
gmod_f_by_f <- gam(suvr_annchg_from_last_frontal ~ s(suvr_bl_frontal, bs = bs), data=df, method="REML")

gmod_t_by_p <- gam(suvr_annchg_from_last_temporal ~ s(suvr_bl_parietal, bs = bs), data=df, method="REML")
gmod_t_by_f <- gam(suvr_annchg_from_last_temporal ~ s(suvr_bl_frontal, bs = bs), data=df, method="REML")
gmod_p_by_t <- gam(suvr_annchg_from_last_parietal ~ s(suvr_bl_temporal, bs = bs), data=df, method="REML")
gmod_p_by_f <- gam(suvr_annchg_from_last_parietal ~ s(suvr_bl_frontal, bs = bs), data=df, method="REML")
gmod_f_by_t <- gam(suvr_annchg_from_last_frontal ~ s(suvr_bl_temporal, bs = bs), data=df, method="REML")
gmod_f_by_p <- gam(suvr_annchg_from_last_frontal ~ s(suvr_bl_parietal, bs = bs), data=df, method="REML")

gmod_t_by_tp <- gam(suvr_annchg_from_last_temporal ~ s(suvr_bl_temporal, bs = bs) + s(suvr_bl_parietal, bs = bs), data=df, method="REML")
gmod_t_by_tf <- gam(suvr_annchg_from_last_temporal ~ s(suvr_bl_temporal, bs = bs) + s(suvr_bl_frontal, bs = bs), data=df, method="REML")
gmod_p_by_pt <- gam(suvr_annchg_from_last_parietal ~ s(suvr_bl_parietal, bs = bs) + s(suvr_bl_temporal, bs = bs), data=df, method="REML")
gmod_p_by_pf <- gam(suvr_annchg_from_last_parietal ~ s(suvr_bl_parietal, bs = bs) + s(suvr_bl_frontal, bs = bs), data=df, method="REML")
gmod_f_by_ft <- gam(suvr_annchg_from_last_frontal ~ s(suvr_bl_frontal, bs = bs) + s(suvr_bl_temporal, bs = bs), data=df, method="REML")
gmod_f_by_fp <- gam(suvr_annchg_from_last_frontal ~ s(suvr_bl_frontal, bs = bs) + s(suvr_bl_parietal, bs = bs), data=df, method="REML")

AIC(
  gmod_t_by_t, gmod_t_by_p, gmod_t_by_f, gmod_t_by_tp, gmod_t_by_tf,
  gmod_p_by_p, gmod_p_by_t, gmod_p_by_f, gmod_p_by_pt, gmod_p_by_pf,
  gmod_f_by_f, gmod_f_by_t, gmod_f_by_p, gmod_f_by_ft, gmod_f_by_fp
)

# Don't run these model comparison tests
# if using select=T or a shrinkage smooth
anova.gam(gmod_t_by_t, gmod_t_by_tp, test = "F") #*
anova.gam(gmod_t_by_t, gmod_t_by_tf, test = "F")
anova.gam(gmod_p_by_p, gmod_p_by_pt, test = "F") #*
anova.gam(gmod_p_by_p, gmod_p_by_pf, test = "F")
anova.gam(gmod_f_by_f, gmod_f_by_ft, test = "F")
anova.gam(gmod_f_by_f, gmod_f_by_fp, test = "F") #*

summary(gmod_t_by_t)
summary(gmod_t_by_tp)
summary(gmod_t_by_tf)
summary(gmod_p_by_p)
summary(gmod_p_by_pt)
summary(gmod_p_by_pf)
summary(gmod_f_by_f)
summary(gmod_f_by_ft)
summary(gmod_f_by_fp)

p1 <- (
  draw(
  gmod_t_by_t,
  constant=fixef(gmod_t_by_t),
  rug=F,
  smooth_col="#2E45B8",
  ci_col="#2E45B8",
  nrow=1,
  ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p2 <- (
  draw(
    gmod_t_by_p,
    constant=fixef(gmod_t_by_p),
    rug=F,
    smooth_col="#3EBCD2",
    ci_col="#3EBCD2",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p3 <- (
  draw(
    gmod_t_by_f,
    constant=fixef(gmod_t_by_f),
    rug=F,
    smooth_col="#FF4983",
    ci_col="#FF4983",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p1 + p2 + p3 + plot_layout(ncol = 3)

p4 <- (
  draw(
    gmod_p_by_t,
    constant=fixef(gmod_p_by_t),
    rug=F,
    smooth_col="#2E45B8",
    ci_col="#2E45B8",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p5 <- (
  draw(
    gmod_p_by_p,
    constant=fixef(gmod_p_by_p),
    rug=F,
    smooth_col="#3EBCD2",
    ci_col="#3EBCD2",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p6 <- (
  draw(
    gmod_p_by_f,
    constant=fixef(gmod_p_by_f),
    rug=F,
    smooth_col="#FF4983",
    ci_col="#FF4983",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p4 + p5 + p6 + plot_layout(ncol = 3)

p7 <- (
  draw(
    gmod_f_by_t,
    constant=fixef(gmod_f_by_t),
    rug=F,
    smooth_col="#2E45B8",
    ci_col="#2E45B8",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p8 <- (
  draw(
    gmod_f_by_p,
    constant=fixef(gmod_f_by_p),
    rug=F,
    smooth_col="#3EBCD2",
    ci_col="#3EBCD2",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p9 <- (
  draw(
    gmod_f_by_f,
    constant=fixef(gmod_f_by_f),
    rug=F,
    smooth_col="#FF4983",
    ci_col="#FF4983",
    nrow=1,
    ncol=1
  ) +
    scale_x_continuous(limits=c(1, 4), labels=seq(1, 4, 1), expand=c(0, 0)) +
    scale_y_continuous(limits=c(-0.25, 0.25), labels=seq(-0.2, 0.2, 0.1), expand=c(0, 0))
)
p7 + p8 + p9 + plot_layout(ncol = 3)

p1 + p5 + p9 + plot_layout(ncol = 3)

# Predict!!
### Plot baseline SUVR ~ change in SUVR: Regional
df <- as.tibble(tau_suvrs_wide)
metarois <- c("temporal", "parietal", "frontal")
bs = "cs"

cutoffs <- tibble()#n = integer(0), xmin = numeric(0), xmax = numeric(0), ymin = numeric(0), ymax = numeric(0))
pred <- tibble()
for (roi in metarois) {
  # Create variable names based on the current ROI
  dv <- str_c("suvr_annchg_from_last_", roi, sep="")
  iv <- str_c("suvr_bl_", roi, sep="")
  smooth_iv <- str_c("s(suvr_bl_", roi, ")", sep="")
  formula_str <- str_c(dv, " ~ ", "s(", iv, ", bs=bs)", sep="")
  formula_obj <- as.formula(formula_str)
  xmod <- gam(formula_obj, data=df, method="REML")
  mod_terms <- predict.gam(xmod, type="terms")
  icpt <- fixef(xmod)["(Intercept)"]

  # Get model residuals
  dat <- as.tibble(df)
  dat$mod <- str_c("gmod_", str_sub(roi, 1, 1), "_by_", str_sub(roi, 1, 1), sep="")
  dat$resid <- residuals.gam(xmod, type = "response")
  dat$partial <- (dat$resid + icpt + mod_terms[,c(smooths(xmod))])

  # Calculate cutoffs for the current ROI
  xmod_names <- paste(names(fixef(xmod)), smooths(xmod))
  cutoffs_roi <- df %>%
    summarize(
      roi = roi,
      n = n(),
      xmin = round(max(min(!!sym(iv)), quantile(!!sym(iv), 0.25) - (1.5 * IQR(!!sym(iv)))), 1),
      xmax = round(min(max(!!sym(iv)), quantile(!!sym(iv), 0.75) + (1.5 * IQR(!!sym(iv)))), 1),
      ymin = round(max(min(!!sym(dv)), quantile(!!sym(dv), 0.25) - (1.5 * IQR(!!sym(dv)))), 2),
      ymax = round(min(max(!!sym(dv)), quantile(!!sym(dv), 0.75) + (1.5 * IQR(!!sym(dv)))), 2)
    )
  xvals <- seq(cutoffs_roi$xmin, cutoffs_roi$xmax, by = 0.01)

  # Append the results to the cutoffs data frame
  cutoffs <- rbind(cutoffs, cutoffs_roi)

  # Get model estimates
  pred_ <- with(df, expand.grid(xvals))
  names(pred_) <- c(iv)
  pred_ <- cbind(
    pred_,
    data.frame(
      predict(
        xmod,
        newdata = pred_,
        type = "response",
        se.fit = TRUE,
        terms = c("(Intercept)", smooth_iv),
        newdata.guaranteed = TRUE
      )
    )
  )
  pred_$roi <- roi
  pred_ <- pred_ %>%
    rename(
      suvr_bl = .data[[iv]]
    ) %>%
    mutate(
      upper = fit + (2 * se.fit),
      lower = fit - (2 * se.fit)
    ) %>%
    tibble()
  pred <- rbind(pred, pred_)
}

# Save the output.
pred_outputf <- file.path(data_dir, "tau-metarois-gam_preds_all-eoad-gt1ftp_2023-10-24.csv")
write_csv(pred, pred_outputf)
print(str_c("Saved", pred_outputf))


# Visualize data
esquisser(tau_suvrs)
esquisser(tau_mms)
esquisser(tau_suvrs2)

tau_mms %>%
 filter(parc %in% "metarois") %>%
 ggplot() +
 aes(x = icpt, y = time, colour = ftp_visits) +
 geom_point(shape = "circle", size = 2L) +
 scale_color_viridis_c(option = "viridis", direction = 1) +
 theme_linedraw() +
 facet_wrap(vars(roi))

View(tau_suvrs)

df <- tau_suvrs
contrasts(tau_suvrs$sex) <- "contr.sum"
roi <- "frontal"
df_roi <- df[df$roi == roi, ]
model <- lmer(suvr ~ ftp_yrs_from_bl + (1 | subj) + (0 + ftp_yrs_from_bl | subj),
              data = df_roi)
